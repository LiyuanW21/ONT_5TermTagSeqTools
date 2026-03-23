#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict, Counter

import pysam


def parse_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
        elif " " in part:
            key, val = part.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def parse_gff(gff_path):
    gene_info = {}
    with open(gff_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attrs = parts
            if feature != "CDS":
                continue
            start_i = int(start)
            end_i = int(end)
            attr_map = parse_attributes(attrs)
            gene_id = (
                attr_map.get("Parent")
                or attr_map.get("gene_id")
                or attr_map.get("ID")
                or attr_map.get("Name")
            )
            if not gene_id:
                continue
            gene_id = gene_id.split(",")[0]
            gene_key = f"{seqid}|{gene_id}|{strand}"
            if gene_key not in gene_info:
                gene_info[gene_key] = {
                    "seqid": seqid,
                    "gene_id": gene_id,
                    "strand": strand,
                    "start": start_i,
                    "end": end_i,
                }
            else:
                gene_info[gene_key]["start"] = min(gene_info[gene_key]["start"], start_i)
                gene_info[gene_key]["end"] = max(gene_info[gene_key]["end"], end_i)
    for info in gene_info.values():
        if info["strand"] == "+":
            info["cds_start"] = info["start"]
        else:
            info["cds_start"] = info["end"]
    return gene_info


def build_bins(gene_info, bin_size):
    bins = defaultdict(lambda: defaultdict(list))
    for gene_key, info in gene_info.items():
        start = info["start"]
        end = info["end"]
        for bin_id in range(start // bin_size, end // bin_size + 1):
            bins[info["seqid"]][bin_id].append(gene_key)
    return bins


def find_best_gene(seqid, aln_start, aln_end, bins, gene_info, bin_size):
    if seqid not in bins:
        return None
    start_bin = aln_start // bin_size
    end_bin = aln_end // bin_size
    candidates = set()
    for bin_id in range(start_bin, end_bin + 1):
        candidates.update(bins[seqid].get(bin_id, []))
    best = None
    best_olap = 0
    for gene_key in candidates:
        info = gene_info[gene_key]
        if aln_end < info["start"] or aln_start > info["end"]:
            continue
        olap = min(aln_end, info["end"]) - max(aln_start, info["start"]) + 1
        if olap > best_olap:
            best_olap = olap
            best = gene_key
    return best


def read_5p(read):
    if read.is_reverse:
        return read.reference_end
    return read.reference_start + 1


def passes_cds_overlap(read, info, min_overlap):
    aln_start = read.reference_start + 1
    aln_end = read.reference_end
    cds = info["cds_start"]
    if not (aln_start <= cds <= aln_end):
        return False
    if info["strand"] == "+":
        window_start = cds
        window_end = cds + min_overlap - 1
        if aln_end < window_end:
            return False
    else:
        window_start = cds - min_overlap + 1
        window_end = cds
        if aln_start > window_start:
            return False
    return True


def nad_tss_keep(read_5p_pos, strand, tss_pos, downstream_window):
    if strand == "+":
        return read_5p_pos <= tss_pos + downstream_window
    return read_5p_pos >= tss_pos - downstream_window


def compute_nad_tss(
    tagged_bam,
    gene_info,
    bins,
    bin_size,
    mapq,
    min_overlap,
    threads,
):
    tss_counts = defaultdict(Counter)
    with pysam.AlignmentFile(tagged_bam, "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < mapq:
                continue
            gene_key = find_best_gene(
                bam.get_reference_name(read.reference_id),
                read.reference_start + 1,
                read.reference_end,
                bins,
                gene_info,
                bin_size,
            )
            if not gene_key:
                continue
            info = gene_info[gene_key]
            if not passes_cds_overlap(read, info, min_overlap):
                continue
            tss_counts[gene_key][read_5p(read)] += 1
    tss_pos = {}
    for gene_key, counter in tss_counts.items():
        info = gene_info[gene_key]
        max_count = max(counter.values())
        candidates = [pos for pos, cnt in counter.items() if cnt == max_count]
        if info["strand"] == "+":
            tss_pos[gene_key] = min(candidates)
        else:
            tss_pos[gene_key] = max(candidates)
    return tss_pos


def filter_bam(
    bam_path,
    out_path,
    gene_info,
    bins,
    bin_size,
    mapq,
    min_overlap,
    nad_tss=None,
    downstream_window=50,
    apply_nad_tss=False,
    threads=1,
):
    stats = defaultdict(int)
    with pysam.AlignmentFile(
        bam_path, "rb", threads=threads
    ) as bam, pysam.AlignmentFile(
        out_path, "wb", template=bam, threads=threads
    ) as out_bam:
        for read in bam.fetch(until_eof=True):
            stats["total"] += 1
            if read.is_unmapped:
                continue
            stats["mapped"] += 1
            if read.is_secondary or read.is_supplementary:
                continue
            stats["primary"] += 1
            if read.mapping_quality < mapq:
                continue
            stats["mapq"] += 1
            gene_key = find_best_gene(
                bam.get_reference_name(read.reference_id),
                read.reference_start + 1,
                read.reference_end,
                bins,
                gene_info,
                bin_size,
            )
            if not gene_key:
                continue
            info = gene_info[gene_key]
            if not passes_cds_overlap(read, info, min_overlap):
                continue
            stats["cds_overlap"] += 1
            if apply_nad_tss and nad_tss is not None and gene_key in nad_tss:
                if not nad_tss_keep(
                    read_5p(read), info["strand"], nad_tss[gene_key], downstream_window
                ):
                    continue
            if apply_nad_tss:
                stats["nad_tss"] += 1
            out_bam.write(read)
            stats["kept"] += 1
    pysam.index(out_path)
    return stats


def main():
    parser = argparse.ArgumentParser(description="Filter complete reads for NAD-tagSeq.")
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--mapq", type=int, default=10)
    parser.add_argument("--cds-overlap", type=int, default=50)
    parser.add_argument("--nad-tss-window", type=int, default=50)
    parser.add_argument("--bin-size", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    gene_info = parse_gff(args.gff)
    bins = build_bins(gene_info, args.bin_size)

    bam_files = [
        f
        for f in os.listdir(args.bam_dir)
        if f.endswith(".sorted.bam") and os.path.isfile(os.path.join(args.bam_dir, f))
    ]

    samples = defaultdict(dict)
    for bam in bam_files:
        if "_tagged" in bam:
            sample = bam.replace("_tagged", "").replace(".sorted.bam", "")
            samples[sample]["tagged"] = bam
        elif "_nontagged" in bam:
            sample = bam.replace("_nontagged", "").replace(".sorted.bam", "")
            samples[sample]["nontagged"] = bam

    stats_rows = []
    for sample, group in sorted(samples.items()):
        tagged_bam = group.get("tagged")
        nontagged_bam = group.get("nontagged")
        if not tagged_bam and not nontagged_bam:
            continue

        nad_tss = {}
        if tagged_bam:
            nad_tss = compute_nad_tss(
                os.path.join(args.bam_dir, tagged_bam),
                gene_info,
                bins,
                args.bin_size,
                args.mapq,
                args.cds_overlap,
                args.threads,
            )

        if tagged_bam:
            out_path = os.path.join(
                args.out_dir, tagged_bam.replace(".sorted.bam", ".complete.bam")
            )
            stats = filter_bam(
                os.path.join(args.bam_dir, tagged_bam),
                out_path,
                gene_info,
                bins,
                args.bin_size,
                args.mapq,
                args.cds_overlap,
                nad_tss=nad_tss,
                downstream_window=args.nad_tss_window,
                apply_nad_tss=False,
                threads=args.threads,
            )
            stats_rows.append(
                {
                    "sample": sample,
                    "group": "tagged",
                    **stats,
                }
            )

        if nontagged_bam:
            out_path = os.path.join(
                args.out_dir, nontagged_bam.replace(".sorted.bam", ".complete.bam")
            )
            stats = filter_bam(
                os.path.join(args.bam_dir, nontagged_bam),
                out_path,
                gene_info,
                bins,
                args.bin_size,
                args.mapq,
                args.cds_overlap,
                nad_tss=nad_tss,
                downstream_window=args.nad_tss_window,
                apply_nad_tss=True,
                threads=args.threads,
            )
            stats_rows.append(
                {
                    "sample": sample,
                    "group": "nontagged",
                    **stats,
                }
            )

    stats_path = os.path.join(args.out_dir, "complete_reads_stats.tsv")
    with open(stats_path, "w", encoding="utf-8") as handle:
        header = [
            "sample",
            "group",
            "total",
            "mapped",
            "primary",
            "mapq",
            "cds_overlap",
            "nad_tss",
            "kept",
        ]
        handle.write("\t".join(header) + "\n")
        for row in stats_rows:
            handle.write(
                "\t".join(str(row.get(col, 0)) for col in header) + "\n"
            )


if __name__ == "__main__":
    main()
