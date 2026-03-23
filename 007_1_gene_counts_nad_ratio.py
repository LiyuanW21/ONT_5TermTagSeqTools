#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def ensure_bam_index(bam_path: Path) -> None:
    bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
    if not bai_path.exists():
        raise FileNotFoundError(f"Missing BAM index: {bai_path}")


def parse_bam_name(bam_path: Path) -> Tuple[str, str]:
    name = bam_path.name
    if not name.endswith(".complete.bam"):
        raise ValueError(f"Unexpected BAM name: {name}")
    core = name[: -len(".complete.bam")]
    for suffix in ("_tagged", "_nontagged"):
        if core.endswith(suffix):
            return core[: -len(suffix)], suffix[1:]
    raise ValueError(f"Missing tagged/nontagged suffix in {name}")


def replicate_key(sample_base: str, regex_pattern: str) -> str:
    match = re.match(regex_pattern, sample_base)
    return match.group(1) if match else sample_base


def run_bedtools_counts(bam_path: Path, bed_path: Path, out_path: Path) -> None:
    cmd = ["bedtools", "coverage", "-a", str(bed_path), "-b", str(bam_path), "-counts"]
    with out_path.open("w", newline="") as out_handle:
        subprocess.run(cmd, check=True, stdout=out_handle)


def read_counts(counts_path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    with counts_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 7:
                continue
            gene_id = row[3]
            try:
                count = int(row[6])
            except ValueError:
                count = 0
            counts[gene_id] = count
    return counts


def load_gene_ids(bed_path: Path) -> List[str]:
    gene_ids: List[str] = []
    with bed_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 4:
                gene_ids.append(row[3])
    return gene_ids


def write_sample_counts(
    out_path: Path,
    gene_ids: Iterable[str],
    tagged: Dict[str, int],
    nontagged: Dict[str, int],
) -> None:
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "tagged_reads", "nontagged_reads", "total_reads", "tagging_ratio"])
        for gene_id in gene_ids:
            tagged_count = tagged.get(gene_id, 0)
            nontagged_count = nontagged.get(gene_id, 0)
            total = tagged_count + nontagged_count
            ratio = tagged_count / total if total else 0.0
            writer.writerow([gene_id, tagged_count, nontagged_count, total, f"{ratio:.6f}"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Gene-level tagged/nontagged counts and tagging ratio")
    parser.add_argument("--bam-dir", required=True, help="Directory containing *.complete.bam files")
    parser.add_argument("--bed", required=True, help="Gene BED file for counting")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--replicate-keys",
        nargs="+",
        default=None,
        help="Ordered replicate keys for fixed-column outputs, e.g. rep1 rep2 rep3",
    )
    parser.add_argument(
        "--replicate-regex",
        default=r"^(.+?)(?:_pass.*)?$",
        help="Regex with one capturing group to extract replicate key from sample base",
    )
    parser.add_argument("--min-tagged-reads", type=int, default=2, help="Minimum tagged reads")
    parser.add_argument("--min-tagging-ratio", type=float, default=0.002, help="Minimum tagging ratio")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum replicate passes")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    bam_dir = Path(args.bam_dir)
    bed_path = Path(args.bed)
    output_dir = Path(args.output_dir)

    if not bam_dir.exists():
        raise FileNotFoundError(f"BAM directory not found: {bam_dir}")
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    counts_dir = output_dir / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)

    gene_ids = load_gene_ids(bed_path)
    bam_paths = sorted(bam_dir.glob("*.complete.bam"))
    if not bam_paths:
        raise FileNotFoundError(f"No BAMs found in {bam_dir}")

    counts_index: Dict[Tuple[str, str], Dict[str, int]] = {}

    for bam_path in bam_paths:
        ensure_bam_index(bam_path)
        sample_base, group = parse_bam_name(bam_path)
        out_path = counts_dir / f"{sample_base}_{group}.counts.tsv"
        if not out_path.exists():
            run_bedtools_counts(bam_path, bed_path, out_path)
        counts_index[(sample_base, group)] = read_counts(out_path)

    sample_bases = sorted({key[0] for key in counts_index.keys()})
    rep_to_sample: Dict[str, str] = {}
    for sample_base in sample_bases:
        rep = replicate_key(sample_base, args.replicate_regex)
        if rep not in rep_to_sample:
            rep_to_sample[rep] = sample_base

    all_counts_path = output_dir / "gene_counts_all_samples.tsv"
    with all_counts_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header: List[str] = ["gene_id"]
        for sample_base in sample_bases:
            sample_label = replicate_key(sample_base, args.replicate_regex)
            header.extend(
                [
                    f"{sample_label}_tagged_reads",
                    f"{sample_label}_nontagged_reads",
                    f"{sample_label}_total_reads",
                    f"{sample_label}_tagging_ratio",
                ]
            )
        writer.writerow(header)
        for gene_id in gene_ids:
            row: List[str] = [gene_id]
            for sample_base in sample_bases:
                tagged = counts_index.get((sample_base, "tagged"), {})
                nontagged = counts_index.get((sample_base, "nontagged"), {})
                tagged_count = tagged.get(gene_id, 0)
                nontagged_count = nontagged.get(gene_id, 0)
                total = tagged_count + nontagged_count
                ratio = tagged_count / total if total else 0.0
                row.extend([str(tagged_count), str(nontagged_count), str(total), f"{ratio:.6f}"])
            writer.writerow(row)

    per_sample_dir = output_dir / "per_sample"
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    for sample_base in sample_bases:
        tagged = counts_index.get((sample_base, "tagged"), {})
        nontagged = counts_index.get((sample_base, "nontagged"), {})
        out_path = per_sample_dir / f"{sample_base}.gene_counts.tsv"
        write_sample_counts(out_path, gene_ids, tagged, nontagged)

    filtered_path = output_dir / "tagging_genes_by_sample.tsv"
    with filtered_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "replicate_key",
                "gene_id",
                "tagged_reads",
                "nontagged_reads",
                "total_reads",
                "tagging_ratio",
            ]
        )
        for sample_base in sample_bases:
            sample_label = replicate_key(sample_base, args.replicate_regex)
            tagged = counts_index.get((sample_base, "tagged"), {})
            nontagged = counts_index.get((sample_base, "nontagged"), {})
            rep_key = replicate_key(sample_base, args.replicate_regex)
            for gene_id in gene_ids:
                tagged_count = tagged.get(gene_id, 0)
                nontagged_count = nontagged.get(gene_id, 0)
                total = tagged_count + nontagged_count
                ratio = tagged_count / total if total else 0.0
                if tagged_count >= args.min_tagged_reads and ratio > args.min_tagging_ratio:
                    writer.writerow(
                        [
                            sample_label,
                            rep_key,
                            gene_id,
                            tagged_count,
                            nontagged_count,
                            total,
                            f"{ratio:.6f}",
                        ]
                    )

    by_sample_dir = output_dir / "tagging_by_sample"
    by_sample_dir.mkdir(parents=True, exist_ok=True)
    for sample_base in sample_bases:
        sample_label = replicate_key(sample_base, args.replicate_regex)
        tagged = counts_index.get((sample_base, "tagged"), {})
        nontagged = counts_index.get((sample_base, "nontagged"), {})
        out_path = by_sample_dir / f"{sample_label}.tagging.tsv"
        with out_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene_id", "tagged_reads", "nontagged_reads", "total_reads", "tagging_ratio"])
            for gene_id in gene_ids:
                tagged_count = tagged.get(gene_id, 0)
                nontagged_count = nontagged.get(gene_id, 0)
                total = tagged_count + nontagged_count
                ratio = tagged_count / total if total else 0.0
                if tagged_count >= args.min_tagged_reads and ratio > args.min_tagging_ratio:
                    writer.writerow([gene_id, tagged_count, nontagged_count, total, f"{ratio:.6f}"])

    keys = args.replicate_keys
    if keys is None:
        keys = sorted(rep_to_sample.keys())

    replicate_summary: Dict[str, List[float]] = {gid: [] for gid in gene_ids}
    replicate_tagged: Dict[str, List[int]] = {gid: [] for gid in gene_ids}
    replicate_total: Dict[str, List[int]] = {gid: [] for gid in gene_ids}
    replicate_passes: Dict[str, int] = {gid: 0 for gid in gene_ids}

    for sample_base in sample_bases:
        rep_key = replicate_key(sample_base, args.replicate_regex)
        if rep_key not in keys:
            continue
        tagged = counts_index.get((sample_base, "tagged"), {})
        nontagged = counts_index.get((sample_base, "nontagged"), {})
        for gene_id in gene_ids:
            tagged_count = tagged.get(gene_id, 0)
            nontagged_count = nontagged.get(gene_id, 0)
            total = tagged_count + nontagged_count
            ratio = tagged_count / total if total else 0.0
            replicate_summary[gene_id].append(ratio)
            replicate_tagged[gene_id].append(tagged_count)
            replicate_total[gene_id].append(total)
            if tagged_count >= args.min_tagged_reads and ratio > args.min_tagging_ratio:
                replicate_passes[gene_id] += 1

    replicate_out = output_dir / "tagging_genes_replicates.tsv"
    with replicate_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        dynamic_cols = []
        for rep in keys:
            dynamic_cols.extend([f"{rep}_ratio", f"{rep}_tagged_reads", f"{rep}_total_reads"])
        writer.writerow(
            [
                "gene_id",
                "replicate_passes",
                "mean_ratio",
                "mean_tagged_reads",
                "mean_total_reads",
                *dynamic_cols,
            ]
        )
        for gene_id in gene_ids:
            ratios = replicate_summary.get(gene_id, [])
            mean_ratio = sum(ratios) / len(ratios) if ratios else 0.0
            tagged_vals = replicate_tagged.get(gene_id, [])
            total_vals = replicate_total.get(gene_id, [])
            mean_tagged = sum(tagged_vals) / len(tagged_vals) if tagged_vals else 0.0
            mean_total = sum(total_vals) / len(total_vals) if total_vals else 0.0
            if replicate_passes[gene_id] >= args.min_replicates:
                per_sample = []
                for rep in keys:
                    sample_base = rep_to_sample.get(rep)
                    if not sample_base:
                        per_sample.extend(["0.000000", "0", "0"])
                        continue
                    tagged = counts_index.get((sample_base, "tagged"), {})
                    nontagged = counts_index.get((sample_base, "nontagged"), {})
                    tagged_count = tagged.get(gene_id, 0)
                    nontagged_count = nontagged.get(gene_id, 0)
                    total = tagged_count + nontagged_count
                    ratio = tagged_count / total if total else 0.0
                    per_sample.extend([f"{ratio:.6f}", str(tagged_count), str(total)])
                writer.writerow(
                    [
                        gene_id,
                        replicate_passes[gene_id],
                        f"{mean_ratio:.6f}",
                        f"{mean_tagged:.2f}",
                        f"{mean_total:.2f}",
                        *per_sample,
                    ]
                )

    top10_out = output_dir / "tagging_ratio_top10_replicates.tsv"
    with top10_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        dynamic_cols = []
        for rep in keys:
            dynamic_cols.extend([f"{rep}_ratio", f"{rep}_tagged_reads", f"{rep}_total_reads"])
        writer.writerow(
            [
                "gene_id",
                "replicate_passes",
                "mean_ratio",
                "mean_tagged_reads",
                "mean_total_reads",
                *dynamic_cols,
            ]
        )
        rows: List[Tuple[str, int, float, float, float, List[str]]] = []
        for gene_id in gene_ids:
            ratios = replicate_summary.get(gene_id, [])
            mean_ratio = sum(ratios) / len(ratios) if ratios else 0.0
            tagged_vals = replicate_tagged.get(gene_id, [])
            total_vals = replicate_total.get(gene_id, [])
            mean_tagged = sum(tagged_vals) / len(tagged_vals) if tagged_vals else 0.0
            mean_total = sum(total_vals) / len(total_vals) if total_vals else 0.0
            if replicate_passes[gene_id] >= args.min_replicates:
                per_sample = []
                for rep in keys:
                    sample_base = rep_to_sample.get(rep)
                    if not sample_base:
                        per_sample.extend(["0.000000", "0", "0"])
                        continue
                    tagged = counts_index.get((sample_base, "tagged"), {})
                    nontagged = counts_index.get((sample_base, "nontagged"), {})
                    tagged_count = tagged.get(gene_id, 0)
                    nontagged_count = nontagged.get(gene_id, 0)
                    total = tagged_count + nontagged_count
                    ratio = tagged_count / total if total else 0.0
                    per_sample.extend([f"{ratio:.6f}", str(tagged_count), str(total)])
                rows.append((gene_id, replicate_passes[gene_id], mean_ratio, mean_tagged, mean_total, per_sample))
        rows.sort(key=lambda item: (item[2], item[3]), reverse=True)
        for gene_id, passes, mean_ratio, mean_tagged, mean_total, per_sample in rows[:10]:
            writer.writerow(
                [
                    gene_id,
                    passes,
                    f"{mean_ratio:.6f}",
                    f"{mean_tagged:.2f}",
                    f"{mean_total:.2f}",
                    *per_sample,
                ]
            )

    print(f"Outputs written to {output_dir}")


if __name__ == "__main__":
    main()
