#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
from collections import defaultdict


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


def write_cds_bed(gff_path, bed_path):
    transcripts = {}
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
            attr_map = parse_attributes(attrs)
            parent = (
                attr_map.get("Parent")
                or attr_map.get("transcript_id")
                or attr_map.get("gene_id")
                or attr_map.get("ID")
            )
            if not parent:
                continue
            parent = parent.split(",")[0]
            start_i = int(start)
            end_i = int(end)
            key = f"{seqid}|{parent}|{strand}"
            if key not in transcripts:
                transcripts[key] = {
                    "seqid": seqid,
                    "strand": strand,
                    "start": start_i,
                    "end": end_i,
                    "name": parent,
                }
            else:
                transcripts[key]["start"] = min(transcripts[key]["start"], start_i)
                transcripts[key]["end"] = max(transcripts[key]["end"], end_i)

    with open(bed_path, "w", encoding="utf-8") as out:
        for info in sorted(transcripts.values(), key=lambda x: (x["seqid"], x["start"])):
            # BED is 0-based, half-open
            out.write(
                "\t".join(
                    [
                        info["seqid"],
                        str(info["start"] - 1),
                        str(info["end"]),
                        info["name"],
                        "0",
                        info["strand"],
                    ]
                )
                + "\n"
            )


def collect_files(data_dir, suffix):
    return sorted(
        [
            os.path.join(data_dir, f)
            for f in os.listdir(data_dir)
            if f.endswith(suffix) and os.path.isfile(os.path.join(data_dir, f))
        ]
    )


def sample_label(path, suffix):
    base = os.path.basename(path)
    if suffix and base.endswith(suffix):
        base = base[: -len(suffix)]
    return base


def run(cmd):
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Metagene coverage using deepTools.")
    parser.add_argument("--bw-dir", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--bins-cds", type=int, default=100)
    parser.add_argument("--bins-flank", type=int, default=25)
    parser.add_argument("--suffix", default=".bw")
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--skip-zeros", action="store_true")
    args = parser.parse_args()

    if shutil.which("computeMatrix") is None or shutil.which("plotProfile") is None:
        raise SystemExit("deepTools not found: please ensure computeMatrix and plotProfile are in PATH.")

    os.makedirs(args.out_dir, exist_ok=True)

    if args.flank % args.bins_flank != 0:
        raise SystemExit("flank must be divisible by bins-flank to define an integer bin size.")

    bin_size = args.flank // args.bins_flank
    region_body_length = args.bins_cds * bin_size

    bed_path = os.path.join(args.out_dir, "cds_regions.bed")
    write_cds_bed(args.gff, bed_path)

    bw_paths = collect_files(args.bw_dir, args.suffix)
    if not bw_paths:
        raise SystemExit("No BigWig files found in bw-dir with given suffix.")

    labels = [sample_label(p, args.suffix) for p in bw_paths]

    matrix_path = os.path.join(args.out_dir, "metagene_matrix.gz")
    profile_plot = os.path.join(args.out_dir, "metagene_profile.png")
    profile_data = os.path.join(args.out_dir, "metagene_profile.tsv")

    cmd_matrix = [
        "computeMatrix",
        "scale-regions",
        "-S",
        *bw_paths,
        "-R",
        bed_path,
        "-b",
        str(args.flank),
        "-a",
        str(args.flank),
        "--regionBodyLength",
        str(region_body_length),
        "--binSize",
        str(bin_size),
        "-p",
        str(args.threads),
        "-o",
        matrix_path,
        "--samplesLabel",
        *labels,
    ]
    if args.skip_zeros:
        cmd_matrix.append("--skipZeros")

    cmd_plot = [
        "plotProfile",
        "-m",
        matrix_path,
        "-out",
        profile_plot,
        "--outFileNameData",
        profile_data,
        "--perGroup",
        "--plotTitle",
        "Metagene coverage (CDS scaled with flanks)",
    ]

    run(cmd_matrix)
    run(cmd_plot)


if __name__ == "__main__":
    main()
