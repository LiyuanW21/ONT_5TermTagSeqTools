#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def count_fastq_reads(file_path: Path) -> int:
    if not file_path.exists() or file_path.stat().st_size == 0:
        return 0
    with file_path.open("rb") as f:
        line_count = sum(1 for _ in f)
    return line_count // 4


def run_cutadapt(
    input_fq: Path,
    tag_seq: str,
    tag_rc: str,
    error_rate: float,
    overlap: int,
    threads: int,
    tmp_dir: Path,
) -> int:
    tagged_fq = tmp_dir / f"{input_fq.stem}_e{error_rate}_O{overlap}.tagged.fq"
    nontagged_fq = tmp_dir / f"{input_fq.stem}_e{error_rate}_O{overlap}.nontagged.fq"

    cmd = [
        "cutadapt",
        "-j",
        str(threads),
        "-g",
        tag_seq,
        "-a",
        tag_rc,
        "-e",
        str(error_rate),
        "-O",
        str(overlap),
        "-o",
        str(tagged_fq),
        "--untrimmed-output",
        str(nontagged_fq),
        str(input_fq),
    ]

    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    tagged_count = count_fastq_reads(tagged_fq)

    if tagged_fq.exists():
        tagged_fq.unlink()
    if nontagged_fq.exists():
        nontagged_fq.unlink()

    return tagged_count


def load_tag_file(tag_file: Path) -> tuple[str, str]:
    if not tag_file.exists():
        raise FileNotFoundError(f"Tag file not found: {tag_file}")

    entries = []
    kv = {}
    with tag_file.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            entries.append(line)
            if "=" in line:
                key, val = line.split("=", 1)
                kv[key.strip().upper()] = val.strip()

    if "TAG_SEQ" in kv and "TAG_RC" in kv:
        return kv["TAG_SEQ"], kv["TAG_RC"]

    if len(entries) >= 2:
        return entries[0], entries[1]

    raise ValueError(
        "Tag file format invalid. Use either:\n"
        "  TAG_SEQ=<sequence>\n"
        "  TAG_RC=<reverse_complement>\n"
        "or two-line format: first line TAG_SEQ, second line TAG_RC."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cutadapt grid tuning for multiple samples")
    parser.add_argument("--input-dir", required=True, help="Directory containing input FASTQ files")
    parser.add_argument("--output-dir", required=True, help="Output directory for CSV and plot")
    parser.add_argument("--tag-file", required=True, help="Tag config file containing TAG_SEQ and TAG_RC")
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Input sample file names in input-dir. If omitted, auto-discover FASTQ files.",
    )
    parser.add_argument("--threads", type=int, default=8, help="cutadapt -j threads")
    parser.add_argument(
        "--errors",
        nargs="+",
        type=float,
        default=[0.15, 0.2, 0.3],
        help="Mismatch rates, e.g., 0.15 0.2 0.3",
    )
    parser.add_argument("--overlap-min", type=int, default=5, help="Minimum overlap (-O)")
    parser.add_argument("--overlap-max", type=int, default=39, help="Maximum overlap (-O)")
    parser.add_argument(
        "--plot-prefix",
        default="cutadapt_grid_tuning",
        help="Prefix for generated CSV/PNG files",
    )
    return parser.parse_args()


def discover_samples(input_dir: Path, samples: list[str] | None) -> list[str]:
    if samples:
        return samples
    patterns = ("*.fq", "*.fastq", "*.fq.gz", "*.fastq.gz")
    discovered: list[Path] = []
    for pattern in patterns:
        discovered.extend(sorted(input_dir.glob(pattern)))
    seen = set()
    result = []
    for path in discovered:
        if path.name not in seen:
            result.append(path.name)
            seen.add(path.name)
    return result


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    tag_file = Path(args.tag_file)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    overlaps = list(range(args.overlap_min, args.overlap_max + 1))
    errors = args.errors
    tag_seq, tag_rc = load_tag_file(tag_file)

    samples = discover_samples(input_dir, args.samples)
    if not samples:
        raise FileNotFoundError("No input FASTQ files found.")

    results = []
    with tempfile.TemporaryDirectory(dir=output_dir) as tmp:
        tmp_dir = Path(tmp)

        for sample in samples:
            input_fq = input_dir / sample
            if not input_fq.exists():
                print(f"[WARN] Missing: {input_fq}", file=sys.stderr)
                continue

            for error_rate in errors:
                for overlap in overlaps:
                    try:
                        tagged = run_cutadapt(
                            input_fq=input_fq,
                            tag_seq=tag_seq,
                            tag_rc=tag_rc,
                            error_rate=error_rate,
                            overlap=overlap,
                            threads=args.threads,
                            tmp_dir=tmp_dir,
                        )
                    except subprocess.CalledProcessError:
                        tagged = None

                    results.append(
                        {
                            "sample": sample,
                            "error_rate": error_rate,
                            "overlap": overlap,
                            "tagged_reads": tagged,
                        }
                    )

    if not results:
        raise RuntimeError("No tuning results generated.")

    df = pd.DataFrame(results)
    csv_path = output_dir / f"{args.plot_prefix}.csv"
    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(14, 8))
    cmap = plt.get_cmap("tab10")
    sample_colors = {sample: cmap(i % 10) for i, sample in enumerate(samples)}
    style_cycle = ["-", "--", ":", "-."]
    error_styles = {err: style_cycle[i % len(style_cycle)] for i, err in enumerate(errors)}

    for sample in samples:
        sub = df[df["sample"] == sample]
        if len(sub) == 0:
            continue
        for error_rate in errors:
            e_sub = sub[sub["error_rate"] == error_rate]
            if len(e_sub) == 0:
                continue
            plt.plot(
                e_sub["overlap"],
                e_sub["tagged_reads"],
                marker="o",
                linewidth=1.5,
                color=sample_colors.get(sample),
                linestyle=error_styles.get(error_rate, "-"),
                label=f"{sample} | e={error_rate}",
            )

    plt.title("Cutadapt grid tuning (all samples)")
    plt.xlabel("Overlap (-O)")
    plt.ylabel("Tagged reads")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.legend(ncol=2, fontsize=9)
    plt.tight_layout()
    plot_path = output_dir / f"{args.plot_prefix}.png"
    plt.savefig(plot_path, dpi=300)

    print(f"Saved CSV: {csv_path}")
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
