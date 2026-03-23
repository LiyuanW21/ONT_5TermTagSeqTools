#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 006_3_merge_bigwig.sh \
    --ref-fasta <reference.fa> \
    --bw-root <bw_root_with_raw_complete> \
    --out-root <merged_bw_output> \
    [--sample-pattern '_rep[0-9]+$']

Description:
  Merge replicate BigWig files by average signal.
  Expected subdirs under --bw-root: raw/ and complete/.

Grouping rule:
  key = <sample_base>_<tag_group>
  tag_group inferred from filename: tagged / nontagged / all

Sample base is derived by removing trailing replicate suffix matched by --sample-pattern.
EOF
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

REF_FASTA=""
BW_ROOT=""
OUT_ROOT=""
SAMPLE_PATTERN='_[Rr]ep[0-9]+$|_[0-9]+$'

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref-fasta)
      REF_FASTA="$2"
      shift 2
      ;;
    --bw-root)
      BW_ROOT="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --sample-pattern)
      SAMPLE_PATTERN="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$REF_FASTA" || -z "$BW_ROOT" || -z "$OUT_ROOT" ]]; then
  echo "--ref-fasta, --bw-root, and --out-root are required" >&2
  usage
  exit 1
fi

if [[ ! -f "$REF_FASTA" ]]; then
  echo "Reference FASTA not found: $REF_FASTA" >&2
  exit 1
fi

require_cmd bigWigMerge
require_cmd bedGraphToBigWig
require_cmd samtools
require_cmd awk
require_cmd sed

mkdir -p "$OUT_ROOT/raw" "$OUT_ROOT/complete" "$OUT_ROOT/tmp"

samtools faidx "$REF_FASTA"
CHROM_SIZES="$OUT_ROOT/chrom.sizes"
cut -f1,2 "${REF_FASTA}.fai" > "$CHROM_SIZES"

merge_dir() {
  local in_dir="$1"
  local out_dir="$2"
  declare -A groups

  for bw in "$in_dir"/*.bw; do
    [[ -f "$bw" ]] || continue
    local base
    base=$(basename "$bw" .bw)
    base=${base%.sorted}

    local tag="all"
    if [[ "$base" == *"_tagged"* ]]; then
      tag="tagged"
    elif [[ "$base" == *"_nontagged"* ]]; then
      tag="nontagged"
    fi

    local sample_base="$base"
    sample_base=$(echo "$sample_base" | sed -E "s/_(tagged|nontagged)$//")
    sample_base=$(echo "$sample_base" | sed -E "s/(${SAMPLE_PATTERN})$//")

    local key="${sample_base}_${tag}"
    groups["$key"]="${groups[$key]-} $bw"
  done

  for key in "${!groups[@]}"; do
    local files=( ${groups[$key]} )
    local count=${#files[@]}
    [[ $count -gt 0 ]] || continue

    local tmp_bg="$OUT_ROOT/tmp/${key}.bedGraph"
    local tmp_avg="$OUT_ROOT/tmp/${key}.avg.bedGraph"
    local out_bw="$out_dir/${key}.merged.bw"

    if [[ $count -eq 1 ]]; then
      cp "${files[0]}" "$out_bw"
      continue
    fi

    bigWigMerge "${files[@]}" "$tmp_bg"
    awk -v n="$count" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/n}' "$tmp_bg" > "$tmp_avg"
    bedGraphToBigWig "$tmp_avg" "$CHROM_SIZES" "$out_bw"
    rm -f "$tmp_bg" "$tmp_avg"
  done
}

merge_dir "$BW_ROOT/raw" "$OUT_ROOT/raw"
merge_dir "$BW_ROOT/complete" "$OUT_ROOT/complete"

echo "[DONE] Merged BigWigs: $OUT_ROOT"
