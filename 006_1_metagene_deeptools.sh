#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 006_1_metagene_deeptools.sh \
    --bam-dir <mapping_dir> \
    --gff <annotation.gff> \
    [--complete-bam-subdir complete_reads] \
    [--out-dir <metagene_dir>] \
    [--threads 8] \
    [--normalize BPM] \
    [--skip-zeros]

Description:
  1) Convert BAM to BigWig for raw + complete BAM sets.
  2) Run metagene profile via 006_2_metagene_coverage.py on complete BigWigs.
EOF
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THREADS=8
BAM_DIR=""
GFF=""
OUT_DIR=""
COMPLETE_BAM_SUBDIR="complete_reads"
NORMALIZE="BPM"
SKIP_ZEROS=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir)
      BAM_DIR="$2"
      shift 2
      ;;
    --gff)
      GFF="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --complete-bam-subdir)
      COMPLETE_BAM_SUBDIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --normalize)
      NORMALIZE="$2"
      shift 2
      ;;
    --skip-zeros)
      SKIP_ZEROS=true
      shift
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

if [[ -z "$BAM_DIR" || -z "$GFF" ]]; then
  echo "--bam-dir and --gff are required" >&2
  usage
  exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="${BAM_DIR}/metagene_deeptools"
fi

if [[ ! -d "$BAM_DIR" ]]; then
  echo "BAM directory not found: $BAM_DIR" >&2
  exit 1
fi

if [[ ! -f "$GFF" ]]; then
  echo "GFF not found: $GFF" >&2
  exit 1
fi

require_cmd bamCoverage
require_cmd python

BW_DIR="${OUT_DIR}/bw_${NORMALIZE,,}"
mkdir -p "${BW_DIR}/raw" "${BW_DIR}/complete"

for bam in "${BAM_DIR}"/*.sorted.bam; do
  [[ -f "$bam" ]] || continue
  base=$(basename "$bam" .bam)
  bamCoverage -b "$bam" -o "${BW_DIR}/raw/${base}.bw" --normalizeUsing "$NORMALIZE" -p "$THREADS"
done

COMPLETE_DIR="${BAM_DIR}/${COMPLETE_BAM_SUBDIR}"
if [[ -d "$COMPLETE_DIR" ]]; then
  for bam in "${COMPLETE_DIR}"/*.complete.bam; do
    [[ -f "$bam" ]] || continue
    base=$(basename "$bam" .bam)
    bamCoverage -b "$bam" -o "${BW_DIR}/complete/${base}.bw" --normalizeUsing "$NORMALIZE" -p "$THREADS"
  done
fi

PROFILE_OUT="${OUT_DIR}/profile_complete"
CMD=(
  python "${SCRIPT_DIR}/006_2_metagene_coverage.py"
  --bw-dir "${BW_DIR}/complete"
  --gff "$GFF"
  --out-dir "$PROFILE_OUT"
  --threads "$THREADS"
)

if [[ "$SKIP_ZEROS" == true ]]; then
  CMD+=(--skip-zeros)
fi

"${CMD[@]}"

echo "[DONE] Metagene outputs in: ${OUT_DIR}"
