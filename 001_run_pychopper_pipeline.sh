#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 001_run_pychopper_pipeline.sh \
    --input-fq <reads.fastq> \
    [--output-dir <out_dir>] \
    [--threads 8] \
    [--kit PCS109] \
    [--min-qual 7.0] \
    [--min-len 50]

Description:
  Two-stage pychopper pipeline for full-length cDNA extraction.
  Stage1: edlib on full input
  Stage2: phmm on stage1 unclassified reads

Outputs:
  <sample>_final_full_length.fq
  <sample>_final_unclassified.fq
  <sample>_summary_report.txt
EOF
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

INPUT_FQ=""
OUTPUT_DIR=""
THREADS=8
KIT="PCS109"
MIN_QUAL=7.0
MIN_LEN=50

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-fq)
      INPUT_FQ="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --kit)
      KIT="$2"
      shift 2
      ;;
    --min-qual)
      MIN_QUAL="$2"
      shift 2
      ;;
    --min-len)
      MIN_LEN="$2"
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

if [[ -z "$INPUT_FQ" ]]; then
  echo "--input-fq is required" >&2
  usage
  exit 1
fi

if [[ ! -f "$INPUT_FQ" ]]; then
  echo "Input FASTQ not found: $INPUT_FQ" >&2
  exit 1
fi

require_cmd pychopper
require_cmd grep
require_cmd bc

INPUT_FQ_PATH="$(realpath "$INPUT_FQ")"
SAMPLE="$(basename "$INPUT_FQ_PATH")"
SAMPLE="${SAMPLE%.*}"
SAMPLE="${SAMPLE%.*}"

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="$(pwd)/pychopper_output_${SAMPLE}_$(date +%Y%m%d_%H%M%S)"
else
  OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"
fi

mkdir -p "$OUTPUT_DIR"

FINAL_FULL_LENGTH="${OUTPUT_DIR}/${SAMPLE}_final_full_length.fq"
FINAL_UNCLASSIFIED="${OUTPUT_DIR}/${SAMPLE}_final_unclassified.fq"
FINAL_SUMMARY_REPORT="${OUTPUT_DIR}/${SAMPLE}_summary_report.txt"

EDLIB_FULL_LENGTH="${OUTPUT_DIR}/${SAMPLE}_edlib_full_length.fq"
EDLIB_UNCLASSIFIED="${OUTPUT_DIR}/${SAMPLE}_edlib_unclassified.fq"
EDLIB_REPORT="${OUTPUT_DIR}/${SAMPLE}_edlib_report.pdf"
EDLIB_STATS="${OUTPUT_DIR}/${SAMPLE}_edlib_stats.txt"
EDLIB_PER_READ_STATS="${OUTPUT_DIR}/${SAMPLE}_edlib_per_read_stats.tsv"

PHMM_FULL_LENGTH="${OUTPUT_DIR}/${SAMPLE}_phmm_full_length.fq"
PHMM_UNCLASSIFIED="${OUTPUT_DIR}/${SAMPLE}_phmm_unclassified.fq"
PHMM_REPORT="${OUTPUT_DIR}/${SAMPLE}_phmm_report.pdf"
PHMM_STATS="${OUTPUT_DIR}/${SAMPLE}_phmm_stats.txt"
PHMM_PER_READ_STATS="${OUTPUT_DIR}/${SAMPLE}_phmm_per_read_stats.tsv"

echo "[INFO] Input:  ${INPUT_FQ_PATH}"
echo "[INFO] Output: ${OUTPUT_DIR}"

INPUT_READS=$(grep -c '^@' "$INPUT_FQ_PATH" || true)

echo "[INFO] Stage 1: pychopper edlib"
pychopper \
  -k "$KIT" \
  -m edlib \
  -t "$THREADS" \
  -Q "$MIN_QUAL" \
  -z "$MIN_LEN" \
  -r "$EDLIB_REPORT" \
  -u "$EDLIB_UNCLASSIFIED" \
  -S "$EDLIB_STATS" \
  -D "$EDLIB_PER_READ_STATS" \
  "$INPUT_FQ_PATH" \
  "$EDLIB_FULL_LENGTH"

EDLIB_COUNT=$(grep -c '^@' "$EDLIB_FULL_LENGTH" || echo "0")
EDLIB_UN_COUNT=$(grep -c '^@' "$EDLIB_UNCLASSIFIED" || echo "0")

if [[ "$EDLIB_UN_COUNT" -gt 0 ]]; then
  echo "[INFO] Stage 2: pychopper phmm"
  pychopper \
    -k "$KIT" \
    -m phmm \
    -t "$THREADS" \
    -Q "$MIN_QUAL" \
    -z "$MIN_LEN" \
    -r "$PHMM_REPORT" \
    -u "$PHMM_UNCLASSIFIED" \
    -S "$PHMM_STATS" \
    -D "$PHMM_PER_READ_STATS" \
    "$EDLIB_UNCLASSIFIED" \
    "$PHMM_FULL_LENGTH"
else
  echo "[INFO] Stage 2 skipped: no unclassified reads"
  : > "$PHMM_FULL_LENGTH"
  cp "$EDLIB_UNCLASSIFIED" "$PHMM_UNCLASSIFIED"
fi

PHMM_COUNT=$(grep -c '^@' "$PHMM_FULL_LENGTH" || echo "0")
PHMM_UN_COUNT=$(grep -c '^@' "$PHMM_UNCLASSIFIED" || echo "0")

cat "$EDLIB_FULL_LENGTH" "$PHMM_FULL_LENGTH" > "$FINAL_FULL_LENGTH"
cp "$PHMM_UNCLASSIFIED" "$FINAL_UNCLASSIFIED"

FINAL_COUNT=$(grep -c '^@' "$FINAL_FULL_LENGTH" || echo "0")
FINAL_PERCENT="0"
if [[ "$INPUT_READS" -gt 0 ]]; then
  FINAL_PERCENT=$(echo "scale=2; $FINAL_COUNT * 100 / $INPUT_READS" | bc)
fi

cat > "$FINAL_SUMMARY_REPORT" <<EOF
Pychopper Two-stage Summary
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Input FASTQ: $INPUT_FQ_PATH
Output Dir:  $OUTPUT_DIR

Input reads:                 $INPUT_READS
Stage1 full-length (edlib):  $EDLIB_COUNT
Stage2 rescued (phmm):       $PHMM_COUNT
Final full-length reads:     $FINAL_COUNT
Final unclassified reads:    $PHMM_UN_COUNT
Final full-length ratio(%):  $FINAL_PERCENT
EOF

echo "[DONE] Final full-length: $FINAL_FULL_LENGTH"
echo "[DONE] Final unclassified: $FINAL_UNCLASSIFIED"
echo "[DONE] Summary report: $FINAL_SUMMARY_REPORT"
