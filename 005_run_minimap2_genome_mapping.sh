#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 005_run_minimap2_genome_mapping.sh \
    --input-dir <cutadapt_fastq_dir> \
    --ref-genome <reference.fa> \
    [--output-dir <mapping_out_dir>] \
    [--threads 8] \
    [--preset map-ont]

Description:
  Map FASTQ files to genome by minimap2, then sort/index BAM and summarize mapping.
EOF
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

THREADS=8
INPUT_DIR=""
REF_GENOME=""
OUT_DIR=""
PRESET="map-ont"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    --ref-genome)
      REF_GENOME="$2"
      shift 2
      ;;
    --output-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --preset)
      PRESET="$2"
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

if [[ -z "$INPUT_DIR" || -z "$REF_GENOME" ]]; then
  echo "--input-dir and --ref-genome are required" >&2
  usage
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

if [[ ! -f "$REF_GENOME" ]]; then
  echo "Reference genome not found: $REF_GENOME" >&2
  exit 1
fi

require_cmd minimap2
require_cmd samtools
require_cmd awk

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="${INPUT_DIR}/mapping_out"
fi
mkdir -p "$OUT_DIR"

SUMMARY_TSV="${OUT_DIR}/mapping_rate_summary.tsv"
echo -e "sample\tinput_fq\ttotal_reads\tprimary_mapped_reads\tprimary_mapping_rate(%)" > "$SUMMARY_TSV"

shopt -s nullglob
FILES=("${INPUT_DIR}"/*.fq "${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fq.gz "${INPUT_DIR}"/*.fastq.gz)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No FASTQ files found in: ${INPUT_DIR}" >&2
  exit 1
fi

for fq in "${FILES[@]}"; do
  base=$(basename "$fq")
  sample=${base%.*}
  sample=${sample%.*}

  sam="${OUT_DIR}/${sample}.sam"
  bam="${OUT_DIR}/${sample}.bam"
  sort_bam="${OUT_DIR}/${sample}.sorted.bam"
  flagstat_txt="${OUT_DIR}/${sample}.flagstat.txt"
  log_txt="${OUT_DIR}/${sample}.minimap2.log"

  echo "[INFO] mapping: ${base}"

  minimap2 -t "$THREADS" -ax "$PRESET" "$REF_GENOME" "$fq" > "$sam" 2> "$log_txt"
  samtools view -@ "$THREADS" -bS "$sam" > "$bam"
  samtools sort -@ "$THREADS" -o "$sort_bam" "$bam"
  samtools index -@ "$THREADS" "$sort_bam"
  samtools flagstat -@ "$THREADS" "$sort_bam" > "$flagstat_txt"

  total_reads=$(awk 'NR==1 {print $1}' "$flagstat_txt")
  primary_mapped_reads=$(awk '/primary mapped/ {print $1; exit}' "$flagstat_txt")
  primary_mapping_rate=$(awk '/primary mapped/ {gsub(/[()%]/,"",$5); print $5; exit}' "$flagstat_txt")

  echo -e "${sample}\t${fq}\t${total_reads}\t${primary_mapped_reads}\t${primary_mapping_rate}" >> "$SUMMARY_TSV"

  rm -f "$sam" "$bam"
done

echo "[DONE] Summary: ${SUMMARY_TSV}"
