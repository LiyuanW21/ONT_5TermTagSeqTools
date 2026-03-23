#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 003_run_cutadapt_tag_pipeline.sh \
    --output-dir <dir> \
    --tag-file <tag_config.txt> \
    [--threads 8] \
    [--error-rate 0.2] \
    [--overlap 12] \
    <input1.fq> [input2.fq ...]

Description:
  Separate tagged and nontagged reads using cutadapt.
EOF
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

OUT_DIR=""
TAG_FILE=""
TAG_SEQ=""
TAG_RC=""
THREADS=8
ERROR_RATE=0.2
OVERLAP=12

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --tag-file)
      TAG_FILE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --error-rate)
      ERROR_RATE="$2"
      shift 2
      ;;
    --overlap)
      OVERLAP="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if [[ -z "$OUT_DIR" || -z "$TAG_FILE" ]]; then
  echo "--output-dir and --tag-file are required" >&2
  usage
  exit 1
fi

if [[ ! -f "$TAG_FILE" ]]; then
  echo "Tag file not found: $TAG_FILE" >&2
  exit 1
fi

while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%%#*}"
  line="${line%$'\r'}"
  [[ -z "$line" ]] && continue
  if [[ "$line" == TAG_SEQ=* ]]; then
    TAG_SEQ="${line#TAG_SEQ=}"
  elif [[ "$line" == TAG_RC=* ]]; then
    TAG_RC="${line#TAG_RC=}"
  fi
done < "$TAG_FILE"

if [[ -z "$TAG_SEQ" || -z "$TAG_RC" ]]; then
  mapfile -t _tag_lines < <(grep -v '^#' "$TAG_FILE" | sed '/^$/d')
  if [[ ${#_tag_lines[@]} -ge 2 ]]; then
    TAG_SEQ="${_tag_lines[0]}"
    TAG_RC="${_tag_lines[1]}"
  fi
fi

if [[ -z "$TAG_SEQ" || -z "$TAG_RC" ]]; then
  echo "Tag file format invalid. Provide TAG_SEQ=... and TAG_RC=..., or first two non-empty lines as TAG_SEQ then TAG_RC." >&2
  exit 1
fi

if [[ ${#POSITIONAL[@]} -lt 1 ]]; then
  echo "At least one input FASTQ is required" >&2
  usage
  exit 1
fi

require_cmd cutadapt
require_cmd grep
require_cmd sed

mkdir -p "$OUT_DIR"
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"

echo -e "Sample\tInput_File\tTagged_Reads\tNonTagged_Reads\tTrim_Percentage" > "$OUT_DIR/tag_summary.tsv"

for INPUT_FQ in "${POSITIONAL[@]}"; do
  if [[ ! -f "$INPUT_FQ" ]]; then
    echo "Input FASTQ not found: $INPUT_FQ" >&2
    exit 1
  fi

  base=$(basename "$INPUT_FQ")
  sample=${base%%.*}
  tagged_fq="$OUT_DIR/${sample}_tagged.fq"
  nontagged_fq="$OUT_DIR/${sample}_nontagged.fq"
  log="$LOG_DIR/${sample}.cutadapt.log"

  cutadapt \
    -j "$THREADS" \
    -g "$TAG_SEQ" \
    -a "$TAG_RC" \
    -e "$ERROR_RATE" \
    -O "$OVERLAP" \
    -o "$tagged_fq" \
    --untrimmed-output "$nontagged_fq" \
    "$INPUT_FQ" > "$log" 2>&1

  tagged_count=$(grep -c '^@' "$tagged_fq" || true)
  nontagged_count=$(grep -c '^@' "$nontagged_fq" || true)
  trim_percent=$(grep 'Reads with adapters:' "$log" | sed -E 's/.*\(([^)]*)\).*/\1/' || echo '0%')

  echo -e "${sample}\t${INPUT_FQ}\t${tagged_count}\t${nontagged_count}\t${trim_percent}" >> "$OUT_DIR/tag_summary.tsv"
  echo "[DONE] ${sample}: tagged=${tagged_count}, nontagged=${nontagged_count}, trim=${trim_percent}"
done

echo "[DONE] Summary: $OUT_DIR/tag_summary.tsv"
