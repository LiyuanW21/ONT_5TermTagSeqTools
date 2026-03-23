# ONT 5'Tag-Seq Tools (Reusable SOP)

Reusable analysis toolkit for ONT 5' tag-seq data.

This repository supports NAD-tagSeq and other 5' tagging designs by parameterizing:
- tag sequence and reverse-complement sequence
- reference genome and annotation
- input FASTQ data and sample set

---

## What this repo now provides

- 7-step SOP with one script per step
- hardcoded paths removed from major scripts
- sample assumptions removed from counting and plotting workflows
- example demo dataset extracted from your request: 2,000 reads

Example data generated at:

`example_data/1306_1_pass.2k.fq`

---

## Pipeline steps (logical execution order)

You requested step numbering as:
1) pychopper 2) cutadapt tuning 3) cutadapt 4) symmetric filtering 5) minimap2 6) deepTools+R 7) tagging genes+Venn.

Because symmetric filtering requires BAM alignment, the runnable order is:

`1 -> 2 -> 3 -> 5 -> 4 -> 6 -> 7`

---

## Script map

| Step | Script | Notes |
|---|---|---|
| 1 | `001_run_pychopper_pipeline.sh` | Two-stage full-length cDNA filtering (edlib + phmm) |
| 2 | `002_cutadapt_grid_tuning.py` | Grid-search cutadapt `-e/-O` over samples |
| 3 | `003_run_cutadapt_tag_pipeline.sh` | Split tagged / nontagged reads |
| 4 | `004_complete_read_filtering.py` | Symmetric complete-read filtering (uses BAM + GFF) |
| 5 | `005_run_minimap2_genome_mapping.sh` | FASTQ -> sorted/indexed BAM + mapping summary |
| 6 | `006_1_metagene_deeptools.sh` + `006_2_metagene_coverage.py` + `006_3_merge_bigwig.sh` + `006_4_plot_metagene_profile.R` | BigWig generation, metagene profile, merge, plotting |
| 7 | `007_1_gene_counts_nad_ratio.py` + `007_2_plot_nad_capping_venn.R` | Tagging genes, tagging ratio, top genes, Venn plots |

---

## Dependencies

### CLI tools
- pychopper
- cutadapt
- minimap2
- samtools
- bedtools
- deepTools (`bamCoverage`, `computeMatrix`, `plotProfile`)
- UCSC tools (`bigWigMerge`, `bedGraphToBigWig`)

### Python
- pysam
- pandas
- matplotlib

### R
- ggplot2
- dplyr
- stringr
- tibble
- VennDiagram
- scales

---

## Inputs you need

- raw FASTQ(s)
- reference genome FASTA
- annotation GFF (for filtering/metagene)
- gene BED (for gene-level counting)
- tag config file (`TAG_SEQ` + `TAG_RC`)

---

## Quickstart with your 2k demo data

Below is a minimal run template. Replace reference/annotation paths with your own.

```bash
PROJECT=/home/liyuan/001_analysis/009_NAD_git/ONT_5TermTagSeqTools
DEMO_FQ=${PROJECT}/example_data/1306_1_pass.2k.fq

REF_FA=/path/to/reference.fa
REF_GFF=/path/to/annotation.gff
GENE_BED=/path/to/genes.bed

TAG_FILE=${PROJECT}/example_data/tag_config.example.txt

OUT=${PROJECT}/demo_run
mkdir -p ${OUT}

# Step 1
bash ${PROJECT}/001_run_pychopper_pipeline.sh \
  --input-fq ${DEMO_FQ} \
  --output-dir ${OUT}/001_pychopper \
  --threads 8 --kit PCS109 --min-qual 7.0 --min-len 50

# Step 2
python ${PROJECT}/002_cutadapt_grid_tuning.py \
  --input-dir ${OUT}/001_pychopper \
  --output-dir ${OUT}/002_tuning \
  --tag-file ${TAG_FILE} \
  --threads 8 \
  --errors 0.15 0.2 0.3 \
  --overlap-min 5 --overlap-max 39

# Step 3 (example uses one full-length FASTQ from step1)
bash ${PROJECT}/003_run_cutadapt_tag_pipeline.sh \
  --output-dir ${OUT}/003_tag \
  --tag-file ${TAG_FILE} \
  --threads 8 --error-rate 0.2 --overlap 12 \
  ${OUT}/001_pychopper/1306_1_pass.2k_final_full_length.fq

# Step 5 (before step 4)
bash ${PROJECT}/005_run_minimap2_genome_mapping.sh \
  --input-dir ${OUT}/003_tag \
  --ref-genome ${REF_FA} \
  --output-dir ${OUT}/004_mapping \
  --threads 8 --preset map-ont

# Step 4
python ${PROJECT}/004_complete_read_filtering.py \
  --bam-dir ${OUT}/004_mapping \
  --gff ${REF_GFF} \
  --out-dir ${OUT}/004_complete_reads \
  --mapq 10 --cds-overlap 50 --nad-tss-window 50 --threads 4

# Step 6
bash ${PROJECT}/006_1_metagene_deeptools.sh \
  --bam-dir ${OUT}/004_mapping \
  --gff ${REF_GFF} \
  --complete-bam-subdir ../004_complete_reads \
  --out-dir ${OUT}/005_metagene \
  --threads 8 --normalize BPM

# Optional merge
bash ${PROJECT}/006_3_merge_bigwig.sh \
  --ref-fasta ${REF_FA} \
  --bw-root ${OUT}/005_metagene/bw_bpm \
  --out-root ${OUT}/005_metagene/bw_bpm_merged

# Plot metagene
Rscript ${PROJECT}/006_4_plot_metagene_profile.R \
  --raw-tsv ${OUT}/005_metagene/profile_complete/metagene_profile.tsv \
  --complete-tsv ${OUT}/005_metagene/profile_complete/metagene_profile.tsv \
  --out-dir ${OUT}/005_metagene/plots \
  --sample-ids 1306_1_pass.2k_final_full_length

# Step 7
python ${PROJECT}/007_1_gene_counts_nad_ratio.py \
  --bam-dir ${OUT}/004_complete_reads \
  --bed ${GENE_BED} \
  --output-dir ${OUT}/007_tagging \
  --min-tagged-reads 2 --min-tagging-ratio 0.002 --min-replicates 1

Rscript ${PROJECT}/007_2_plot_nad_capping_venn.R \
  --input-dir ${OUT}/007_tagging/tagging_by_sample \
  --output-dir ${OUT}/007_tagging/plots \
  --samples 1306_1_pass.2k_final_full_length
```

---

## Per-script interfaces

### 001_run_pychopper_pipeline.sh

```bash
bash 001_run_pychopper_pipeline.sh \
  --input-fq <fastq> \
  [--output-dir <dir>] [--threads 8] [--kit PCS109] \
  [--min-qual 7.0] [--min-len 50]
```

### 002_cutadapt_grid_tuning.py

```bash
python 002_cutadapt_grid_tuning.py \
  --input-dir <dir> --output-dir <dir> \
  --tag-file <tag_config.txt> \
  [--samples s1.fq s2.fq ...] [--threads 8] \
  [--errors 0.15 0.2 0.3] [--overlap-min 5] [--overlap-max 39]
```

### 003_run_cutadapt_tag_pipeline.sh

```bash
bash 003_run_cutadapt_tag_pipeline.sh \
  --output-dir <dir> --tag-file <tag_config.txt> \
  [--threads 8] [--error-rate 0.2] [--overlap 12] \
  <input1.fq> [input2.fq ...]
```

### 005_run_minimap2_genome_mapping.sh

```bash
bash 005_run_minimap2_genome_mapping.sh \
  --input-dir <dir> --ref-genome <fa> \
  [--output-dir <dir>] [--threads 8] [--preset map-ont]
```

### 004_complete_read_filtering.py

```bash
python 004_complete_read_filtering.py \
  --bam-dir <dir> --gff <gff> --out-dir <dir> \
  [--mapq 10] [--cds-overlap 50] [--nad-tss-window 50] [--threads 4]
```

### 006_1_metagene_deeptools.sh

```bash
bash 006_1_metagene_deeptools.sh \
  --bam-dir <dir> --gff <gff> \
  [--complete-bam-subdir complete_reads] [--out-dir <dir>] \
  [--threads 8] [--normalize BPM] [--skip-zeros]
```

### 006_3_merge_bigwig.sh

```bash
bash 006_3_merge_bigwig.sh \
  --ref-fasta <fa> --bw-root <dir> --out-root <dir> \
  [--sample-pattern '_[Rr]ep[0-9]+$|_[0-9]+$']
```

### 006_4_plot_metagene_profile.R

```bash
Rscript 006_4_plot_metagene_profile.R \
  --raw-tsv <tsv> --complete-tsv <tsv> --out-dir <dir> \
  [--out-prefix name] [--sample-ids id1,id2,id3] \
  [--bins-flank 25] [--bins-cds 100]
```

### 007_1_gene_counts_nad_ratio.py

```bash
python 007_1_gene_counts_nad_ratio.py \
  --bam-dir <dir> --bed <bed> --output-dir <dir> \
  [--replicate-keys rep1 rep2 rep3] \
  [--replicate-regex '^(.+?)(?:_pass.*)?$'] \
  [--min-tagged-reads 2] [--min-tagging-ratio 0.002] [--min-replicates 2]
```

### 007_2_plot_nad_capping_venn.R

```bash
Rscript 007_2_plot_nad_capping_venn.R \
  --input-dir <dir> --output-dir <dir> --samples s1,s2,s3 \
  [--suffix .tagging.tsv] [--gene-col gene_id] [--out-prefix tagging_genes_venn]
```

### tag config file format

Recommended format:

```text
TAG_SEQ=CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT
TAG_RC=AGGTTCAGGTTCAGGTTCAGGTTCAGGTTCAGGTTCAGG
```

Also supported:

```text
<TAG_SEQ>
<TAG_RC>
```

---

## Output highlights (step 7)

- `gene_counts_all_samples.tsv`
- `tagging_genes_by_sample.tsv`
- `tagging_genes_replicates.tsv`
- `tagging_ratio_top10_replicates.tsv`
- `tagging_by_sample/*.tagging.tsv`
- Venn plots (`png/pdf/eps/svg`)

---

## GitHub release checklist

1. `git add .`
2. `git commit -m "refactor: parameterize ONT 5' tag-seq SOP scripts"`
3. `git push`
4. add release tag (optional): `git tag v0.2.0 && git push --tags`
