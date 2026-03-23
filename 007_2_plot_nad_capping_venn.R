#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    install.packages("VennDiagram", repos = "https://cloud.r-project.org")
  }
  library(VennDiagram)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    input_dir = NULL,
    output_dir = NULL,
    samples = NULL,
    suffix = ".tagging.tsv",
    file_gene_col = "gene_id",
    out_prefix = "tagging_genes_venn"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--input-dir") {
      out$input_dir <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--output-dir") {
      out$output_dir <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--samples") {
      out$samples <- strsplit(args[[i + 1]], ",")[[1]]
      i <- i + 2
    } else if (key == "--suffix") {
      out$suffix <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--gene-col") {
      out$file_gene_col <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--out-prefix") {
      out$out_prefix <- args[[i + 1]]
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  Rscript 007_plot_nad_capping_venn.R \\\n",
        "    --input-dir <tagging_by_sample_dir> \\\n",
        "    --output-dir <plot_dir> \\\n",
        "    --samples sample1,sample2,sample3 \\\n",
        "    [--suffix .tagging.tsv] \\\n",
        "    [--gene-col gene_id] \\\n",
        "    [--out-prefix tagging_genes_venn]\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }

  if (is.null(out$input_dir) || is.null(out$output_dir) || is.null(out$samples)) {
    stop("--input-dir, --output-dir, --samples are required")
  }
  out
}

read_gene_list <- function(input_dir, sample_id, suffix, gene_col) {
  file_path <- file.path(input_dir, paste0(sample_id, suffix))
  if (!file.exists(file_path)) {
    stop("Missing input file: ", file_path)
  }
  df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!(gene_col %in% colnames(df))) {
    stop("Missing column '", gene_col, "' in ", file_path)
  }
  unique(df[[gene_col]])
}

args <- parse_args()

if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

sets <- lapply(args$samples, function(s) {
  read_gene_list(args$input_dir, s, args$suffix, args$file_gene_col)
})
names(sets) <- args$samples

png_path <- file.path(args$output_dir, paste0(args$out_prefix, ".png"))
pdf_path <- file.path(args$output_dir, paste0(args$out_prefix, ".pdf"))
eps_path <- file.path(args$output_dir, paste0(args$out_prefix, ".eps"))
svg_path <- file.path(args$output_dir, paste0(args$out_prefix, ".svg"))

fill_colors <- scales::hue_pal()(length(args$samples))

venn <- venn.diagram(
  x = sets,
  filename = NULL,
  fill = fill_colors,
  alpha = 0.45,
  cex = 1.1,
  cat.cex = 1.1,
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica",
  margin = 0.08
)

venn_eps <- venn.diagram(
  x = sets,
  filename = NULL,
  fill = fill_colors,
  alpha = 1,
  cex = 1.1,
  cat.cex = 1.1,
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica",
  margin = 0.08
)

png(png_path, width = 1600, height = 1600, res = 200)
grid::grid.draw(venn)
dev.off()

pdf(pdf_path, width = 6, height = 6)
grid::grid.draw(venn)
dev.off()

postscript(
  eps_path,
  width = 6,
  height = 6,
  horizontal = FALSE,
  paper = "special",
  family = "Helvetica"
)
grid::grid.draw(venn_eps)
dev.off()

svg(svg_path, width = 6, height = 6)
grid::grid.draw(venn)
dev.off()

cat("Venn plot written to:\n", png_path, "\n", pdf_path, "\n", eps_path, "\n", svg_path, "\n", sep = "")
