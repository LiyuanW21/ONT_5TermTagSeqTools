#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    raw_tsv = NULL,
    complete_tsv = NULL,
    out_dir = NULL,
    out_prefix = "metagene_profile_merged",
    sample_ids = NULL,
    bins_flank = 25,
    bins_cds = 100,
    plot_width = 11,
    plot_height = 4.5
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--raw-tsv") {
      out$raw_tsv <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--complete-tsv") {
      out$complete_tsv <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--out-dir") {
      out$out_dir <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--out-prefix") {
      out$out_prefix <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--sample-ids") {
      out$sample_ids <- strsplit(args[[i + 1]], ",")[[1]]
      i <- i + 2
    } else if (key == "--bins-flank") {
      out$bins_flank <- as.integer(args[[i + 1]])
      i <- i + 2
    } else if (key == "--bins-cds") {
      out$bins_cds <- as.integer(args[[i + 1]])
      i <- i + 2
    } else if (key == "--plot-width") {
      out$plot_width <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else if (key == "--plot-height") {
      out$plot_height <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  Rscript 006_plot_metagene_profile.R \\\n",
        "    --raw-tsv <raw metagene_profile.tsv> \\\n",
        "    --complete-tsv <complete metagene_profile.tsv> \\\n",
        "    --out-dir <plot output dir> \\\n",
        "    [--out-prefix metagene_profile_merged] \\\n",
        "    [--sample-ids id1,id2,id3] \\\n",
        "    [--bins-flank 25] [--bins-cds 100] \\\n",
        "    [--plot-width 11] [--plot-height 4.5]\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }

  if (is.null(out$raw_tsv) || is.null(out$complete_tsv) || is.null(out$out_dir)) {
    stop("--raw-tsv, --complete-tsv, --out-dir are required")
  }
  out
}

read_profile <- function(tsv_path, dataset_label) {
  if (!file.exists(tsv_path)) {
    stop("Missing file: ", tsv_path)
  }

  lines <- readLines(tsv_path)
  split_cols <- strsplit(lines, "\t", fixed = TRUE)
  max_cols <- max(vapply(split_cols, length, integer(1)))
  split_cols <- lapply(split_cols, function(x) {
    length(x) <- max_cols
    x
  })
  df <- as.data.frame(do.call(rbind, split_cols), stringsAsFactors = FALSE)

  bin_row <- suppressWarnings(as.numeric(df[2, ]))
  numeric_cols <- which(!is.na(bin_row))

  if (length(numeric_cols) == 0) {
    stop("No numeric bins found in: ", tsv_path)
  }
  if (nrow(df) < 3) {
    stop("No data rows found in: ", tsv_path)
  }

  data_list <- lapply(3:nrow(df), function(i) {
    sample_name <- df[i, 1]
    group_name <- df[i, 2]
    values <- suppressWarnings(as.numeric(df[i, numeric_cols]))

    tibble(
      sample = sample_name,
      group = group_name,
      x = seq_along(values),
      value = values
    )
  })

  bind_rows(data_list) %>%
    mutate(
      sample = str_replace(sample, "\\.merged$", ""),
      tag = case_when(
        str_detect(sample, "nontagged") ~ "nontagged",
        str_detect(sample, "tagged") ~ "tagged",
        TRUE ~ "unknown"
      ),
      sample_id = str_replace(sample, "_(tagged|nontagged)$", ""),
      dataset = dataset_label
    )
}

args <- parse_args()

tis_bin <- args$bins_flank
tes_bin <- args$bins_flank + args$bins_cds
end_bin <- args$bins_flank * 2 + args$bins_cds

raw_data <- read_profile(args$raw_tsv, "raw")
complete_data <- read_profile(args$complete_tsv, "complete")

all_data <- bind_rows(raw_data, complete_data)
if (!is.null(args$sample_ids)) {
  all_data <- all_data %>% filter(sample_id %in% args$sample_ids)
}

sample_levels <- if (is.null(args$sample_ids)) sort(unique(all_data$sample_id)) else args$sample_ids

all_data <- all_data %>%
  mutate(
    sample_label = factor(sample_id, levels = sample_levels),
    tag = factor(tag, levels = c("tagged", "nontagged")),
    color_key = paste(sample_label, tag, sep = "_"),
    dataset = factor(dataset, levels = c("raw", "complete"))
  )

keys <- unique(all_data$color_key)
pal <- scales::hue_pal()(length(keys))
color_map <- setNames(pal, keys)
legend_labels <- setNames(gsub("_", " ", keys), keys)

p <- ggplot(all_data, aes(x = x, y = value, color = color_key)) +
  facet_wrap(~dataset, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = c(tis_bin, tes_bin), linetype = "dashed", linewidth = 0.5, color = "#555555") +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(
    limits = c(1, end_bin),
    breaks = c(1, tis_bin, tes_bin, end_bin),
    labels = c("-flank", "TIS", "TES", "+flank"),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = color_map, labels = legend_labels) +
  labs(y = "Normalized coverage", title = "Metagene profile") +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(linewidth = 0.8, colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_rect(fill = "#F5F5F5", color = "black", linewidth = 0.6)
  )

if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

png_path <- file.path(args$out_dir, paste0(args$out_prefix, ".png"))
eps_path <- file.path(args$out_dir, paste0(args$out_prefix, ".eps"))

ggsave(png_path, p, width = args$plot_width, height = args$plot_height, dpi = 300)

if (capabilities("cairo")) {
  ggsave(eps_path, p, width = args$plot_width, height = args$plot_height, device = cairo_ps)
} else {
  ggsave(eps_path, p, width = args$plot_width, height = args$plot_height, device = "eps")
}

message("Saved: ", png_path)
message("Saved: ", eps_path)
