#!/usr/bin/env Rscript
# Copyright 2026 Martin Sikora <martin.sikora@sund.ku.dk>
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(scico)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "Windowed coverage table (chromosome, pop_id, window_start, window_end, avg_n_seg)"
)

parser$add_argument("-c", "--color_file",
  action = "store",
  dest = "color_file",
  help = "Population color map table with pop_id and color/fill (used for ordering)"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output plot basename; writes one PNG per chromosome as <basename>.<chrom>.png"
)

parser$add_argument("--value_col",
  action = "store",
  dest = "value_col",
  default = "avg_n_seg",
  help = "Value column to plot as heatmap fill [default %(default)s]"
)

parser$add_argument("--min_cluster_n",
  action = "store",
  dest = "min_cluster_n",
  type = "integer",
  default = 4L,
  help = "Minimum number of individuals required per cluster/pop_id [default %(default)s]"
)

parser$add_argument("--chrom_list",
  action = "store",
  dest = "chrom_list",
  default = NULL,
  help = "Optional comma-separated chromosome list to force output files"
)

args <- parser$parse_args()
if (args$min_cluster_n < 1) {
  stop("--min_cluster_n must be >= 1")
}


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

d_cov <- read_tsv(args$in_file, show_col_types = FALSE)
color_map <- read_tsv(args$color_file, show_col_types = FALSE)

if (!all(c("pop_id") %in% colnames(color_map))) {
  stop("color_file must include at least: pop_id")
}
if (!(args$value_col %in% colnames(d_cov))) {
  stop("value column not found in input: ", args$value_col)
}


## --------------------------------------------------
## prepare plot data

cat("__ preparing plot data __\n")

pop_levels <- color_map |>
  pull(pop_id) |>
  as.character()

if ("norm_factor" %in% colnames(d_cov)) {
  if ("group_n" %in% colnames(d_cov)) {
    pop_n <- d_cov |>
      transmute(
        pop_id = as.character(pop_id),
        group_n = suppressWarnings(as.integer(group_n))
      ) |>
      filter(!is.na(pop_id), !is.na(group_n), group_n > 0) |>
      group_by(pop_id) |>
      summarise(cluster_n = max(group_n), .groups = "drop")
  } else if ("norm_unit" %in% colnames(d_cov) && any(as.character(d_cov$norm_unit) == "n_pairs", na.rm = TRUE)) {
    pop_n <- d_cov |>
      transmute(
        pop_id = as.character(pop_id),
        norm_factor = suppressWarnings(as.double(norm_factor))
      ) |>
      filter(!is.na(pop_id), !is.na(norm_factor), norm_factor > 0) |>
      group_by(pop_id) |>
      summarise(
        cluster_n = as.integer(round((1 + sqrt(1 + 8 * max(norm_factor))) / 2)),
        .groups = "drop"
      )
  } else {
    pop_n <- d_cov |>
      transmute(
        pop_id = as.character(pop_id),
        norm_factor = suppressWarnings(as.integer(norm_factor))
      ) |>
      filter(!is.na(pop_id), !is.na(norm_factor), norm_factor > 0) |>
      group_by(pop_id) |>
      summarise(cluster_n = max(norm_factor), .groups = "drop")
  }
} else {
  warning("Input has no norm_factor column; assuming cluster_n = 1 for all clusters")
  pop_n <- tibble::tibble(pop_id = unique(pop_levels), cluster_n = 1L)
}

eligible_pops <- pop_n |>
  filter(cluster_n >= args$min_cluster_n) |>
  pull(pop_id)

pop_levels <- pop_levels[pop_levels %in% eligible_pops]
if (length(pop_levels) == 0) {
  stop("No clusters pass --min_cluster_n filter")
}

d <- d_cov |>
  transmute(
    chromosome = as.character(chromosome),
    pop_id = as.character(pop_id),
    window_start = as.double(window_start),
    window_end = as.double(window_end),
    value = as.double(.data[[args$value_col]])
  ) |>
  filter(
    !is.na(chromosome),
    !is.na(pop_id),
    !is.na(window_start),
    !is.na(window_end),
    !is.na(value)
  ) |>
  mutate(
    pop_id = factor(pop_id, levels = rev(pop_levels)),
    window_mid_mb = (window_start + window_end) / 2e6
  ) |>
  filter(!is.na(pop_id))

n_pop <- d |>
  pull(pop_id) |>
  unique() |>
  length()

chr_levels <- d |>
  pull(chromosome) |>
  unique()

target_chroms <- if (!is.null(args$chrom_list) && nzchar(args$chrom_list)) {
  unique(trimws(strsplit(args$chrom_list, ",", fixed = TRUE)[[1]]))
} else {
  chr_levels
}
target_chroms <- target_chroms[nzchar(target_chroms)]
if (length(target_chroms) == 0) {
  stop("No chromosomes available for plotting")
}

chr_info <- d |>
  group_by(chromosome) |>
  summarise(
    xmin = min(window_mid_mb, na.rm = TRUE),
    xmax = max(window_mid_mb, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    chromosome = factor(chromosome, levels = chr_levels),
    span = pmax(xmax - xmin, 1e-6)
  ) |>
  arrange(chromosome)

if (nrow(chr_info) > 0) {
  rel_widths <- chr_info$span / min(chr_info$span, na.rm = TRUE)
  rel_widths[!is.finite(rel_widths) | rel_widths <= 0] <- 1
  rel_width_map <- stats::setNames(rel_widths, as.character(chr_info$chromosome))
} else {
  rel_width_map <- setNames(numeric(0), character(0))
}
n_chr <- length(target_chroms)

h <- max(4, n_pop %/% 7 + 2)
dpi <- 300
if (!is.finite(h) || is.na(h)) {
  h <- 4
}
height_px <- as.integer(round(h * dpi))
if (!is.finite(height_px) || is.na(height_px) || height_px < 1) {
  height_px <- as.integer(4 * dpi)
}

sanitize_chr <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
chrom_out_file <- function(out_file, chr) {
  ext <- tools::file_ext(out_file)
  stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", out_file) else out_file
  if (!nzchar(ext)) {
    ext <- "png"
  }
  paste0(stem, ".", sanitize_chr(chr), ".", ext)
}

th <- theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.1),
    panel.grid.major = element_line(
      linetype = "dotted",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.1),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    panel.spacing = grid::unit(0.1, "lines")
  )


## --------------------------------------------------
## plot heatmap

cat("__ generating heatmap __\n")

max_png_px <- 30000L

for (i in seq_len(n_chr)) {
  chr <- target_chroms[i]
  d_chr <- d |>
    filter(chromosome == chr)

  out_chr <- chrom_out_file(args$out_file, chr)

  w_rel <- as.numeric(rel_width_map[chr])
  if (length(w_rel) == 0 || is.na(w_rel) || !is.finite(w_rel) || w_rel <= 0) {
    w_rel <- 1
  }
  w_chr <- max(4, w_rel * 5)
  if (!is.finite(w_chr) || is.na(w_chr)) {
    w_chr <- 4
  }
  width_px <- as.integer(round(w_chr * dpi))
  if (!is.finite(width_px) || is.na(width_px) || width_px < 1) {
    width_px <- as.integer(4 * dpi)
  }
  local_h_px <- height_px
  if (isTRUE(width_px > max_png_px) || isTRUE(local_h_px > max_png_px)) {
    scale_down <- max(width_px / max_png_px, local_h_px / max_png_px)
    width_px <- as.integer(max(1, floor(width_px / scale_down)))
    local_h_px <- as.integer(max(1, floor(local_h_px / scale_down)))
  }

  png(out_chr, width = width_px, height = local_h_px, res = dpi)

  if (nrow(d_chr) == 0) {
    p_chr <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No data", size = 4) +
      xlim(0, 1) +
      ylim(0, 1) +
      xlab("Position (Mb)") +
      ylab("") +
      labs(title = chr) +
      th +
      theme(
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
      )
  } else {
    p_chr <- ggplot(d_chr, aes(
      x = window_mid_mb,
      y = pop_id,
      fill = value
    )) +
      geom_tile() +
      scale_fill_scico(palette = "davos", direction = -1) +
      scale_x_continuous(labels = label_number(accuracy = 0.1)) +
      xlab("Position (Mb)") +
      ylab("") +
      labs(title = chr, fill = args$value_col) +
      th +
      theme(
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.position = "right"
      )
  }

  print(p_chr)
  dev.off()
}

cat("__ done! __\n")
