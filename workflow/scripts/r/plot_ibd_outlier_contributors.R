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
  library(tidyr)
  library(ggplot2)
  library(scales)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "Contributor table from ibd_window_outliers.R --contrib_file"
)

parser$add_argument("-c", "--color_file",
  action = "store",
  dest = "color_file",
  help = "Population color map table with pop_id and optional color/fill/shape"
)

parser$add_argument("-x", "--x_range_file",
  action = "store",
  dest = "x_range_file",
  default = NULL,
  help = "Optional window table used to set full chromosome x-range"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output plot basename; writes one PNG per chromosome as <basename>.<chrom>.png"
)

parser$add_argument("--chrom_list",
  action = "store",
  dest = "chrom_list",
  default = NULL,
  help = "Optional comma-separated chromosome list to force output files"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

d <- read_tsv(args$in_file, show_col_types = FALSE)
color_map <- read_tsv(args$color_file, show_col_types = FALSE)
x_range_tbl <- if (!is.null(args$x_range_file)) read_tsv(args$x_range_file, show_col_types = FALSE) else NULL

req_cols <- c(
  "chromosome", "window_start", "window_end",
  "pop_id", "contribution_frac", "contribution_rank", "score"
)
missing_cols <- setdiff(req_cols, colnames(d))
if (length(missing_cols) > 0) {
  stop("contrib table missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (!("pop_id" %in% colnames(color_map))) {
  stop("color_file must include at least: pop_id")
}


## --------------------------------------------------
## prepare plot data

cat("__ preparing plot data __\n")

if ("shape" %in% colnames(color_map)) {
  color_map <- color_map |>
    mutate(shape = suppressWarnings(as.integer(shape)))
} else {
  color_map <- color_map |>
    mutate(shape = 16L)
}

if (!("color" %in% colnames(color_map))) {
  color_map <- color_map |>
    mutate(color = "#2c3e50")
}

pop_levels <- color_map |>
  pull(pop_id) |>
  as.character()

d_win <- d |>
  transmute(
    chromosome = as.character(chromosome),
    window_start = as.integer(window_start),
    window_end = as.integer(window_end),
    pop_id = as.character(pop_id),
    contribution_frac = as.double(contribution_frac),
    contribution_rank = as.integer(contribution_rank),
    score = as.double(score)
  ) |>
  filter(
    !is.na(chromosome),
    !is.na(window_start),
    !is.na(window_end),
    !is.na(pop_id),
    !is.na(contribution_rank),
    !is.na(score)
  ) |>
  mutate(
    contribution_frac = if_else(is.na(contribution_frac), 0, contribution_frac),
    window_key = paste(chromosome, window_start, window_end, sep = ":")
  )

if (nrow(d_win) == 0) {
  stop("No contributor rows available for plotting")
}

windows <- d_win |>
  distinct(window_key, chromosome, window_start, window_end, score) |>
  arrange(desc(score), chromosome, window_start, window_end)

d_plot <- d_win |>
  inner_join(select(windows, window_key), by = "window_key") |>
  mutate(
    window_mid_mb = (window_start + window_end) / 2e6,
    pop_id = factor(pop_id, levels = rev(pop_levels))
  ) |>
  filter(!is.na(pop_id))

if (nrow(d_plot) == 0) {
  if (is.null(args$chrom_list) || !nzchar(args$chrom_list)) {
    stop("No contributor rows overlap populations in color_map")
  }
}

shape_values <- color_map |>
  transmute(
    pop_id = as.character(pop_id),
    shape = if_else(is.na(shape), 16L, shape)
  ) |>
  tibble::deframe()

color_values <- color_map |>
  transmute(
    pop_id = as.character(pop_id),
    color = if_else(is.na(color), "#2c3e50", color)
  ) |>
  tibble::deframe()

chr_levels <- if (nrow(d_plot) > 0) {
  d_plot |>
    pull(chromosome) |>
    unique()
} else {
  character(0)
}

target_chroms <- if (!is.null(args$chrom_list) && nzchar(args$chrom_list)) {
  unique(trimws(strsplit(args$chrom_list, ",", fixed = TRUE)[[1]]))
} else {
  chr_levels
}
target_chroms <- target_chroms[nzchar(target_chroms)]
if (length(target_chroms) == 0) {
  stop("No chromosomes available for plotting")
}

chr_info <- d_plot |>
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

if (!is.null(x_range_tbl)) {
  if (!all(c("chromosome", "window_start", "window_end") %in% colnames(x_range_tbl))) {
    stop("x_range_file must include columns: chromosome, window_start, window_end")
  }
  chr_extent <- x_range_tbl |>
    transmute(
      chromosome = as.character(chromosome),
      window_start = as.double(window_start),
      window_end = as.double(window_end),
      x_min = window_start / 1e6,
      x_max = window_end / 1e6
    ) |>
    filter(!is.na(chromosome), !is.na(x_min), !is.na(x_max)) |>
    group_by(chromosome) |>
    summarise(
      xmin = min(x_min, na.rm = TRUE),
      xmax = max(x_max, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(span = pmax(xmax - xmin, 1e-6))

  if (nrow(chr_extent) > 0) {
    chr_info <- chr_extent |>
      mutate(chromosome = factor(chromosome, levels = target_chroms)) |>
      filter(!is.na(chromosome)) |>
      arrange(chromosome)
  }
}

if (nrow(chr_info) > 0) {
  rel_widths <- chr_info$span / min(chr_info$span, na.rm = TRUE)
  rel_widths[!is.finite(rel_widths) | rel_widths <= 0] <- 1
  rel_width_map <- stats::setNames(rel_widths, as.character(chr_info$chromosome))
} else {
  rel_width_map <- setNames(numeric(0), character(0))
}
n_chr <- length(target_chroms)

n_pop <- d_plot |>
  pull(pop_id) |>
  unique() |>
  length()

h <- max(2.6, n_pop / 14 + 1.2)
dpi <- 300
if (!is.finite(h) || is.na(h)) {
  h <- 2.6
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
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    panel.spacing = grid::unit(0.1, "lines"),
    legend.position = "right"
  )


## --------------------------------------------------
## plot

cat("__ generating contributor plot __\n")

max_png_px <- 30000L

for (i in seq_len(n_chr)) {
  chr <- target_chroms[i]
  d_chr <- d_plot |>
    filter(chromosome == chr)
  chr_limits <- chr_info |>
    filter(as.character(chromosome) == chr)

  out_chr <- chrom_out_file(args$out_file, chr)

  w_rel <- as.numeric(rel_width_map[chr])
  if (length(w_rel) == 0 || is.na(w_rel) || !is.finite(w_rel) || w_rel <= 0) {
    w_rel <- 1
  }
  w_chr <- max(8, w_rel * 11)
  if (!is.finite(w_chr) || is.na(w_chr)) {
    w_chr <- 8
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
      xlab("") +
      ylab("Clusters (pop_id)") +
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
      y = pop_id
    )) +
      geom_tile(
        aes(fill = contribution_frac),
        color = "grey85",
        linewidth = 0.1
      ) +
      geom_point(
        aes(shape = pop_id, color = pop_id),
        size = 1.8,
        stroke = 0.2
      ) +
      scale_fill_gradient(
        low = "#f7fbff",
        high = "#08519c",
        labels = label_percent(accuracy = 1)
      ) +
      scale_shape_manual(values = shape_values, guide = "none") +
      scale_color_manual(values = color_values, guide = "none") +
      xlab("") +
      ylab("Clusters (pop_id)") +
      labs(title = chr, fill = "Contribution share") +
      th +
      theme(
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.position = "right"
      )
  }

  if (nrow(chr_limits) > 0 && is.finite(chr_limits$xmin[1]) && is.finite(chr_limits$xmax[1]) && chr_limits$xmax[1] > chr_limits$xmin[1]) {
    p_chr <- p_chr + scale_x_continuous(limits = c(chr_limits$xmin[1], chr_limits$xmax[1]))
  }

  print(p_chr)
  dev.off()
}

cat("__ done! __\n")
