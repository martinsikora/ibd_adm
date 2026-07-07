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
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
  nargs = "+",
  help = "One or more outlier score tables from ibd_window_outliers.R --score_file"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output plot filename (PNG)"
)

parser$add_argument("--point_size",
  action = "store",
  dest = "point_size",
  type = "double",
  default = 0.8,
  help = "Point size for outlier points [default %(default)s]"
)

parser$add_argument("--min_cluster_n",
  action = "store",
  dest = "min_cluster_n",
  type = "integer",
  default = 4L,
  help = "Minimum number of individuals required per cluster/pop_id [default %(default)s]"
)

args <- parser$parse_args()
if (args$min_cluster_n < 1) {
  stop("--min_cluster_n must be >= 1")
}


## --------------------------------------------------
## read and combine

cat("__ reading score tables __\n")

d <- lapply(args$files, function(f) {
  r <- read_tsv(f, show_col_types = FALSE)
  if (!all(c("chromosome", "window_start", "window_end", "score", "threshold", "is_outlier", "method") %in% colnames(r))) {
    stop("score table missing required columns in file: ", f)
  }
  r
}) |>
  bind_rows()


## --------------------------------------------------
## prepare plot data

cat("__ preparing plot data __\n")

d1 <- d |>
  transmute(
    chromosome = as.character(chromosome),
    window_start = as.double(window_start),
    window_end = as.double(window_end),
    score = as.double(score),
    threshold = as.double(threshold),
    is_outlier = as.logical(is_outlier),
    method = as.character(method),
    score_min_cluster_n = if ("min_cluster_n" %in% colnames(d)) as.integer(min_cluster_n) else NA_integer_,
    window_mid_mb = (window_start + window_end) / 2e6
  ) |>
  filter(
    !is.na(chromosome),
    !is.na(window_mid_mb),
    !is.na(score),
    !is.na(threshold),
    !is.na(method)
  )

if (any(!is.na(d1$score_min_cluster_n))) {
  d1 <- d1 |>
    filter(!is.na(score_min_cluster_n), score_min_cluster_n >= args$min_cluster_n)
}
if (nrow(d1) == 0) {
  stop("No score rows pass --min_cluster_n filter")
}

n_methods <- d1 |>
  pull(method) |>
  unique() |>
  length()

n_chr <- d1 |>
  pull(chromosome) |>
  unique() |>
  length()

w <- max(10, n_chr * 4)
h <- max(4, n_methods * 2.8)
dpi <- 300
width_px <- as.integer(round(w * dpi))
height_px <- as.integer(round(h * dpi))
max_png_px <- 30000L
if (width_px > max_png_px || height_px > max_png_px) {
  scale_down <- max(width_px / max_png_px, height_px / max_png_px)
  width_px <- as.integer(max(1, floor(width_px / scale_down)))
  height_px <- as.integer(max(1, floor(height_px / scale_down)))
  warning("Requested image size was too large for Cairo PNG device; downscaled to ", width_px, "x", height_px, " px")
}

th <- theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.1),
    panel.grid.major = element_line(
      linetype = "dotted",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    panel.spacing = grid::unit(0.1, "lines")
  )


## --------------------------------------------------
## plot

cat("__ generating outlier score plot __\n")

png(args$out_file, width = width_px, height = height_px, res = dpi)
p <- ggplot(d1, aes(
  x = window_mid_mb,
  y = score
))
print(p +
  geom_line(linewidth = 0.3, color = "grey30") +
  geom_hline(
    aes(yintercept = threshold),
    linewidth = 0.25,
    linetype = "dashed",
    color = "tomato"
  ) +
  geom_point(
    data = d1 |>
      filter(is_outlier),
    color = "tomato",
    size = args$point_size
  ) +
  facet_grid(method ~ chromosome, scales = "free", space = "free_x") +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  xlab("Position (Mb)") +
  ylab("Outlier score") +
  th)
dev.off()

cat("__ done! __\n")
