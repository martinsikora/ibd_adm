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

## --------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(argparse)
})


## --------------------------------------------------
## functions

window_average_coverage <- function(d_cov, window_size = 1e6) {
  ws <- as.integer(window_size)

  ## vectorized window index computation
  i0 <- floor(d_cov$pos_start / ws)
  i1 <- floor((d_cov$pos_end - 1L) / ws)
  n_windows <- as.integer(i1 - i0 + 1L)

  ## expand each segment to its overlapping windows
  row_idx <- rep(seq_len(nrow(d_cov)), times = n_windows)
  win_offset <- sequence(n_windows) - 1L
  win_i <- rep(i0, times = n_windows) + win_offset

  w_start <- as.integer(win_i * ws)
  w_end <- as.integer((win_i + 1L) * ws)
  ov_start <- pmax(d_cov$pos_start[row_idx], w_start)
  ov_end <- pmin(d_cov$pos_end[row_idx], w_end)
  ov_bp <- ov_end - ov_start

  keep <- ov_bp > 0

  d_expanded <- tibble(
    chromosome = d_cov$chromosome[row_idx[keep]],
    pop_id = d_cov$pop_id[row_idx[keep]],
    window_start = w_start[keep],
    window_end = w_end[keep],
    overlap_bp = ov_bp[keep],
    weighted_n_seg = ov_bp[keep] * d_cov$n_seg[row_idx[keep]]
  )

  d_expanded |>
    group_by(chromosome, pop_id, window_start, window_end) |>
    summarise(
      covered_bp = sum(overlap_bp),
      avg_n_seg = sum(weighted_n_seg) / sum(overlap_bp),
      .groups = "drop"
    ) |>
    arrange(pop_id, chromosome, window_start)
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "Input overlap coverage table (pop_id, chromosome, pos_start, pos_end, n_seg)"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output table with windowed average IBD coverage"
)

parser$add_argument("-w", "--window_size",
  action = "store",
  dest = "window_size",
  type = "integer",
  default = 1000000L,
  help = "Window size in bp [default %(default)s]"
)

parser$add_argument("-s", "--sample_file",
  action = "store",
  dest = "sample_file",
  default = NULL,
  help = "Optional sample map with sample_id/pop_id used for normalization"
)

parser$add_argument("--norm_mode",
  action = "store",
  dest = "norm_mode",
  default = "pop_size",
  help = "Normalization mode: pop_size (unique pairs), none [default %(default)s]"
)

args <- parser$parse_args()
if (!(args$norm_mode %in% c("none", "pop_size"))) {
  stop("--norm_mode must be one of: none, pop_size")
}


## --------------------------------------------------
## read and process

cat("__ reading coverage data __\n")
d_cov <- read_tsv(args$in_file, show_col_types = FALSE)

req_cols <- c("pop_id", "chromosome", "pos_start", "pos_end", "n_seg")
missing_cols <- setdiff(req_cols, colnames(d_cov))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

d_cov <- d_cov |>
  mutate(
    chromosome = as.character(chromosome),
    pop_id = as.character(pop_id),
    pos_start = as.integer(pos_start),
    pos_end = as.integer(pos_end),
    n_seg = as.double(n_seg)
  ) |>
  filter(
    !is.na(pop_id),
    !is.na(chromosome),
    !is.na(pos_start),
    !is.na(pos_end),
    !is.na(n_seg),
    pos_end > pos_start
  )

cat("__ computing windowed averages __\n")
d_win <- window_average_coverage(d_cov, window_size = args$window_size)

if (args$norm_mode == "pop_size") {
  if (is.null(args$sample_file)) {
    stop("--sample_file is required when --norm_mode pop_size")
  }

  cat("__ normalizing by number of unique within-group pairs __\n")
  sample_map <- read_tsv(args$sample_file, show_col_types = FALSE)
  if (!all(c("sample_id", "pop_id") %in% colnames(sample_map))) {
    stop("sample_file must include columns: sample_id, pop_id")
  }

  if ("group" %in% colnames(sample_map)) {
    pop_n <- sample_map |>
      filter(group != "exclude", pop_id != "exclude") |>
      count(pop_id, name = "group_n")
  } else {
    pop_n <- sample_map |>
      filter(pop_id != "exclude") |>
      count(pop_id, name = "group_n")
  }

  d_win <- d_win |>
    left_join(pop_n, by = "pop_id") |>
    mutate(
      group_n = if_else(is.na(group_n) | group_n <= 0, 1L, group_n),
      norm_factor = as.integer(group_n * (group_n - 1L) / 2L),
      norm_factor = if_else(norm_factor <= 0L, 1L, norm_factor),
      avg_n_seg_raw = avg_n_seg,
      avg_n_seg = avg_n_seg_raw / norm_factor,
      norm_mode = "pop_size",
      norm_unit = "n_pairs"
    ) |>
    select(chromosome, pop_id, window_start, window_end, covered_bp, avg_n_seg, avg_n_seg_raw, norm_factor, group_n, norm_mode, norm_unit)
} else {
  d_win <- d_win |>
    mutate(
      avg_n_seg_raw = avg_n_seg,
      norm_factor = 1L,
      group_n = 1L,
      norm_mode = "none",
      norm_unit = "none"
    ) |>
    select(chromosome, pop_id, window_start, window_end, covered_bp, avg_n_seg, avg_n_seg_raw, norm_factor, group_n, norm_mode, norm_unit)
}

cat("__ writing output __\n")
write_tsv(d_win, args$out_file)

cat("__ done! __\n")
