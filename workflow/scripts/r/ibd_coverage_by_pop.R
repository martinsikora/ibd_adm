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

## Compute disjoint chromosome segments with overlap counts.
## Input must contain: pop_id, chromosome, pos_start, pos_end.
## Coordinates are treated as half-open [pos_start, pos_end) unless inclusive_end=TRUE.
compute_overlap_coverage <- function(d, inclusive_end = FALSE) {
  events <- d |>
    transmute(
      pop_id,
      chromosome,
      pos = pos_start,
      delta = 1L
    ) |>
    bind_rows(
      d |>
        transmute(
          pop_id,
          chromosome,
          pos = if (inclusive_end) pos_end + 1L else pos_end,
          delta = -1L
        )
    ) |>
    group_by(pop_id, chromosome, pos) |>
    summarise(delta = sum(delta), .groups = "drop") |>
    arrange(pop_id, chromosome, pos)

  events |>
    group_by(pop_id, chromosome) |>
    arrange(pos, .by_group = TRUE) |>
    mutate(
      n_seg = cumsum(delta),
      pos_next = lead(pos)
    ) |>
    filter(
      !is.na(pos_next),
      n_seg > 0L,
      pos_next > pos
    ) |>
    transmute(
      pop_id,
      chromosome,
      pos_start = pos,
      pos_end = pos_next,
      n_seg
    ) |>
    ungroup()
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--ibd_file",
  action = "store",
  dest = "ibd_file",
  help = "Input file with pairwise IBD tracts"
)

parser$add_argument("-s", "--sample_file",
  action = "store",
  dest = "sample_file",
  help = "File with sample_id and pop_id columns"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output filename with overlap coverage per pop_id/chromosome"
)

parser$add_argument("--sample_col",
  action = "store",
  dest = "sample_col",
  default = "sample2",
  choices = c("sample1", "sample2"),
  help = "Which individual of each pair is mapped to pop_id [default %(default)s]"
)

parser$add_argument("--min_l_cm",
  action = "store",
  dest = "min_l_cm",
  type = "double",
  default = 0,
  help = "Minimum tract length in cM [default %(default)s]"
)

parser$add_argument("--min_lod",
  action = "store",
  dest = "min_lod",
  type = "double",
  default = 0,
  help = "Minimum LOD score [default %(default)s]"
)

parser$add_argument("--inclusive_end",
  action = "store_true",
  default = FALSE,
  help = "Treat end coordinate as inclusive (BED-like off by default)"
)

args <- parser$parse_args()


## --------------------------------------------------
## read data

cat("__ reading IBD data __\n")
## columns are read by position (1-6 and 9), not by header name
ibd <- read_tsv(args$ibd_file,
  skip = 1, col_names = FALSE, col_select = c(1:6, 9),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)
colnames(ibd) <- c("sample1", "sample2", "chromosome", "pos_start", "pos_end", "lod", "l_cm")
ibd <- ibd |>
  mutate(across(c(pos_start, pos_end), as.numeric), across(c(lod, l_cm), as.numeric))

cat("__ reading metadata __\n")
sample_map <- read_tsv(args$sample_file, show_col_types = FALSE)

if (!"sample_id" %in% colnames(sample_map) || !"pop_id" %in% colnames(sample_map)) {
  stop("sample_file must include columns: sample_id, pop_id")
}

## --------------------------------------------------
## preprocess

cat("__ processing tracts __\n")

ibd1 <- ibd |>
  filter(l_cm >= args$min_l_cm, lod >= args$min_lod)

if ("group" %in% colnames(sample_map)) {
  valid_samples <- sample_map |>
    filter(group != "exclude", pop_id != "exclude") |>
    pull(sample_id)
} else {
  valid_samples <- sample_map |>
    filter(pop_id != "exclude") |>
    pull(sample_id)
}

d <- ibd1 |>
  filter(.data[[args$sample_col]] %in% valid_samples) |>
  transmute(
    sample_id = .data[[args$sample_col]],
    chromosome = as.character(chromosome),
    pos_start = as.integer(pos_start),
    pos_end = as.integer(pos_end)
  ) |>
  filter(
    !is.na(sample_id),
    !is.na(chromosome),
    !is.na(pos_start),
    !is.na(pos_end),
    pos_end > pos_start
  ) |>
  left_join(
    sample_map |>
      select(sample_id, pop_id),
    by = c("sample_id" = "sample_id")
  ) |>
  filter(!is.na(pop_id), pop_id != "exclude") |>
  select(pop_id, chromosome, pos_start, pos_end)


## --------------------------------------------------
## compute and write output

cat("__ computing overlap coverage __\n")
cov_by_pop <- compute_overlap_coverage(d, inclusive_end = args$inclusive_end) |>
  arrange(pop_id, chromosome, pos_start)

cat("__ writing output __\n")
write_tsv(cov_by_pop, args$out_file)

cat("__ done! __\n")
