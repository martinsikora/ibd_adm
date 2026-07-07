#!/usr/bin/env Rscript
# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>
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

## Stage 1: build the dense symmetric IBD feature matrix (m_raw) from the
## per-chromosome IBD tables. Output is shared across every distance config and
## cut height.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(data.table)
})

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(sub("^--file=", "", f[1])) else "."
}
source(file.path(get_script_dir(), "cluster_lib.R"))

parser <- ArgumentParser()
parser$add_argument("files", nargs = "+",
  help = "Per-chromosome IBD sharing tables")
parser$add_argument("-s", "--sample_file", dest = "sample_file",
  help = "Sample-to-group mapping table")
parser$add_argument("--out", dest = "out_file",
  help = "Output feature-matrix .rds")
parser$add_argument("-t", "--threads", dest = "threads", type = "integer",
  default = 1L, help = "Number of threads [default %(default)s]")
args <- parser$parse_args()

cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file, col_types = "ccc")
inds <- sample_label |>
  filter(group != "exclude") |>
  pull(sample_id)

cat("__ reading IBD data + building feature matrix __\n")
m_raw <- read_ibd_matrix(args$files, inds, threads = args$threads)

cat("__ writing feature matrix __\n")
dir.create(dirname(args$out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(m_raw, args$out_file)

cat("__ done! __\n")
