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

## Stage 3: hierarchical clustering of the cluster_full individuals. Also
## precomputes the tree-only hierarchy strings (cl_hier / cl_ids_expand) so the
## per-height cut stage is fast. Shared across all cut heights.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(data.tree)
  library(fastcluster)
})

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(sub("^--file=", "", f[1])) else "."
}
source(file.path(get_script_dir(), "cluster_lib.R"))

parser <- ArgumentParser()
parser$add_argument("--in", dest = "in_file",
  help = "Input distance-matrix .rds (m_d)")
parser$add_argument("-s", "--sample_file", dest = "sample_file",
  help = "Sample-to-group mapping table")
parser$add_argument("--out_hc", dest = "out_hc",
  help = "Output hclust .rds (res_hc)")
parser$add_argument("--out_hier", dest = "out_hier",
  help = "Output tree-hierarchy .rds (cl_hier + cl_ids_expand)")
parser$add_argument("--clust_method", dest = "clust_method", type = "character",
  default = "ward.D2", help = "Clustering method for hclust [default %(default)s]")
args <- parser$parse_args()

cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file, col_types = "ccc")
inds <- sample_label |>
  filter(group != "exclude") |>
  pull(sample_id)
inds_cl_full <- sample_label |>
  filter(group == "cluster_full") |>
  pull(sample_id)

cat("__ loading distance matrix __\n")
m_d <- readRDS(args$in_file)

cat("__ hierarchical clustering __\n")
## fastcluster::hclust is a drop-in for stats::hclust (identical ward.D2 result)
## with much lower memory/time on the ~15.7k cluster_full individuals
res_hc <- fastcluster::hclust(
  d = as.dist(m_d[inds_cl_full, inds_cl_full]),
  method = args$clust_method
)
rm(m_d)
gc()

cat("__ building tree hierarchy strings __\n")
tree_hier <- build_tree_hier(res_hc, inds)

cat("__ writing outputs __\n")
dir.create(dirname(args$out_hc), showWarnings = FALSE, recursive = TRUE)
saveRDS(res_hc, args$out_hc)
saveRDS(tree_hier, args$out_hier)

cat("__ done! __\n")
