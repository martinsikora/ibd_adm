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

## Stage 2: apply the feature transforms and compute the pairwise distance
## matrix (m_d). Depends only on the distance config (dist_method + transforms),
## so it is shared across all cut heights and clustering methods.

suppressPackageStartupMessages({
  library(argparse)
  library(parallelDist)
})

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(sub("^--file=", "", f[1])) else "."
}
source(file.path(get_script_dir(), "cluster_lib.R"))

parser <- ArgumentParser()
parser$add_argument("--in", dest = "in_file",
  help = "Input feature-matrix .rds (m_raw)")
parser$add_argument("--out", dest = "out_file",
  help = "Output distance-matrix .rds (m_d)")
parser$add_argument("--dist_method", dest = "dist_method", type = "character",
  default = "cosine", help = "Distance method for parDist [default %(default)s]")
parser$add_argument("--normalize_ibd_vectors", dest = "normalize_ibd_vectors",
  action = "store_true", default = FALSE,
  help = "L2-normalize per-sample IBD vectors before distance calculation")
parser$add_argument("--standardize_features", dest = "standardize_features",
  action = "store_true", default = FALSE,
  help = "Z-score standardize each feature before distance calculation")
parser$add_argument("--scale_features", dest = "scale_features",
  action = "store_true", default = FALSE,
  help = "Scale each feature by its SD WITHOUT centring (safe under cosine)")
parser$add_argument("-t", "--threads", dest = "threads", type = "integer",
  default = 1L, help = "Number of threads [default %(default)s]")
args <- parser$parse_args()

cat("__ loading feature matrix __\n")
m_use <- readRDS(args$in_file)

cat("__ applying feature transforms __\n")
m_use <- apply_feature_transforms(
  m_use,
  standardize_features = args$standardize_features,
  normalize_ibd_vectors = args$normalize_ibd_vectors,
  scale_features = args$scale_features
)

cat("__ computing distance matrix __\n")
m_d <- parDist(m_use, method = args$dist_method, threads = args$threads) |>
  as.matrix()
rm(m_use)
gc()

cat("__ writing distance matrix __\n")
dir.create(dirname(args$out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(m_d, args$out_file)

cat("__ done! __\n")
