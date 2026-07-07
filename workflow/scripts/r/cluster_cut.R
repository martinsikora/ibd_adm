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

## Stage 4: adaptive tree cut at a single height, map dynamic clusters onto the
## dendrogram, assign readable labels, add cluster_min_dist individuals, and
## write the clusters.tsv. This is the only stage that depends on the cut
## parameters (height / cl_size / deep_split).

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(dynamicTreeCut)
})

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(sub("^--file=", "", f[1])) else "."
}
source(file.path(get_script_dir(), "cluster_lib.R"))

parser <- ArgumentParser()
parser$add_argument("--hc", dest = "hc_file", help = "hclust .rds (res_hc)")
parser$add_argument("--hier", dest = "hier_file",
  help = "tree-hierarchy .rds (cl_hier + cl_ids_expand)")
parser$add_argument("--dist", dest = "dist_file",
  help = "distance-matrix .rds (m_d)")
parser$add_argument("-s", "--sample_file", dest = "sample_file",
  help = "Sample-to-group mapping table")
parser$add_argument("--out_tsv", dest = "out_file_tsv", help = "Output clusters.tsv")
parser$add_argument("--height", dest = "height", type = "double",
  help = "Cut height for adaptive tree cutting")
parser$add_argument("--cl_size", dest = "cl_size", type = "integer", default = 2L,
  help = "Minimum cluster size for adaptive tree cutting [default %(default)s]")
parser$add_argument("--deep_split", dest = "deep_split", type = "integer",
  default = 3L, help = "Deep split parameter for adaptive tree cutting [default %(default)s]")
parser$add_argument("--knn", dest = "knn", type = "integer", default = 1L,
  help = paste("k for k-NN majority-vote assignment of cluster_min_dist samples;",
    "1 = original single-nearest-neighbour rule [default %(default)s]"))
args <- parser$parse_args()

cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file, col_types = "ccc")
inds_cl_full <- sample_label |>
  filter(group == "cluster_full") |>
  pull(sample_id)

cat("__ loading clustering artifacts __\n")
res_hc <- readRDS(args$hc_file)
tree_hier <- readRDS(args$hier_file)
cl_hier <- tree_hier$cl_hier
cl_ids_expand <- tree_hier$cl_ids_expand
m_d <- readRDS(args$dist_file)

cat("__ adaptive tree cut __\n")
distM_full <- m_d[inds_cl_full, inds_cl_full]
r1 <- cutreeDynamic(res_hc,
  distM = distM_full,
  method = "hybrid",
  deepSplit = args$deep_split,
  minClusterSize = args$cl_size,
  cutHeight = args$height
)
res_hc_cut <- tibble(
  sample_id = inds_cl_full,
  cluster_terminal = r1,
  cut_height = args$height
) |>
  arrange(cluster_terminal, sample_id)
rm(distM_full)
gc()

## --------------------------------------------------
## reformat: map each dynamic cluster onto its plurality dendrogram node
cl_expand <- cl_hier |>
  left_join(res_hc_cut, by = "sample_id") |>
  left_join(cl_ids_expand, by = "cluster_id", relationship = "many-to-many") |>
  arrange(cut_height, cluster_id, cluster_level)

cl_term <- cl_expand |>
  filter(cluster_terminal != 0) |>
  count(cut_height, cluster_terminal, cluster_level, cluster_id_anc) |>
  group_by(cut_height, cluster_terminal) |>
  slice_max(n) |>
  slice_max(cluster_level)

cl_1 <- cl_expand |>
  semi_join(cl_term) |>
  select(sample_id, cut_height, cluster_id_anc, cluster_level) |>
  rename("cluster_id" = "cluster_id_anc")

cl_2 <- cl_hier |>
  left_join(cl_ids_expand, by = "cluster_id", relationship = "many-to-many") |>
  filter(cluster_id == cluster_id_anc) |>
  mutate(cut_height = -1, cluster_id) |>
  select(sample_id, cut_height, cluster_id, cluster_level)

## samples left unassigned by cutreeDynamic (cluster_terminal == 0) have no
## corresponding dendrogram node; label them explicitly instead of collapsing
## them onto the tree root (which previously produced a spurious diverse "C0")
cl_unassigned <- res_hc_cut |>
  filter(cluster_terminal == 0) |>
  transmute(
    sample_id,
    cut_height,
    cluster_id = "unassigned",
    cluster_level = NA_integer_
  )

cl_final <- bind_rows(cl_2, cl_1, cl_unassigned) |>
  left_join(sample_label, by = "sample_id") |>
  arrange(sample_id)

## add min-dist clustering individuals to the nearest cluster.
## Assignment uses a k-nearest-neighbour majority vote over cluster_full samples
## (k = args$knn); k = 1 reproduces the original single-nearest-neighbour rule.
## Voting on the cut-level cluster prevents a single mis-placed "bridge" sample in
## an isolated, low-sharing cluster from capturing many unrelated recipients (the
## min_dist "sink" effect): a poorly-characterised sample now has to be close to
## a plurality of a cluster's members, not just one accidental nearest neighbour.
inds_cl_min <- sample_label |>
  filter(group == "cluster_min_dist") |>
  pull(sample_id)

if (length(inds_cl_min) > 0) {
  ## cut-level cluster of each cluster_full sample (one row per sample)
  cut_cluster_by_full <- cl_final |>
    filter(cut_height != -1) |>
    distinct(sample_id, cluster_id) |>
    (\(d) setNames(d$cluster_id, d$sample_id))()
  k <- min(args$knn, length(inds_cl_full))

  cl_final_inds_cl_min <- map_dfr(inds_cl_min, ~ {
    ## k nearest cluster_full samples, in ascending distance
    neigh <- inds_cl_full[order(m_d[.x, inds_cl_full])[seq_len(k)]]
    neigh_cl <- cut_cluster_by_full[neigh]
    tab <- table(neigh_cl)
    winners <- names(tab)[tab == max(tab)]
    ## tie-break: winning cluster of the closest neighbour among the tied set
    win <- neigh_cl[neigh_cl %in% winners][1]
    ## representative = closest neighbour in the winning cluster; the min-dist
    ## sample inherits its (deep + cut-level) cluster assignment rows
    rep_id <- neigh[which(neigh_cl == win)[1]]
    cl_final |>
      filter(sample_id == rep_id) |>
      select(sample_id:cluster_level) |>
      mutate(sample_id = .x) |>
      left_join(sample_label, by = "sample_id")
  })
  cl_final <- bind_rows(cl_final, cl_final_inds_cl_min)
}
rm(m_d)
gc()

## add hierarchy-aware readable cluster labels + parent labels
cl_final <- add_cluster_labels(cl_final)

cat("__ writing output __\n")
dir.create(dirname(args$out_file_tsv), showWarnings = FALSE, recursive = TRUE)
write_tsv(cl_final, file = args$out_file_tsv)

cat("__ done! __\n")
