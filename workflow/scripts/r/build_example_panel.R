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

## Build / refresh a custom aggregation panel from your own cluster metadata.
##
## This is a template for turning a hierarchical-clustering assignment into a
## custom panel (config/panels/<name>/). The compound pop_id is
## cluster_label:cluster_alias, with a matching per-cluster colour map, derived
## from two metadata tables you supply:
##   aggregate.tsv  <- sample_info: sample_id, pop_id (= cluster_label:alias),
##                     group (donor_recipient / recipient, derived from the
##                     cluster type column: cluster_full -> donor_recipient,
##                     everything else -> recipient)
##   color_map.tsv  <- cluster_info: pop_id (= cluster_label:alias), color,
##                     fill, shape
## The two files are written in lock-step so their pop_id sets stay identical.
## Samples with no cluster assignment (not part of the clustering) are dropped.
## The catch-all "unassigned" pseudo cluster (no alias, absent from
## cluster_info) is kept as a single "unassigned" pop of recipients, drawn as
## neutral grey circles.
##
## If you prefer, skip this script entirely and hand-write aggregate.tsv +
## color_map.tsv directly (see docs/CONFIGURATION.md for the schema).

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(tibble)
})

## the "unassigned" pseudo cluster and its neutral grey-circle styling
UNASSIGNED_LABEL <- "unassigned"
UNASSIGNED_COLOR <- "#808080"
UNASSIGNED_SHAPE <- "16" # solid circle


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("--sample_info",
  action = "store",
  dest = "sample_info",
  default = "config/metadata/sample_info.tsv",
  help = "Sample metadata file (per-sample cluster_label, alias, type) [default %(default)s]" # nolint
)

parser$add_argument("--cluster_info",
  action = "store",
  dest = "cluster_info",
  default = "config/metadata/cluster_info.tsv",
  help = "Cluster metadata file (per-cluster alias + colour scheme) [default %(default)s]" # nolint
)

parser$add_argument("--alias_col",
  action = "store",
  dest = "alias_col",
  default = "cluster_alias",
  help = "column used as the alias part of the compound pop_id [default %(default)s]" # nolint
)

parser$add_argument("--group_col",
  action = "store",
  dest = "group_col",
  default = "cluster_assignment",
  help = "sample_info column giving the cluster type; cluster_full -> donor_recipient, else recipient [default %(default)s]" # nolint
)

parser$add_argument("--out_dir",
  action = "store",
  dest = "out_dir",
  default = "config/panels/example_panel",
  help = "Output panel directory [default %(default)s]"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

sample_info <- read_tsv(args$sample_info, col_types = cols(.default = "c"))
cluster_info <- read_tsv(args$cluster_info, col_types = cols(.default = "c"))

req_si <- c("sample_id", "cluster_label", args$alias_col, args$group_col)
if (!all(req_si %in% colnames(sample_info))) {
  stop(sprintf(
    "sample_info must contain columns: %s",
    paste(req_si, collapse = ", ")
  ))
}
req_ci <- c("cluster_label", args$alias_col, "color", "fill", "shape")
if (!all(req_ci %in% colnames(cluster_info))) {
  stop(sprintf(
    "cluster_info must contain columns: %s",
    paste(req_ci, collapse = ", ")
  ))
}

is_blank <- function(x) is.na(x) | x == "" | x == "NA"


## --------------------------------------------------
## build aggregate rows from sample_info

cat("__ building aggregate.tsv __\n")

base <- sample_info |>
  select(
    sample_id,
    cluster_label,
    alias = all_of(args$alias_col),
    ctype = all_of(args$group_col)
  ) |>
  ## keep only clustered samples; drop NA-assignment individuals that are not
  ## part of the clustering at all
  filter(!is_blank(ctype), !is_blank(cluster_label))

## real clusters (must carry an alias) vs the "unassigned" pseudo-cluster
real <- base |> filter(cluster_label != UNASSIGNED_LABEL)
unassigned <- base |> filter(cluster_label == UNASSIGNED_LABEL)

if (any(is_blank(real$alias))) {
  stop(sprintf(
    "%d clustered samples have a missing %s value",
    sum(is_blank(real$alias)), args$alias_col
  ))
}

## the alias must be constant within each cluster, else the compound is ambiguous
inconsistent <- real |>
  group_by(cluster_label) |>
  summarise(n_alias = n_distinct(alias), .groups = "drop") |>
  filter(n_alias > 1)
if (nrow(inconsistent) > 0) {
  stop(sprintf(
    "%s is not constant within %d cluster(s); compound pop_id would be ambiguous",
    args$alias_col, nrow(inconsistent)
  ))
}

## group: donor_recipient for the tree-clustered core (cluster_full), recipient
## for everything else (cluster_min_dist); the unassigned pseudo-cluster is
## folded in as recipients under the bare "unassigned" pop_id
agg <- real |>
  transmute(
    sample_id,
    pop_id = paste(cluster_label, alias, sep = ":"),
    group = ifelse(ctype == "cluster_full", "donor_recipient", "recipient")
  )
if (nrow(unassigned) > 0) {
  agg <- bind_rows(
    agg,
    unassigned |> transmute(sample_id, pop_id = UNASSIGNED_LABEL, group = "recipient")
  )
}


## --------------------------------------------------
## build colour map from cluster_info, restricted to panel clusters

cat("__ building color_map.tsv __\n")

panel_clusters <- unique(real$cluster_label)

missing_color <- setdiff(panel_clusters, cluster_info$cluster_label)
if (length(missing_color) > 0) {
  stop(sprintf(
    "%d panel cluster(s) missing from cluster_info (e.g. %s)",
    length(missing_color), paste(head(missing_color, 5), collapse = ", ")
  ))
}

cm <- cluster_info |>
  filter(cluster_label %in% panel_clusters) |>
  transmute(
    pop_id = paste(cluster_label, .data[[args$alias_col]], sep = ":"),
    color, fill, shape
  ) |>
  distinct(pop_id, .keep_all = TRUE)

## unassigned individuals: single neutral grey filled-circle entry
if (nrow(unassigned) > 0) {
  cm <- bind_rows(
    cm,
    tibble(
      pop_id = UNASSIGNED_LABEL,
      color = UNASSIGNED_COLOR,
      fill = UNASSIGNED_COLOR,
      shape = UNASSIGNED_SHAPE
    )
  )
}


## --------------------------------------------------
## final consistency check: identical pop_id sets

if (!setequal(unique(agg$pop_id), cm$pop_id)) {
  only_agg <- setdiff(unique(agg$pop_id), cm$pop_id)
  only_cm <- setdiff(cm$pop_id, unique(agg$pop_id))
  stop(sprintf(
    "aggregate and color_map pop_id sets differ (%d only in aggregate, %d only in color_map)", # nolint
    length(only_agg), length(only_cm)
  ))
}


## --------------------------------------------------
## write panel files

cat("__ writing panel files __\n")

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(agg, file.path(args$out_dir, "aggregate.tsv"))
write_tsv(cm, file.path(args$out_dir, "color_map.tsv"))

cat(sprintf(
  "wrote %s: %d samples, %d clusters; %s: %d entries\n",
  file.path(args$out_dir, "aggregate.tsv"), nrow(agg), n_distinct(agg$pop_id),
  file.path(args$out_dir, "color_map.tsv"), nrow(cm)
))
cat("group (cluster type) tally:\n")
print(table(agg$group))

cat("__ done! __\n")
