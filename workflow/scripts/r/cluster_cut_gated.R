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

## Sharing-gated adaptive cut.
##
## Combines a COARSE and a FINE cut of the same tree into one labelling, taking
## the fine labels only where the coarse cluster has enough IBD to support the
## finer distinctions. The gate is a DATA-SUFFICIENCY test, not a split-quality
## test -- see the caveat below.
##
## Motivation. A single cut height cannot suit the whole panel. Going from
## h0.5 to h0.3 on the ho_20260806 nscale tree resolves populations in the
## Americas (Mayan/Zapotec, Wayku/Guarani, Andean coast/highland: 7 balanced
## splits, 0 arbitrary) and in later Europe (18 informative / 4 arbitrary), but
## in Africa and South Asia the same step mostly shreds single sampled
## populations (YRI 1 -> 3 clusters, Yoruba 1 -> 3, Naidu, Vysya, Pulliyar,
## JuHoan), dropping African population cohesion from 0.559 to 0.399.
##
## What separates those cases is how much IBD the cluster carries: the median
## per-sample genome-wide total is ~600k cM for Americas clusters that split
## informatively and ~15-18k cM for African ones that do not. Of the statistics
## tested on a 139-split benchmark, this was the best (AUC 0.775); a permutation
## gap test (0.502), odd/even-genome reproducibility (0.441-0.613) and the
## direct cross-daughter sharing ratio (0.753) all did worse.
##
## CAVEAT, and it matters. Sharing magnitude separates low-sharing REGIONS from
## the rest; it does not distinguish good splits from bad ones *within* them.
## African informative vs arbitrary splits sit at 17.9k vs 14.8k median sharing,
## South Asian ones at 125-188k vs 123-161k -- overlapping. So the gate freezes
## low-sharing clusters wholesale, suppressing ~13 true distinctions
## (Sugali/Adi_Dravider, Tiwari/Bhumihar, ...) along with the spurious ones. It
## buys Americas/Europe resolution at the price of no refinement anywhere the
## data are thin. A metadata-based guard scored better (AUC 0.885) but was
## rejected as it makes the clustering depend on group_label quality.
##
## Non-nesting. The two cuts are not strictly nested -- dynamicTreeCut's PAM
## stage assigns by distance, not purely by topology, so on this tree 16 fine
## clusters (681 cluster_full samples) span more than one coarse cluster. The
## gate is applied per COARSE cluster to its own samples, so such a fine cluster
## can end up represented by only the subset whose coarse cluster passed. Every
## emitted group is still a well-defined set of samples; the count of affected
## fine clusters is reported so the effect stays visible.
##
## Usage:
##   cluster_cut_gated.R --base BASE.clusters.tsv --fine FINE.clusters.tsv \
##     --matrix m_raw.rds --min_sharing 250000 --out OUT.clusters.tsv
##   --base_height / --fine_height select the rows to use (default: read from
##   the file, which must contain exactly one non-negative cut_height).

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
})

parser <- ArgumentParser()
parser$add_argument("--base", dest = "base", help = "Coarse clusters.tsv")
parser$add_argument("--fine", dest = "fine", help = "Fine clusters.tsv")
parser$add_argument("--matrix", dest = "matrix",
  help = "Feature matrix .rds (m_raw); row sums give per-sample total IBD")
parser$add_argument("--min_sharing", dest = "min_sharing", type = "double",
  default = 250000,
  help = "Median per-sample total IBD (cM) a coarse cluster needs before its finer splits are used [default %(default)s]") # nolint
parser$add_argument("--base_height", dest = "base_height", default = NULL)
parser$add_argument("--fine_height", dest = "fine_height", default = NULL)
## NOTE: there is no minimum-daughter-size option here. Tiny
## (1-2 member) sub-clusters are a problem, but they are handled by
## dynamicTreeCut's own minClusterSize -- `clustering.cl_size` in config.yml,
## passed to cluster_cut.R. Raising it from 2 to 3 removes every sub-3 cluster
## at the CUT stage (151 -> 0 at h0.5), which is strictly better than reverting
## them here: a post-hoc revert only relocates the problem, because the rest of
## the coarse cluster has already taken fine labels, leaving the coarse label
## holding the 1-2 reverted stragglers. Measured: such a guard reverted 92
## samples and 28 landed in a newly-tiny group.
parser$add_argument("--out", dest = "out", help = "Output clusters.tsv")
args <- parser$parse_args()

pick_height <- function(df, want) {
  hs <- unique(df$cut_height)
  hs <- hs[as.numeric(hs) >= 0]
  if (!is.null(want)) {
    sel <- hs[abs(as.numeric(hs) - as.numeric(want)) < 1e-9]
    if (length(sel) != 1) stop("cut height ", want, " not found; have: ", paste(hs, collapse = ", "))
    return(sel)
  }
  if (length(hs) != 1) {
    stop("file has several cut heights (", paste(hs, collapse = ", "),
         "); pass --base_height / --fine_height")
  }
  hs
}

cat("__ reading cuts __\n")
b_all <- read_tsv(args$base, col_types = cols(.default = "c"))
f_all <- read_tsv(args$fine, col_types = cols(.default = "c"))
bh <- pick_height(b_all, args$base_height)
fh <- pick_height(f_all, args$fine_height)
b <- b_all |> filter(cut_height == bh)
f <- f_all |> filter(cut_height == fh)
cat(sprintf("   coarse h=%s: %d rows, %d clusters | fine h=%s: %d rows, %d clusters\n",
  bh, nrow(b), n_distinct(b$cluster_label), fh, nrow(f), n_distinct(f$cluster_label)))
if (!setequal(b$sample_id, f$sample_id)) stop("the two cuts cover different samples")

cat("__ per-sample total IBD __\n")
m <- readRDS(args$matrix)
tot <- rowSums(m)
rm(m); invisible(gc())

## coarse cluster -> median total IBD over its cluster_full members (fall back to
## all members if a cluster has none), then the gate decision per coarse cluster
bf <- b |> filter(group == "cluster_full")
med_of <- function(df) {
  df |>
    mutate(t = unname(tot[sample_id])) |>
    filter(!is.na(t)) |>
    group_by(cluster_label) |>
    summarise(med = median(t), n = n(), .groups = "drop")
}
med <- med_of(bf)
missing <- setdiff(unique(b$cluster_label), med$cluster_label)
if (length(missing) > 0) med <- bind_rows(med, med_of(b |> filter(cluster_label %in% missing)))
pass <- med$cluster_label[med$med >= args$min_sharing]
pass <- setdiff(pass, "unassigned")   # never subdivide the catch-all bin
cat(sprintf("   %d of %d coarse clusters pass the %.0f cM floor (%d of %d samples)\n",
  length(pass), nrow(med), args$min_sharing,
  sum(b$cluster_label %in% pass), nrow(b)))

## emit: fine row where the sample's coarse cluster passed, coarse row otherwise
take_fine <- b$sample_id[b$cluster_label %in% pass]

out <- bind_rows(
  f |> filter(sample_id %in% take_fine),
  b |> filter(!sample_id %in% take_fine)
) |>
  ## Stamp a single, DISTINCT cut_height for the gated labelling, e.g. "0.5g0.2".
  ## It has to be one value because every consumer filters on it, and it has to
  ## differ from the base height so this panel is addressable separately from the
  ## plain base-height panel in the workflow's {height} wildcards. Consumers
  ## compare it as a string when it does not parse as a number
  ## (make_agg_panel_from_clusters.py: height_match; build_cluster_assignment.py:
  ## --cut), so a non-numeric tag is safe.
  mutate(cut_height = paste0(bh, "g", fh)) |>
  arrange(match(sample_id, b$sample_id))

## how many fine clusters were partially taken (the non-nesting effect)
split_fine <- f |>
  mutate(taken = sample_id %in% take_fine) |>
  group_by(cluster_label) |>
  summarise(k = n_distinct(taken), .groups = "drop") |>
  filter(k > 1) |>
  nrow()
cat(sprintf("   fine clusters only partially taken (non-nesting): %d\n", split_fine))

cat("__ writing __\n")
write_tsv(out, args$out)
cat(sprintf("wrote %s: %d rows, %d clusters (coarse %d, fine %d)\n",
  args$out, nrow(out), n_distinct(out$cluster_label),
  n_distinct(b$cluster_label), n_distinct(f$cluster_label)))
cat("__ done! __\n")
