#!/usr/bin/env Rscript
# Copyright 2023 Martin Sikora <martin.sikora@sund.ku.dk>
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
  library(tibble)
})


## --------------------------------------------------
## helpers

silhouette_mean <- function(dist_mat, clusters) {
  n <- length(clusters)
  s <- numeric(n)
  for (i in seq_len(n)) {
    same <- which(clusters == clusters[i] & seq_len(n) != i)
    if (length(same) == 0) {
      a <- 0
    } else {
      a <- mean(dist_mat[i, same])
    }
    other_clusters <- unique(clusters[clusters != clusters[i]])
    if (length(other_clusters) == 0) {
      b <- 0
    } else {
      b <- min(sapply(other_clusters, function(c) {
        idx <- which(clusters == c)
        mean(dist_mat[i, idx])
      }))
    }
    if (max(a, b) == 0) {
      s[i] <- 0
    } else {
      s[i] <- (b - a) / max(a, b)
    }
  }
  mean(s)
}

pick_medoids <- function(dist_mat, clusters) {
  meds <- sapply(sort(unique(clusters)), function(c) {
    idx <- which(clusters == c)
    if (length(idx) == 1) {
      return(idx)
    }
    dsub <- dist_mat[idx, idx, drop = FALSE]
    idx[which.min(rowSums(dsub))]
  })
  as.integer(meds)
}

pick_farthest <- function(dist_mat, k) {
  n <- nrow(dist_mat)
  if (k >= n) {
    return(seq_len(n))
  }

  ## deterministic seed point: global medoid
  selected <- which.min(rowSums(dist_mat))
  selected <- as.integer(selected[1])

  while (length(selected) < k) {
    remaining <- setdiff(seq_len(n), selected)
    min_d <- apply(dist_mat[remaining, selected, drop = FALSE], 1, min)
    best <- remaining[which.max(min_d)]
    selected <- c(selected, best)
  }
  as.integer(selected)
}

topology_leaf_dist <- function(hc) {
  n <- length(hc$order)
  total_nodes <- 2L * n - 1L

  parent <- integer(total_nodes)
  parent[] <- 0L

  to_node_id <- function(x) {
    if (x < 0) {
      return(as.integer(-x))
    }
    as.integer(n + x)
  }

  for (i in seq_len(n - 1L)) {
    node_id <- n + i
    ch <- hc$merge[i, ]
    c1 <- to_node_id(ch[1])
    c2 <- to_node_id(ch[2])
    parent[c1] <- node_id
    parent[c2] <- node_id
  }

  root <- 2L * n - 1L
  depth <- integer(total_nodes)
  for (node in seq_len(total_nodes)) {
    d <- 0L
    cur <- node
    while (cur != root && parent[cur] != 0L) {
      d <- d + 1L
      cur <- parent[cur]
    }
    depth[node] <- d
  }

  leaf_anc <- vector("list", n)
  for (leaf in seq_len(n)) {
    anc <- integer(0)
    cur <- leaf
    while (cur != 0L) {
      anc <- c(anc, cur)
      if (cur == root) break
      cur <- parent[cur]
    }
    leaf_anc[[leaf]] <- anc
  }

  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    ai <- leaf_anc[[i]]
    if (i < n) {
      for (j in seq.int(i + 1L, n)) {
        aj <- leaf_anc[[j]]
        lca <- ai[ai %in% aj][1]
        dij <- depth[i] + depth[j] - 2L * depth[lca]
        d[i, j] <- dij
        d[j, i] <- dij
      }
    }
  }
  d
}

## Select, per TVD cluster, the population sitting furthest from its divergence
## point -- i.e. the tip with the longest terminal (pendant) branch on an NJ
## tree. Long pendant = strong population-specific drift with little shared
## post-divergence history => a well-differentiated, "clean" source surrogate.
## A size guard drops singleton/low-support tips whose long branch may be an
## artifact (e.g. low coverage) rather than genuine drift.
pick_differentiated <- function(
    dist_mat,
    pops,
    clusters,
    sizes = NULL,
    min_size = 1L,
    relative = FALSE
) {
  ## work in a canonical order == pops
  dist_mat <- dist_mat[pops, pops, drop = FALSE]
  clusters <- clusters[pops]

  ## additive tree -> pendant (terminal) branch length per tip. NJ can emit
  ## small negative pendants (well-known artifact); floor them at 0.
  tr <- ape::nj(as.dist(dist_mat))
  ntip <- length(tr$tip.label)
  is_term <- tr$edge[, 2] <= ntip
  pend <- numeric(ntip)
  pend[tr$edge[is_term, 2]] <- tr$edge.length[is_term]
  pend <- pmax(pend, 0)
  names(pend) <- tr$tip.label
  pend <- pend[pops]

  if (relative) {
    ## pendant relative to the root-to-tip path: what fraction of a tip's total
    ## divergence is its own terminal (post-split) branch -- "furthest from the
    ## divergence points *within* its branch" rather than globally longest
    dep <- ape::node.depth.edgelength(tr)
    root_to_tip <- dep[seq_len(ntip)]
    names(root_to_tip) <- tr$tip.label
    root_to_tip <- root_to_tip[pops]
    score <- pend / pmax(root_to_tip, 1e-9)
  } else {
    score <- pend
  }

  ## size guard: undersized source clusters get excluded from being chosen
  if (!is.null(sizes)) {
    sz <- as.integer(sizes[pops])
    sz[is.na(sz)] <- 0L
    score[sz < min_size] <- -Inf
  }

  ## one differentiated representative per cluster (longest-pendant tip); if the
  ## size guard emptied a cluster, fall back to its longest pendant regardless
  sel <- integer(0)
  for (cl in sort(unique(clusters))) {
    idx <- which(clusters == cl)
    s <- score[idx]
    if (all(!is.finite(s))) {
      s <- pend[idx]
    }
    sel <- c(sel, idx[which.max(s)])
  }
  as.integer(sel)
}

## Blended picker: partition the tree into clades balanced by TARGET-sample mass
## (repeatedly split the clade holding the most target samples), so densely
## sampled regions are subdivided into many clades and sparse outgroups stay
## lumped -- the clades span where the targets actually are. Pendant length is
## then used only as the *within-clade tiebreak*, choosing the most-differentiated
## (longest terminal branch) representative of each target-spanning clade. This
## yields sources that are differentiated but proximate to the target diversity,
## rather than a set of distal global outgroups.
pick_differentiated_spread <- function(
    dist_mat,
    pops,
    k,
    sizes = NULL,
    min_size = 1L,
    relative = FALSE
) {
  dist_mat <- dist_mat[pops, pops, drop = FALSE]
  n <- length(pops)
  if (k >= n) {
    return(seq_len(n))
  }

  ## pendant (differentiation) score per tip -- same as pick_differentiated
  tr <- ape::nj(as.dist(dist_mat))
  ntip <- length(tr$tip.label)
  is_term <- tr$edge[, 2] <= ntip
  pend <- numeric(ntip)
  pend[tr$edge[is_term, 2]] <- tr$edge.length[is_term]
  pend <- pmax(pend, 0)
  names(pend) <- tr$tip.label
  pend <- pend[pops]
  if (relative) {
    dep <- ape::node.depth.edgelength(tr)
    r2t <- dep[seq_len(ntip)]
    names(r2t) <- tr$tip.label
    r2t <- r2t[pops]
    score <- pend / pmax(r2t, 1e-9)
  } else {
    score <- pend
  }

  ## target mass per tip (sample count; 1 if unknown)
  w <- rep(1, n)
  names(w) <- pops
  if (!is.null(sizes)) {
    s <- as.integer(sizes[pops])
    s[is.na(s)] <- 0L
    w <- pmax(s, 0)
  }

  ## nested clade structure from the dendrogram
  hc <- hclust(as.dist(dist_mat), method = "average")
  merge <- hc$merge
  leaves_cache <- vector("list", nrow(merge))
  node_leaves <- function(node) {
    if (!is.null(leaves_cache[[node]])) {
      return(leaves_cache[[node]])
    }
    res <- integer(0)
    for (ch in merge[node, ]) {
      if (ch < 0) {
        res <- c(res, -ch)
      } else {
        res <- c(res, node_leaves(ch))
      }
    }
    leaves_cache[[node]] <<- res
    res
  }
  clade_mass <- function(cl) {
    if (cl < 0) w[[-cl]] else sum(w[node_leaves(cl)])
  }

  ## top-down, target-mass-balanced partition into k clades
  clades <- list(nrow(merge)) # start at root
  while (length(clades) < k) {
    splittable <- which(vapply(clades, function(x) x > 0, logical(1)))
    if (length(splittable) == 0) break
    m <- vapply(clades[splittable], clade_mass, numeric(1))
    i <- splittable[which.max(m)]
    ch <- merge[clades[[i]], ]
    clades <- c(clades[-i], list(ch[1], ch[2]))
  }

  ## within each clade, pick the longest-pendant size-qualifying tip
  sel <- integer(0)
  for (cl in clades) {
    idx <- if (cl < 0) (-cl) else node_leaves(cl)
    s <- score[idx]
    if (!is.null(sizes)) {
      sz <- as.integer(sizes[pops[idx]])
      sz[is.na(sz)] <- 0L
      s2 <- s
      s2[sz < min_size] <- -Inf
      if (any(is.finite(s2))) s <- s2
    }
    sel <- c(sel, idx[which.max(s)])
  }
  as.integer(unique(sel))
}

pick_tree_spread <- function(
    dist_mat,
    pops,
    k,
    broad_k = NULL,
    max_per_broad_clade = 1L,
    min_tree_dist_quantile = 0,
    label_prefix_parts = 0L,
    max_per_label_prefix = 0L
) {
  n <- nrow(dist_mat)
  if (k >= n) {
    return(seq_len(n))
  }

  hc <- hclust(as.dist(dist_mat), method = "average")
  tree_dist <- topology_leaf_dist(hc)

  ## deterministic start: farthest pair by topology
  tree_dist_upper <- tree_dist
  diag(tree_dist_upper) <- -Inf
  pair_idx <- which(tree_dist_upper == max(tree_dist_upper), arr.ind = TRUE)
  pair_idx <- pair_idx[1, ]
  selected <- as.integer(unique(c(pair_idx[1], pair_idx[2])))

  while (length(selected) < k) {
    remaining <- setdiff(seq_len(n), selected)
    min_tree_d <- apply(tree_dist[remaining, selected, drop = FALSE], 1, min)
    max_min <- max(min_tree_d)
    candidates <- remaining[min_tree_d == max_min]
    if (length(candidates) == 1L) {
      best <- candidates
    } else {
      mean_tree_d <- rowMeans(tree_dist[candidates, selected, drop = FALSE])
      best <- candidates[which.max(mean_tree_d)]
    }
    selected <- c(selected, best)
  }
  as.integer(selected)
}


## --------------------------------------------------
## command line arguments

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "TVD table (pop_id1, pop_id2, tvd)"
)

parser$add_argument("-s", "--sample_file",
  action = "store",
  dest = "sample_file",
  help = "Sample map with sample_id and pop_id"
)

parser$add_argument("-o", "--out",
  action = "store",
  dest = "out_file",
  help = "Output mixture file"
)

parser$add_argument("--k_min",
  action = "store",
  dest = "k_min",
  type = "integer",
  default = 2,
  help = "Minimum k to evaluate"
)

parser$add_argument("--k_max",
  action = "store",
  dest = "k_max",
  type = "integer",
  default = 10,
  help = "Maximum k to evaluate"
)

parser$add_argument("--seed",
  action = "store",
  dest = "seed",
  type = "integer",
  default = 1,
  help = "Random seed"
)

parser$add_argument("--source_pick_method",
  action = "store",
  dest = "source_pick_method",
  type = "character",
  default = "tree_spread",
  help = "Source picker: tree_spread, farthest, or cluster_medoids [default %(default)s]"
)

parser$add_argument("--source_broad_k",
  action = "store",
  dest = "source_broad_k",
  type = "integer",
  default = 0,
  help = "For tree_spread: broad-clade cut k (<=0 means auto)"
)

parser$add_argument("--source_max_per_broad_clade",
  action = "store",
  dest = "source_max_per_broad_clade",
  type = "integer",
  default = 1,
  help = "For tree_spread: preferred max sources per broad clade before spillover [default %(default)s]"
)

parser$add_argument("--source_min_tree_dist_quantile",
  action = "store",
  dest = "source_min_tree_dist_quantile",
  type = "double",
  default = 0.0,
  help = "For tree_spread: minimum pairwise tree-distance quantile to enforce while possible (0 disables; range 0-1)"
)

parser$add_argument("--source_label_prefix_parts",
  action = "store",
  dest = "source_label_prefix_parts",
  type = "integer",
  default = 0,
  help = "For tree_spread: if pop labels are hierarchical (e.g. Cx_y_z), guard by first N parts (0 disables)"
)

parser$add_argument("--source_max_per_label_prefix",
  action = "store",
  dest = "source_max_per_label_prefix",
  type = "integer",
  default = 0,
  help = "For tree_spread: preferred max sources per label-prefix group (used when source_label_prefix_parts>0)"
)

parser$add_argument("--source_relative_pendant",
  action = "store_true",
  dest = "source_relative_pendant",
  default = FALSE,
  help = "For differentiated: score pendant length relative to root-to-tip depth"
)

parser$add_argument("--source_min_cluster_size",
  action = "store",
  dest = "source_min_cluster_size",
  type = "integer",
  default = 1L,
  help = "For differentiated: minimum #samples a source cluster must have to be picked [default %(default)s]"
)

args <- parser$parse_args()

if (!(args$source_pick_method %in% c("tree_spread", "farthest", "cluster_medoids", "differentiated", "differentiated_spread"))) {
  stop("--source_pick_method must be one of: tree_spread, farthest, cluster_medoids, differentiated, differentiated_spread")
}
if (args$source_pick_method %in% c("differentiated", "differentiated_spread") && !requireNamespace("ape", quietly = TRUE)) {
  stop("--source_pick_method differentiated/differentiated_spread requires the 'ape' package")
}
if (args$source_max_per_broad_clade < 1) {
  stop("--source_max_per_broad_clade must be >= 1")
}
if (args$source_min_tree_dist_quantile < 0 || args$source_min_tree_dist_quantile > 1) {
  stop("--source_min_tree_dist_quantile must be in [0, 1]")
}
if (args$source_label_prefix_parts < 0 || args$source_max_per_label_prefix < 0) {
  stop("--source_label_prefix_parts and --source_max_per_label_prefix must be >= 0")
}


## --------------------------------------------------
## read data

cat("__ reading data __\n")

tvd <- read_tsv(args$in_file,
  col_types = "ccd"
)

sample_map <- read_tsv(args$sample_file,
  show_col_types = FALSE
)
if (!all(c("sample_id", "pop_id") %in% colnames(sample_map))) {
  stop("sample_file must contain columns: sample_id, pop_id")
}
sample_map <- sample_map |> select(sample_id, pop_id)

## #samples per (core) population, used by the differentiated picker's size guard
pop_sizes <- table(sample_map$pop_id)

pops <- sort(unique(c(tvd$pop_id1, tvd$pop_id2)))

m <- tvd |>
  mutate(
    pop_id1 = factor(pop_id1, levels = pops),
    pop_id2 = factor(pop_id2, levels = pops)
  ) |>
  pivot_wider(
    names_from = pop_id2,
    values_from = tvd,
    values_fill = 0
  ) |>
  column_to_rownames("pop_id1") |>
  as.matrix()

m <- (m + t(m)) / 2
dist_mat <- as.matrix(as.dist(m))
n <- nrow(dist_mat)

if (n < 2) {
  cat("__ only one population, using first sample as source __\n")
  src <- sample_map |>
    filter(!is.na(pop_id)) |>
    slice(1) |>
    pull(sample_id)
  out <- sample_map |>
    mutate(group = ifelse(sample_id == src, "source", "target"))
  write_tsv(out, args$out_file)
  quit(status = 0)
}


## --------------------------------------------------
## choose k by silhouette

cat("__ selecting k __\n")

hc <- hclust(as.dist(dist_mat), method = "average")

k_min <- max(2, args$k_min)
k_max <- min(args$k_max, n - 1)
if (k_min > k_max) {
  k_min <- k_max
}
ks <- k_min:k_max

sil_scores <- sapply(ks, function(k) {
  cl <- cutree(hc, k = k)
  silhouette_mean(dist_mat, cl)
})

best_k <- ks[which.max(sil_scores)]
clusters <- cutree(hc, k = best_k)


## --------------------------------------------------
## pick sources (medoids)

cat("__ selecting sources __\n")

if (args$source_pick_method == "cluster_medoids") {
  med_idx <- pick_medoids(dist_mat, clusters)
} else if (args$source_pick_method == "farthest") {
  med_idx <- pick_farthest(dist_mat, best_k)
} else if (args$source_pick_method == "differentiated") {
  med_idx <- pick_differentiated(
    dist_mat,
    pops,
    clusters,
    sizes = pop_sizes,
    min_size = args$source_min_cluster_size,
    relative = args$source_relative_pendant
  )
} else if (args$source_pick_method == "differentiated_spread") {
  med_idx <- pick_differentiated_spread(
    dist_mat,
    pops,
    best_k,
    sizes = pop_sizes,
    min_size = args$source_min_cluster_size,
    relative = args$source_relative_pendant
  )
} else {
  broad_k <- ifelse(args$source_broad_k <= 0, NA_integer_, args$source_broad_k)
  med_idx <- pick_tree_spread(
    dist_mat,
    pops,
    best_k,
    broad_k = broad_k,
    max_per_broad_clade = args$source_max_per_broad_clade,
    min_tree_dist_quantile = args$source_min_tree_dist_quantile,
    label_prefix_parts = args$source_label_prefix_parts,
    max_per_label_prefix = args$source_max_per_label_prefix
  )
}
source_pops <- pops[med_idx]

source_samples <- sample_map |>
  filter(pop_id %in% source_pops) |>
  pull(sample_id)


## --------------------------------------------------
## write output

cat("__ writing output __\n")

out <- sample_map |>
  mutate(group = ifelse(sample_id %in% source_samples, "source", "target")) |>
  select(sample_id, group)

write_tsv(out, args$out_file)

cat("__ done! __\n")
