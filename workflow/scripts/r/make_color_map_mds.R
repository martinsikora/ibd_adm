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
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(Rtsne)
})


## --------------------------------------------------
## helpers

scale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(r[1], r[2]))) {
    return(rep(0.5, length(x)))
  }
  (x - r[1]) / (r[2] - r[1])
}

spread_hue_by_range <- function(x, hue_start = 10, hue_span = 330) {
  h <- x
  n <- length(h)
  if (n <= 1) {
    return(rep(hue_start + hue_span / 2, n))
  }
  hs <- sort(h)
  hs2 <- c(hs, hs[1] + 360)
  gaps <- diff(hs2)
  i_gap <- which.max(gaps)[1]
  start <- hs2[i_gap + 1]
  unwrapped <- (h - start + 360) %% 360
  span <- 360 - gaps[i_gap]
  if (span <= 1e-9) {
    return(rep(hue_start + hue_span / 2, n))
  }
  hue_start + hue_span * (unwrapped / span)
}

embedding_to_pca_axes <- function(emb) {
  if (is.null(dim(emb))) {
    emb <- matrix(emb, ncol = 1)
  }
  emb <- as.matrix(emb)
  rot <- stats::prcomp(emb, center = TRUE, scale. = FALSE)$x
  rot <- rot[, seq_len(min(3, ncol(rot))), drop = FALSE]
  if (ncol(rot) < 3) {
    rot <- cbind(rot, matrix(0, nrow = nrow(rot), ncol = 3 - ncol(rot)))
  }
  colnames(rot) <- c("pc1", "pc2", "pc3")
  rot
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "TVD matrix in long format (pop_id1, pop_id2, tvd)"
)

parser$add_argument("-o", "--out",
  action = "store",
  dest = "out_file",
  help = "Output color map file"
)

parser$add_argument("-k", "--k_clusters",
  action = "store",
  dest = "k_clusters",
  type = "integer",
  default = 8,
  help = "Number of shape clusters (default: 8)"
)

parser$add_argument("-s", "--shapes",
  action = "store",
  dest = "shapes",
  default = NULL,
  help = "Comma-separated shape numbers (default: 0:15)"
)

parser$add_argument("-n", "--k_neighbors",
  action = "store",
  dest = "k_neighbors",
  type = "integer",
  default = 5,
  help = "Number of nearest neighbors to differentiate by shape (default: 5)"
)

parser$add_argument("--chroma_min",
  action = "store",
  dest = "chroma_min",
  type = "double",
  default = 20,
  help = "Minimum chroma for HCL colors (default: 20)"
)

parser$add_argument("--chroma_max",
  action = "store",
  dest = "chroma_max",
  type = "double",
  default = 130,
  help = "Maximum chroma for HCL colors (default: 130)"
)

parser$add_argument("--lum_min",
  action = "store",
  dest = "lum_min",
  type = "double",
  default = 15,
  help = "Minimum luminance for HCL colors (default: 15)"
)

parser$add_argument("--lum_max",
  action = "store",
  dest = "lum_max",
  type = "double",
  default = 95,
  help = "Maximum luminance for HCL colors (default: 95)"
)

parser$add_argument("--gamma_c",
  action = "store",
  dest = "gamma_c",
  type = "double",
  default = 0.7,
  help = "Gamma exponent for chroma scaling (default: 0.7)"
)

parser$add_argument("--gamma_l",
  action = "store",
  dest = "gamma_l",
  type = "double",
  default = 0.8,
  help = "Gamma exponent for luminance scaling (default: 0.8)"
)

parser$add_argument("--hue_scale",
  action = "store",
  dest = "hue_scale",
  type = "double",
  default = 1.15,
  help = "Multiplicative hue scaling factor (default: 1.15)"
)

parser$add_argument("--hue_rotate",
  action = "store",
  dest = "hue_rotate",
  type = "double",
  default = 25,
  help = "Additive hue rotation in degrees (default: 25)"
)

parser$add_argument("--hue_spread_mode",
  action = "store",
  dest = "hue_spread_mode",
  type = "character",
  default = "range",
  help = "Hue spreading mode: raw, range, rank (default: range)"
)

parser$add_argument("--lc_spread_mode",
  action = "store",
  dest = "lc_spread_mode",
  type = "character",
  default = "raw",
  help = paste(
    "Chroma/luminance spreading: raw (min-max of the embedding axes) or rank",
    "(quantile-spread, the analogue of hue_spread_mode rank). raw is sensitive",
    "to a skewed embedding -- a few extreme values on the luminance axis squash",
    "every other cluster into the dark end, which kills the hues that need high",
    "luminance (yellow) or high chroma (red). (default: raw)"
  )
)

parser$add_argument("--embedding",
  action = "store",
  dest = "embedding",
  type = "character",
  default = "tsne3",
  help = "Embedding method: tsne3, mds3 (default: tsne3)"
)

parser$add_argument("--mapping",
  action = "store",
  dest = "mapping",
  type = "character",
  default = "radial",
  help = "Color mapping mode: radial, pca_axes (default: radial)"
)

args <- parser$parse_args()
if (args$chroma_max <= args$chroma_min) {
  stop("--chroma_max must be > --chroma_min")
}
if (args$lum_max <= args$lum_min) {
  stop("--lum_max must be > --lum_min")
}
if (args$gamma_c <= 0 || args$gamma_l <= 0) {
  stop("--gamma_c and --gamma_l must be > 0")
}
if (!(args$hue_spread_mode %in% c("raw", "range", "rank"))) {
  stop("--hue_spread_mode must be one of: raw, range, rank")
}
if (!(args$lc_spread_mode %in% c("raw", "rank"))) {
  stop("--lc_spread_mode must be one of: raw, rank")
}
if (!(args$embedding %in% c("tsne3", "mds3"))) {
  stop("--embedding must be one of: tsne3, mds3")
}
if (!(args$mapping %in% c("radial", "pca_axes"))) {
  stop("--mapping must be one of: radial, pca_axes")
}


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

tvd <- read_tsv(args$in_file,
  col_types = "ccd"
)

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


## --------------------------------------------------
## embedding

cat("__ computing embedding __\n")

dist_mat <- as.matrix(as.dist(m))
n_samples <- nrow(dist_mat)
perp <- min(30, max(2, floor((n_samples - 1) / 3)))

set.seed(1)
if (args$embedding == "tsne3") {
  emb <- Rtsne(dist_mat,
    is_distance = TRUE,
    dims = 3,
    perplexity = perp,
    pca = FALSE,
    verbose = FALSE
  )$Y
} else {
  emb <- stats::cmdscale(stats::as.dist(m), k = 3, eig = FALSE)
}

coords <- as_tibble(emb, .name_repair = "minimal")
colnames(coords) <- c("x", "y", "z")
coords$pop_id <- rownames(dist_mat)


## --------------------------------------------------
## colors in HCL

cat("__ mapping colors __\n")

if (args$mapping == "pca_axes") {
  rot_coords <- embedding_to_pca_axes(emb)
  coords <- as_tibble(rot_coords, .name_repair = "minimal")
  colnames(coords) <- c("pc1", "pc2", "pc3")
  coords$pop_id <- rownames(dist_mat)

  h <- 360 * scale01(coords$pc1)
  if (args$hue_spread_mode == "range") {
    h <- spread_hue_by_range(coords$pc1, hue_start = 10, hue_span = 330)
  } else if (args$hue_spread_mode == "rank") {
    rnk <- rank(coords$pc1, ties.method = "average")
    h <- 10 + 330 * ((rnk - 1) / max(1, length(h) - 1))
  }
  c_norm <- scale01(coords$pc2)
  l_norm <- scale01(coords$pc3)
} else {
  h <- (atan2(coords$y, coords$x) * 180 / pi) %% 360
  if (args$hue_spread_mode == "range") {
    h <- spread_hue_by_range(h, hue_start = 10, hue_span = 330)
  } else if (args$hue_spread_mode == "rank") {
    rnk <- rank(h, ties.method = "average")
    h <- 10 + 330 * ((rnk - 1) / max(1, length(h) - 1))
  }
  r <- sqrt(coords$x^2 + coords$y^2)
  c_norm <- scale01(r)
  l_norm <- scale01(coords$z)
}
h <- (h * args$hue_scale + args$hue_rotate) %% 360

## Quantile-spread chroma and luminance. scale01() is min-max, so a skewed
## embedding axis (a few extreme values) leaves most clusters bunched at the
## dark, desaturated end -- on the h0.5 panel that dropped the median luminance
## from 68 to 54 and left 0.4% vivid yellows and no vivid reds, even though the
## red and yellow HUE bins were as populated as before. Ranking uses the full
## chroma/luminance range whatever the axis distribution, while preserving the
## ordering those axes encode.
if (args$lc_spread_mode == "rank") {
  rank01 <- function(x) {
    n <- length(x)
    if (n <= 1) {
      return(rep(0.5, n))
    }
    (rank(x, ties.method = "average") - 1) / (n - 1)
  }
  c_norm <- rank01(c_norm)
  l_norm <- rank01(l_norm)
}

c_norm <- c_norm^args$gamma_c
l_norm <- l_norm^args$gamma_l

chroma <- args$chroma_min + (args$chroma_max - args$chroma_min) * c_norm
lum <- args$lum_min + (args$lum_max - args$lum_min) * l_norm

color <- grDevices::hcl(h, chroma, lum)


## --------------------------------------------------
## shape assignment by nearest-neighbor graph coloring

cat("__ assigning shapes __\n")

if (is.null(args$shapes) || args$shapes == "") {
  shape_vals <- 0:15
} else {
  shape_vals <- as.integer(strsplit(args$shapes, ",")[[1]])
}

n <- nrow(dist_mat)
k_neighbors <- min(args$k_neighbors, n - 1)

neighbors <- lapply(seq_len(n), function(i) {
  if (k_neighbors <= 0) {
    return(integer(0))
  }
  ord <- order(dist_mat[i, ])
  ord[ord != i][seq_len(k_neighbors)]
})

deg <- sapply(neighbors, length)
order_nodes <- order(deg, decreasing = TRUE)
node_shapes <- rep(NA_integer_, n)

for (i in order_nodes) {
  used <- node_shapes[neighbors[[i]]]
  used <- used[!is.na(used)]
  global_counts <- sapply(shape_vals, function(s) sum(node_shapes == s, na.rm = TRUE))
  avail <- setdiff(shape_vals, used)
  if (length(avail) > 0) {
    avail_counts <- global_counts[match(avail, shape_vals)]
    node_shapes[i] <- avail[which.min(avail_counts)]
  } else {
    used_counts <- global_counts[match(shape_vals, shape_vals)]
    node_shapes[i] <- shape_vals[which.min(used_counts)]
  }
}

shape <- node_shapes


## --------------------------------------------------
## write output

cat("__ writing output __\n")

out <- tibble(
  pop_id = coords$pop_id,
  color = color,
  fill = color,
  shape = shape
)

write_tsv(out, args$out_file)

cat("__ done! __\n")
