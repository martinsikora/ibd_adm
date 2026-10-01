## Palette scale helpers for mixmodel_ibd.R (--palette_scale raw|normalized).
##
## normalized: every palette (target individual or source population) is divided by its own total, so
##   sources with a larger total IBD per individual are over-credited (a mixture of normalised palettes
##   weights source k by its ancestry share times its total IBD, not by the share alone).
## raw: sources are mean per-individual palettes in cM, the target palette is fitted up to a free overall
##   scale, and the weights are normalised afterwards, so they are ancestry fractions.

## Model prediction on the proportion scale. Normalised sources mix linearly; raw-cM sources are rescaled
## to sum to 1.
model_pred <- function(source_mat, p, raw = FALSE) {
  q <- as.vector(source_mat %*% p)
  if (raw) q / sum(q) else q
}

## Raw-cM source matrix from the summed palettes of the source individuals: the mean per-individual palette of
## each source population (n_src, named by population). No donor-count correction is needed for palettes from
## aggregate_ibd.py: it sums every individual's within-cluster entry over the n - 1 other members and its
## between-cluster entries over a random subset of n - 1 of the n donors, so a source individual and a target
## are compared on the same number of donors.
raw_source_matrix <- function(sum_mat, n_src) {
  sweep(sum_mat, 2, n_src[colnames(sum_mat)], "/")
}
