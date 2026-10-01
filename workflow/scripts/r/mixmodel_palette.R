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

## Raw-cM source matrix from the summed palettes of the source individuals: mean per individual in each
## source population (n_src, named by population). A source individual cannot share with itself, so its
## own-cluster entry is built from n - 1 donors; that entry is rescaled by n / (n - 1) to the n donors a
## target sees (donor_n, named by population).
raw_source_matrix <- function(sum_mat, n_src, donor_n) {
  m <- sweep(sum_mat, 2, n_src[colnames(sum_mat)], "/")
  for (k in colnames(m)) {
    r <- match(k, rownames(m))
    n <- donor_n[k]
    if (!is.na(r) && !is.na(n) && n > 1) m[r, k] <- m[r, k] * n / (n - 1)
  }
  m
}
