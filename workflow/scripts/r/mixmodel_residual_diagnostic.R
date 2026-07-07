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

## Post-hoc residual diagnostic for the IBD mixture models. Reuses an existing
## model output (no re-fit): reconstructs the target sharing vector y and the
## source palette S exactly as mixmodel_ibd.R does (sum IBD over chromosomes,
## normalize each column to sum 1 over the donor palette all_pops), forms the
## residual r = S %*% p - y per target CLUSTER, and projects the unexplained
## part (y - pred) onto every *unused* population to find what a deep/distal
## source is standing in for. Flags sources that behave as bad deep proxies:
## heavy loading + low residual (absorber) + structured leftover pointing at an
## unused population + large nnls/bayesian discordance.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

## ---- paths (positional args) ----
a <- commandArgs(trailingOnly = TRUE)
POP  <- a[1]   # pop_prof.tsv : pop_id1, pop_id2, sum   (all samples -> palette P)
SRCP <- a[2]   # src_prof.tsv : pop_id1, pop_id2, sum   (source samples -> S)
VALP <- a[3]   # val_prof.tsv : sample1, pop_id2, sum    (validation samples)
SMAP <- a[4]   # sample_map.tsv (sample_id, pop_id)
NNLS <- a[5]   # model nnls tsv
BAYE <- a[6]   # model bayesian tsv
VALI <- a[7]   # val.ids (validation sample ids)
OUT  <- a[8]   # output prefix

## ---- donor palette all_pops (base pop_id, strip _r) ----
smap <- read_tsv(SMAP, col_types = cols(.default = "c"))
all_pops <- smap |> mutate(b = sub("_r$", "", pop_id)) |> distinct(b) |> pull(b) |> sort()

## ---- helpers: long (g,pop2,sum) -> normalized matrix rows=all_pops, cols=g ----
build_mat <- function(df, gcol) {
  wide <- df |>
    filter(pop2 %in% all_pops) |>
    group_by(.data[[gcol]], pop2) |> summarise(v = sum(v), .groups = "drop")
  gs <- sort(unique(wide[[gcol]]))
  M <- matrix(0, nrow = length(all_pops), ncol = length(gs),
              dimnames = list(all_pops, gs))
  M[cbind(match(wide$pop2, all_pops), match(wide[[gcol]], gs))] <- wide$v
  cs <- colSums(M); cs[cs == 0] <- 1
  sweep(M, 2, cs, "/")
}

cat("__ reading profiles __\n")
pp <- read_tsv(POP,  col_names = c("g", "pop2", "v"), col_types = "ccd")
sp <- read_tsv(SRCP, col_names = c("g", "pop2", "v"), col_types = "ccd")
vp <- read_tsv(VALP, col_names = c("g", "pop2", "v"), col_types = "ccd")

P <- build_mat(pp, "g")     # palette: every pop's profile (candidates + target y)
S <- build_mat(sp, "g")     # source_mat: source pops' profiles (from source samples)
V <- build_mat(vp, "g")     # validation samples' profiles
src_pops <- colnames(S)

## ---- model outputs -> cluster-mean p (non-_r targets), and per-sample for validation ----
read_p <- function(f) read_tsv(f, col_types = cols_only(
  sample_id = "c", pop_id = "c", group = "c", source_pop = "c", p = "d", res_norm = "d"))
mn <- read_p(NNLS); mb <- read_p(BAYE)

cluster_p <- function(m) m |>
  filter(!grepl("_r$", pop_id)) |>
  group_by(pop_id, source_pop) |> summarise(p = mean(p), .groups = "drop")
cpn <- cluster_p(mn); cpb <- cluster_p(mb)
resn_cl <- mn |> filter(!grepl("_r$", pop_id)) |> distinct(sample_id, pop_id, res_norm) |>
  group_by(pop_id) |> summarise(model_res_norm = mean(res_norm), n = n(), .groups = "drop")

## wide cluster x source p matrices aligned to src_pops
to_wide <- function(cp) {
  w <- cp |> filter(source_pop %in% src_pops) |>
    pivot_wider(names_from = source_pop, values_from = p, values_fill = 0)
  cl <- w$pop_id; W <- as.matrix(w[, -1]); rownames(W) <- cl
  miss <- setdiff(src_pops, colnames(W)); if (length(miss)) W <- cbind(W, matrix(0, nrow(W), length(miss), dimnames = list(NULL, miss)))
  W[, src_pops, drop = FALSE]
}
Pn <- to_wide(cpn); Pb <- to_wide(cpb)

## ================= VALIDATION: recompute res_norm for a few samples =================
val_ids <- readLines(VALI)
cat("\n=== VALIDATION: reconstructed vs model res_norm (nnls) ===\n")
for (s in val_ids) {
  if (!s %in% colnames(V)) next
  y <- V[, s]
  pv <- mn |> filter(sample_id == s, source_pop %in% src_pops)
  p <- setNames(rep(0, length(src_pops)), src_pops); p[pv$source_pop] <- pv$p
  pred <- as.vector(S %*% p)
  rn <- sqrt(mean((pred - y)^2))
  cat(sprintf("  %-12s recon=%.5f  model=%.5f  (p_sum=%.3f)\n",
              s, rn, unique(pv$res_norm)[1], sum(p)))
}

## ================= CLUSTER RESIDUALS + WHAT'S-MISSING PROJECTION =================
## use bayesian cluster p (sparse; the estimator that exposes proxy behaviour)
targets <- intersect(rownames(Pb), colnames(P))
unused  <- setdiff(colnames(P), src_pops)      # candidate "missing" pops
Pu <- P[, unused, drop = FALSE]
## de-mean candidate columns once for correlation
Pu_c <- scale(Pu, center = TRUE, scale = FALSE)
Pu_ss <- sqrt(colSums(Pu_c^2)); Pu_ss[Pu_ss == 0] <- 1

cat("\n__ computing cluster residuals (", length(targets), " target clusters) __\n", sep = "")
rows <- lapply(targets, function(cl) {
  y <- P[, cl]; p <- Pb[cl, ]
  pred <- as.vector(S %*% p)
  leftover <- y - pred
  rn <- sqrt(mean(leftover^2))
  ## project leftover onto unused pops (Pearson r over donor palette)
  lc <- leftover - mean(leftover); lss <- sqrt(sum(lc^2)); if (lss == 0) lss <- 1
  scores <- as.vector(crossprod(Pu_c, lc)) / (Pu_ss * lss)
  scores[unused == cl] <- NA                       # ignore self
  j <- which.max(scores)
  data.frame(pop_id = cl, res_norm = rn,
             top_source = src_pops[which.max(p)], top_p = max(p),
             miss_pop = unused[j], miss_score = scores[j],
             stringsAsFactors = FALSE)
})
cl_res <- bind_rows(rows)
write_tsv(cl_res, paste0(OUT, ".cluster_residuals.tsv"))

## ================= PER-SOURCE FLAGS =================
deep_re <- "UpperPal|Palaeolithic|Paleolithic|Mesolithic|Epipal|Hoabin|Jomon|HG_|_HG|Neolithic"
disc <- tibble(source_pop = src_pops,
               disc = sapply(src_pops, function(s) mean(abs(Pn[targets, s] - Pb[targets, s]))))
flag <- lapply(src_pops, function(s) {
  pv <- Pb[targets, s]
  loaded <- targets[pv >= 0.15]
  coupling <- suppressWarnings(cor(pv, cl_res$res_norm[match(targets, cl_res$pop_id)]))
  sub <- cl_res |> filter(pop_id %in% loaded)
  mm <- if (nrow(sub)) sub |> count(miss_pop) |> slice_max(n, n = 1, with_ties = FALSE) else tibble(miss_pop = NA, n = 0)
  data.frame(
    source_pop = s,
    is_deep = grepl(deep_re, s),
    n_loaded = length(loaded),
    mean_p_loaded = ifelse(length(loaded), mean(pv[pv >= 0.15]), 0),
    mean_res_loaded = ifelse(nrow(sub), mean(sub$res_norm), NA),
    mean_miss_score = ifelse(nrow(sub), mean(sub$miss_score), NA),
    coupling_p_res = coupling,
    modal_miss_pop = mm$miss_pop[1],
    stringsAsFactors = FALSE)
}) |> bind_rows() |> left_join(disc, by = "source_pop")

## two proxy-abuse failure modes for deep sources with >=3 loaded clusters:
##  (A) "absorber": fits well where loaded (residual bought down: coupling<=0),
##      but redundant/structured -- structured leftover OR estimator discordance.
##  (B) "poor_fit": fits *worse* where loaded (coupling>0) and leaves a large,
##      structured residual (mean_res above the typical cluster) pointing at an
##      unmodeled population -- a deep source forced onto ancestry it can't span.
global_med_res <- median(cl_res$res_norm, na.rm = TRUE)
flag <- flag |>
  mutate(
    absorber = is_deep & n_loaded >= 3 &
      (coupling_p_res <= 0.05 | is.na(coupling_p_res)) &
      (mean_miss_score >= 0.30 | disc >= 0.08),
    poor_fit = is_deep & n_loaded >= 3 &
      !is.na(coupling_p_res) & coupling_p_res > 0.05 &
      mean_miss_score >= 0.30 & mean_res_loaded > global_med_res,
    flag_type = case_when(absorber ~ "absorber", poor_fit ~ "poor_fit", TRUE ~ "")
  ) |>
  arrange(flag_type == "", desc(mean_miss_score))
write_tsv(flag, paste0(OUT, ".source_flags.tsv"))

cat(sprintf("\n(global median cluster res_norm = %.5f; poor_fit threshold)\n", global_med_res))
cat("=== PER-SOURCE PROXY-ABUSE FLAGS (deep sources, flagged first) ===\n")
print(as.data.frame(flag |> filter(is_deep) |>
  transmute(source_pop, n_loaded, mean_p = round(mean_p_loaded, 2),
            res = round(mean_res_loaded, 4), miss_score = round(mean_miss_score, 2),
            couple = round(coupling_p_res, 2), disc = round(disc, 3),
            modal_miss_pop, FLAG = flag_type)), row.names = FALSE)

cat("\n=== WORST-FIT CLUSTERS (top res_norm) with what's-missing ===\n")
print(as.data.frame(cl_res |> arrange(desc(res_norm)) |> head(20) |>
  transmute(pop_id, res = round(res_norm, 4), top_source, top_p = round(top_p, 2),
            miss_pop, miss_score = round(miss_score, 2))), row.names = FALSE)
cat("\n__ wrote", paste0(OUT, ".cluster_residuals.tsv /"), paste0(OUT, ".source_flags.tsv __\n"))
