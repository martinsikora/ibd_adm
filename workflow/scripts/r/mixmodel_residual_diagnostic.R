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
##
## It also runs a PER-STRATUM SINK TEST (see below), which catches a failure
## mode the global flags miss: a source that is not ancestry at all but an
## unfitted intercept, soaking up drift wherever the panel has no proximate
## source. Observed case: Ust'-Ishim/Ranis on the ho_20260806 h05 tier3 panel
## scored coupling_p_res = 0.042 globally (unflagged) while inside South Asia it
## correlated +0.70 with res_norm and +0.54 with cluster endogamy, taking 0.74 of
## Pulliyar and 0.60 of Palliyar -- the most inbred clusters in the panel. Its
## sibling dead-end Tianyuan/AR33K, over the same targets, scored +0.15 and +0.14
## and was tracking a real ancestry axis (r = -0.74 against the West Eurasian
## sources, i.e. the ASI/AASI cline). Both are >40 ka lineages with no
## descendants; only the stratified test tells them apart.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

## ---- paths (positional args) ----
## Without the NNLS table only the Bayesian-based outputs are written (cluster_residuals,
## source_sink_by_stratum). source_flags.tsv needs both estimators because it reports the
## NNLS/Bayesian discordance `disc`, which also feeds the absorber flag.
a <- commandArgs(trailingOnly = TRUE)
POP  <- a[1]   # pop_prof.tsv : pop_id1, pop_id2, sum   (all samples -> palette P)
SRCP <- a[2]   # src_prof.tsv : pop_id1, pop_id2, sum   (source samples -> S)
VALP <- a[3]   # val_prof.tsv : sample1, pop_id2, sum    (validation samples)
SMAP <- a[4]   # sample_map.tsv (sample_id, pop_id)
NNLS <- a[5]   # model nnls tsv, or "-" to skip everything that needs it (source_flags.tsv)
BAYE <- a[6]   # model bayesian tsv
VALI <- a[7]   # val.ids (validation sample ids)
OUT  <- a[8]   # output prefix
## optional 9th arg: quantile of the source distality distribution above which a
## source counts as distal enough for the proxy-abuse flags (default 0.5)
HAVE_NNLS <- nzchar(NNLS) && NNLS != "-"
DISTAL_Q <- if (length(a) >= 9 && nzchar(a[9])) as.numeric(a[9]) else 0.5
if (is.na(DISTAL_Q) || DISTAL_Q < 0 || DISTAL_Q > 1) {
  stop("distal quantile (arg 9) must be in [0, 1]")
}
## optional args for the per-stratum sink test (see the SINK section below)
##  10: number of strata to cut the target clusters into (data-driven, from the
##      palette geometry -- NOT a metadata region column)
##  11: |r| against res_norm or endogamy at which a source counts as a sink
##  12: minimum targets in a stratum before its correlations are trusted
##  13: minimum mean weight in the stratum -- below this a source has no material
##      influence there and a correlation is not worth a flag
##  14: minimum lead over the next-best source in the stratum (uniqueness margin)
N_STRATA   <- if (length(a) >= 10 && nzchar(a[10])) as.integer(a[10]) else 12L
SINK_MIN_R <- if (length(a) >= 11 && nzchar(a[11])) as.numeric(a[11]) else 0.4
SINK_MIN_N <- if (length(a) >= 12 && nzchar(a[12])) as.integer(a[12]) else 15L
SINK_MIN_P <- if (length(a) >= 13 && nzchar(a[13])) as.numeric(a[13]) else 0.01
SINK_MIN_GAP <- if (length(a) >= 14 && nzchar(a[14])) as.numeric(a[14]) else 0.15
if (is.na(N_STRATA) || N_STRATA < 1) stop("n_strata (arg 10) must be >= 1")
if (is.na(SINK_MIN_R) || SINK_MIN_R <= 0 || SINK_MIN_R > 1) {
  stop("sink_min_r (arg 11) must be in (0, 1]")
}
if (is.na(SINK_MIN_N) || SINK_MIN_N < 3) stop("sink_min_n (arg 12) must be >= 3")
if (is.na(SINK_MIN_P) || SINK_MIN_P < 0 || SINK_MIN_P > 1) {
  stop("sink_min_p (arg 13) must be in [0, 1]")
}
if (is.na(SINK_MIN_GAP) || SINK_MIN_GAP < 0 || SINK_MIN_GAP > 2) {
  stop("sink_min_gap (arg 14) must be in [0, 2]")
}

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
mb <- read_p(BAYE)
mn <- if (HAVE_NNLS) read_p(NNLS) else NULL

cluster_p <- function(m) m |>
  filter(!grepl("_r$", pop_id)) |>
  group_by(pop_id, source_pop) |> summarise(p = mean(p), .groups = "drop")
cpb <- cluster_p(mb)
cpn <- if (HAVE_NNLS) cluster_p(mn) else NULL

## wide cluster x source p matrices aligned to src_pops
to_wide <- function(cp) {
  w <- cp |> filter(source_pop %in% src_pops) |>
    pivot_wider(names_from = source_pop, values_from = p, values_fill = 0)
  cl <- w$pop_id; W <- as.matrix(w[, -1]); rownames(W) <- cl
  miss <- setdiff(src_pops, colnames(W)); if (length(miss)) W <- cbind(W, matrix(0, nrow(W), length(miss), dimnames = list(NULL, miss)))
  W[, src_pops, drop = FALSE]
}
Pb <- to_wide(cpb)
Pn <- if (HAVE_NNLS) to_wide(cpn) else NULL

## ================= VALIDATION: recompute res_norm for a few samples =================
val_ids <- readLines(VALI)
if (HAVE_NNLS) {
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

## ================= STRATA + ENDOGAMY (inputs to the sink test) =================
## Strata are cut from the palette geometry itself -- ward.D2 on the TVD between
## target-cluster profiles -- not from a metadata region column. On this panel
## the top-level cut recovers the continental blocks anyway, and deriving them
## from the data keeps the diagnostic working on a default panel whose pop_id is
## a bare cluster label with no geography in it.
cat("__ cutting", length(targets), "target clusters into", N_STRATA, "strata __\n")
k <- min(N_STRATA, max(1L, length(targets) - 1L))
if (k > 1L) {
  ## TVD = half the L1 distance between two normalized profiles
  d_t <- stats::dist(t(P[, targets, drop = FALSE]), method = "manhattan") / 2
  stratum <- stats::cutree(stats::hclust(d_t, method = "ward.D2"), k = k)
  stratum <- setNames(paste0("S", stratum), targets)
} else {
  stratum <- setNames(rep("S1", length(targets)), targets)
}

## Endogamy: how much more of a cluster's IBD stays inside itself than its share
## of the panel would predict. P is column-normalized, so P[cl, cl] is that
## self-share directly; dividing by n_cl / N removes the mechanical growth of the
## diagonal with cluster size (r = +0.26 between raw self-share and n here).
n_cl <- smap |> mutate(b = sub("_r$", "", pop_id)) |> count(b) |> tibble::deframe()
n_tot <- sum(n_cl)
self_share <- vapply(targets, function(cl) {
  if (cl %in% rownames(P)) P[cl, cl] else NA_real_
}, numeric(1))
exp_share <- as.numeric(n_cl[targets]) / n_tot
endogamy <- self_share / exp_share            # 1 = panel-average, >1 = inbred

cl_res <- cl_res |>
  left_join(
    tibble(pop_id = targets, stratum = unname(stratum[targets]),
           n_samples = as.numeric(n_cl[targets]), endogamy = unname(endogamy)),
    by = "pop_id"
  )
write_tsv(cl_res, paste0(OUT, ".cluster_residuals.tsv"))

## ================= PER-SOURCE FLAGS =================
## How distal each source is from the target mass: mean TVD between the source
## palette and the target-cluster palettes, weighted by target sample count.
##
## This replaces a regex over source_pop ("UpperPal|Mesolithic|Jomon|...") that
## was used to decide which sources were "deep" enough to flag. That test only
## worked on panels whose pop_id embeds a descriptive alias: on any default
## panel the pop_id is a bare cluster label (C5_2_1_0), the regex never matched,
## and every flag column silently came out FALSE. Distality is read off the same
## palette geometry the diagnostic already uses, so it behaves identically
## whatever the labelling scheme.
tgt_w <- smap |>
  mutate(b = sub("_r$", "", pop_id)) |>
  count(b) |>
  tibble::deframe()
w_t <- as.numeric(tgt_w[targets])
w_t[is.na(w_t)] <- 1
distality <- vapply(src_pops, function(s) {
  d <- 0.5 * colSums(abs(P[, targets, drop = FALSE] - S[, s]))
  stats::weighted.mean(d, w = w_t)
}, numeric(1))
distal_thr <- stats::quantile(distality, DISTAL_Q, names = FALSE)
cat(sprintf(
  "\n__ distality: median %.4f, threshold %.4f (q=%.2f) -> %d of %d sources distal __\n",
  stats::median(distality), distal_thr, DISTAL_Q,
  sum(distality >= distal_thr), length(src_pops)
))

disc <- tibble(source_pop = src_pops,
               disc = if (HAVE_NNLS) sapply(src_pops, function(s) mean(abs(Pn[targets, s] - Pb[targets, s]))) else NA_real_)
flag <- lapply(src_pops, function(s) {
  pv <- Pb[targets, s]
  loaded <- targets[pv >= 0.15]
  coupling <- suppressWarnings(cor(pv, cl_res$res_norm[match(targets, cl_res$pop_id)]))
  sub <- cl_res |> filter(pop_id %in% loaded)
  mm <- if (nrow(sub)) sub |> count(miss_pop) |> slice_max(n, n = 1, with_ties = FALSE) else tibble(miss_pop = NA, n = 0)
  data.frame(
    source_pop = s,
    distality = distality[[s]],
    is_distal = distality[[s]] >= distal_thr,
    n_loaded = length(loaded),
    mean_p_loaded = ifelse(length(loaded), mean(pv[pv >= 0.15]), 0),
    mean_res_loaded = ifelse(nrow(sub), mean(sub$res_norm), NA),
    mean_miss_score = ifelse(nrow(sub), mean(sub$miss_score), NA),
    coupling_p_res = coupling,
    modal_miss_pop = mm$miss_pop[1],
    stringsAsFactors = FALSE)
}) |> bind_rows() |> left_join(disc, by = "source_pop")

## two proxy-abuse failure modes for distal sources with >=3 loaded clusters:
##  (A) "absorber": fits well where loaded (residual bought down: coupling<=0),
##      but redundant/structured -- structured leftover OR estimator discordance.
##  (B) "poor_fit": fits *worse* where loaded (coupling>0) and leaves a large,
##      structured residual (mean_res above the typical cluster) pointing at an
##      unmodeled population -- a deep source forced onto ancestry it can't span.
## ================= PER-STRATUM SINK TEST =================
## A source can be an ancestry component or an unfitted intercept, and the two
## are indistinguishable panel-wide: a sink is quiet everywhere it has competition
## and only fires inside the one stratum where the panel runs out of proximate
## sources, so its global correlation with the residual is diluted to nothing.
## Within a stratum the signature is unambiguous -- weight that rises with
## res_norm is absorbing misfit, weight that rises with endogamy is absorbing
## drift, and neither is ancestry. A real component correlates with neither, and
## loads on well-fit clusters (see the Tianyuan/Ust'-Ishim contrast in the header).
## for each element, the largest of the OTHER elements (-Inf if there is none)
max_other <- function(v) {
  vapply(seq_along(v), function(i) {
    o <- v[-i]; o <- o[!is.na(o)]
    if (!length(o)) -Inf else max(o)
  }, numeric(1))
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok]))
}
res_by_cl   <- setNames(cl_res$res_norm, cl_res$pop_id)
endog_by_cl <- setNames(log10(cl_res$endogamy), cl_res$pop_id)
strata_ids  <- sort(unique(stratum[targets]))

cat("\n__ per-stratum sink test (", length(strata_ids), " strata; r_res >= ",
    SINK_MIN_R, ", unique with lead >= ", SINK_MIN_GAP, ", min ", SINK_MIN_N,
    " targets, mean_p >= ", SINK_MIN_P, ") __\n", sep = "")
sink_tbl <- lapply(strata_ids, function(st) {
  tin <- targets[stratum[targets] == st]
  ## name the stratum by its largest member so "S3" is identifiable in the output
  exemplar <- tin[which.max(as.numeric(n_cl[tin]))]
  rn <- res_by_cl[tin]; en <- endog_by_cl[tin]
  lapply(src_pops, function(s) {
    pv <- Pb[tin, s]
    data.frame(
      stratum = st, exemplar = exemplar, n_targets = length(tin),
      source_pop = s, n_loaded = sum(pv >= 0.15), mean_p = mean(pv),
      r_res = safe_cor(pv, rn), r_endog = safe_cor(pv, en),
      stringsAsFactors = FALSE)
  }) |> bind_rows()
}) |> bind_rows() |>
  ## UNIQUENESS. A stratum that hides a badly-fit sub-block pulls EVERY distant
  ## source onto it at once, so several sources correlate with res_norm together
  ## and none of them is specifically the intercept. Measured: in the 230-cluster
  ## European stratum EastEurope_Mesolithic (+0.52) and Georgia_UP (+0.49) tie,
  ## and in the AASI sub-block Japan_Jomon (+0.60, at mean_p 0.000), Tianyuan
  ## (+0.56) and SouthAmerica_Paleoindian (+0.50) all clear the bar. A real sink
  ## stands alone: Ust'-Ishim over South Asia is +0.71 with the runner-up at
  ## +0.39, Tagalog over Melanesia +0.72 with the runner-up at 0.00. So require
  ## the source to be the only one over the bar AND clear of the next by a
  ## margin; the shared case is reported as a stratum property instead.
  group_by(stratum) |>
  mutate(
    n_over_r = sum(!is.na(r_res) & r_res >= SINK_MIN_R),
    ## distance to the best OTHER source in this stratum; only meaningful for the
    ## top-ranked one, which is the only row that can pass n_over_r == 1 anyway
    r_res_gap = r_res - max_other(r_res),
    shared_gradient = n_over_r >= 2
  ) |>
  ungroup() |>
  mutate(
    r_res_gap = ifelse(is.finite(r_res_gap), r_res_gap, r_res),
    ## Coupling to res_norm is NECESSARY, and drift coupling only corroborates.
    ## Drift coupling alone is not evidence of a sink: a genuinely proximate
    ## source takes a larger share in the more endogamous members of its own
    ## stratum simply because they are less admixed. Measured here: JuHoan scores
    ## r_endog = +0.57 over the KhoeSan stratum while being exactly the right
    ## source for it (r_res = +0.20, fits fine). Absorbing drift is only a fault
    ## when it is also buying down misfit.
    drift_coupled = !is.na(r_endog) & r_endog >= SINK_MIN_R,
    sink = n_targets >= SINK_MIN_N & n_loaded >= 3 & mean_p >= SINK_MIN_P &
      !is.na(r_res) & r_res >= SINK_MIN_R &
      n_over_r == 1 & r_res_gap >= SINK_MIN_GAP,
    sink_via = case_when(!sink ~ "", drift_coupled ~ "misfit+drift", TRUE ~ "misfit")
  ) |>
  arrange(desc(sink), desc(r_res))
write_tsv(sink_tbl, paste0(OUT, ".source_sink_by_stratum.tsv"))

## ================= CLADE CONFINEMENT =================
## Distality alone cannot tell a BASAL ANCESTOR of the targets from an
## UNANCHORED profile: both sit close to the target mass and both look diffuse.
## The split is at which level the diffuseness lives. A lineage basal to one
## clade shares broadly *inside* that clade and little outside it; a genuinely
## unanchored profile is smeared across clades. Measuring the effective number
## of strata a source's profile spans separates them where distality inverts
## them -- NEO283 (Kotias Klde 25.7 ka, Dzudzuana-related, ancestral to later
## West Eurasian farmers and hunter-gatherers) has the LOWEST distality of any
## source here, 0.374, yet spans only 4.05 strata with 83.5% of its profile on
## West Eurasia, landing beside Satsurblia (4.14) and EastEurope_Mesolithic
## (4.11). Ust'-Ishim, more distal at 0.538, spans 8.65 strata with a top share
## of 0.213 and only 34.4% on West Eurasia. Low distality plus few strata is a
## deep source doing real work; many strata is the intercept.
cat("__ clade confinement per source __\n")
strat_of <- stratum[targets]
clade <- lapply(src_pops, function(s) {
  v <- S[, s]
  pops <- intersect(names(v), targets)
  w <- tapply(v[pops], strat_of[pops], sum)
  w <- w[!is.na(w) & w > 0]
  if (!length(w) || sum(w) <= 0) {
    return(data.frame(source_pop = s, n_eff_strata = NA_real_,
                      top_stratum = NA_character_, top_stratum_share = NA_real_,
                      stringsAsFactors = FALSE))
  }
  q <- w / sum(w)
  data.frame(source_pop = s, n_eff_strata = exp(-sum(q * log(q))),
             top_stratum = names(q)[which.max(q)], top_stratum_share = max(q),
             stringsAsFactors = FALSE)
}) |> bind_rows()

## roll up to one row per source for the flags table
sink_src <- sink_tbl |>
  group_by(source_pop) |>
  summarise(
    sink_strata = sum(sink),
    max_r_res = suppressWarnings(max(r_res, na.rm = TRUE)),
    max_r_endog = suppressWarnings(max(r_endog, na.rm = TRUE)),
    sink_stratum = ifelse(any(sink), exemplar[sink][which.max(r_res[sink])], NA_character_),
    .groups = "drop") |>
  mutate(across(c(max_r_res, max_r_endog), ~ ifelse(is.finite(.x), .x, NA_real_)))

global_med_res <- median(cl_res$res_norm, na.rm = TRUE)
flag <- flag |>
  left_join(sink_src, by = "source_pop") |>
  left_join(clade, by = "source_pop") |>
  mutate(
    sink = sink_strata > 0,
    absorber = is_distal & n_loaded >= 3 &
      (coupling_p_res <= 0.05 | is.na(coupling_p_res)) &
      (mean_miss_score >= 0.30 | disc >= 0.08),
    poor_fit = is_distal & n_loaded >= 3 &
      !is.na(coupling_p_res) & coupling_p_res > 0.05 &
      mean_miss_score >= 0.30 & mean_res_loaded > global_med_res,
    ## sink last so the existing two labels keep their meaning; a source can be
    ## both, and the boolean columns stay independent of flag_type
    flag_type = case_when(absorber ~ "absorber", poor_fit ~ "poor_fit",
                          sink ~ "sink", TRUE ~ "")
  ) |>
  arrange(flag_type == "", desc(mean_miss_score))
if (HAVE_NNLS) write_tsv(flag, paste0(OUT, ".source_flags.tsv"))

cat(sprintf("\n(global median cluster res_norm = %.5f; poor_fit threshold)\n", global_med_res))
cat("=== PER-SOURCE PROXY-ABUSE FLAGS (distal sources, flagged first) ===\n")
print(as.data.frame(flag |> filter(is_distal) |>
  transmute(source_pop, dist = round(distality, 3), strata = round(n_eff_strata, 1),
            n_loaded, mean_p = round(mean_p_loaded, 2),
            res = round(mean_res_loaded, 4), miss_score = round(mean_miss_score, 2),
            couple = round(coupling_p_res, 2), disc = round(disc, 3),
            modal_miss_pop, FLAG = flag_type)), row.names = FALSE)

cat("\n=== PER-STRATUM SINKS (source absorbing misfit/drift inside one stratum) ===\n")
sink_hits <- sink_tbl |> filter(sink)
if (nrow(sink_hits) == 0) {
  cat("  none\n")
} else {
  print(as.data.frame(sink_hits |>
    transmute(stratum, exemplar = substr(exemplar, 1, 34), n_targets, source_pop,
              n_loaded, mean_p = round(mean_p, 3), r_res = round(r_res, 2),
              lead = round(r_res_gap, 2), r_endog = round(r_endog, 2),
              via = sink_via)), row.names = FALSE)
  cat("\n  A sink is not ancestry: its weight is buying down misfit (r_res) in a\n",
      "  stratum with no proximate source, optionally also absorbing drift\n",
      "  (r_endog). Add a closer source there, or read the weight as a\n",
      "  diagnostic, not a proportion.\n", sep = "")
}

## Strata where several sources tie on res_norm: the misfit is a property of a
## sub-block, not of any one source, so no source-level flag is raised. The
## actionable read is the same -- a population is missing -- but it is the
## stratum that is telling you, and cluster_residuals$miss_pop names the target.
shared <- sink_tbl |>
  filter(shared_gradient, !is.na(r_res), r_res >= SINK_MIN_R) |>
  arrange(stratum, desc(r_res))
cat("\n=== STRATA WITH A SHARED RESIDUAL GRADIENT (sub-block missing a source) ===\n")
if (nrow(shared) == 0) {
  cat("  none\n")
} else {
  print(as.data.frame(shared |>
    transmute(stratum, exemplar = substr(exemplar, 1, 34), n_targets, source_pop,
              mean_p = round(mean_p, 3), r_res = round(r_res, 2))), row.names = FALSE)
  cat("\n  >=2 sources track the residual together, so none is singled out as a\n",
      "  sink. Check miss_pop for the worst-fit members of these strata.\n", sep = "")
}

cat("\n=== WORST-FIT CLUSTERS (top res_norm) with what's-missing ===\n")
print(as.data.frame(cl_res |> arrange(desc(res_norm)) |> head(20) |>
  transmute(pop_id, res = round(res_norm, 4), top_source, top_p = round(top_p, 2),
            miss_pop, miss_score = round(miss_score, 2))), row.names = FALSE)
cat("\n__ wrote", paste0(OUT, ".cluster_residuals.tsv /"),
    paste0(OUT, ".source_flags.tsv /"),
    paste0(OUT, ".source_sink_by_stratum.tsv __\n"))
