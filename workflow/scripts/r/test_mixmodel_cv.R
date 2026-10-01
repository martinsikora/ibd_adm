#!/usr/bin/env Rscript
## Tests for mixmodel_cv.R. Run from the repo root:
##   Rscript workflow/scripts/r/test_mixmodel_cv.R
## Pulls the pure helper functions it depends on straight out of mixmodel_ibd.R
## (evaluating only those assignments), so it tests the code that is run.
suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(future); library(furrr); library(lsei)
})
plan(sequential)

pull_fns <- function(file, names) {
  for (e in parse(file)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% names) {
      eval(e, globalenv())
    }
  }
}
pull_fns("workflow/scripts/r/mixmodel_ibd.R",
  c("get_sum_matrix", "normalize_matrix_cols", "res_norm_ex_self", "align_matrix_cols"))
source("workflow/scripts/r/mixmodel_cv.R")

ok <- function(cond, msg) {
  cat(sprintf("%s  %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) quit(status = 1)
}
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

## ---- folds
nm <- read.delim("config/n_markers.tsv"); ch <- as.character(nm$chrom); n <- nm$n
eo <- cv_make_folds("evenodd", ch, n)
ok(length(eo) == 2 && setequal(c(eo[[1]]$test, eo[[2]]$test), ch) &&
   length(intersect(eo[[1]]$test, eo[[2]]$test)) == 0, "evenodd: two disjoint folds covering all chromosomes")
ok(all(as.integer(eo[[1]]$train) %% 2 == 0) && all(as.integer(eo[[1]]$test) %% 2 == 1), "evenodd: train_even fits even, scores odd")
odd_share <- sum(n[as.integer(ch) %% 2 == 1]) / sum(n)
cat(sprintf("     odd chromosomes carry %.1f%% of markers\n", 100 * odd_share))
ok(abs(odd_share - 0.5) < 0.05, "evenodd: halves within 5% of balanced on markers")
ok(length(cv_make_folds("loco", ch, n)) == length(ch), "loco: one fold per chromosome")
k5 <- cv_make_folds("k5", ch, n)
ld <- sapply(k5, function(f) sum(n[ch %in% f$test]))
cat(sprintf("     k5 block marker shares: %s\n", paste(sprintf("%.1f%%", 100 * ld / sum(n)), collapse = " ")))
ok(length(k5) == 5 && setequal(unlist(lapply(k5, `[[`, "test")), ch) && max(ld) / min(ld) < 1.1, "k5: partition, balanced within 10%")
ok(identical(cv_make_folds("k5", ch, n), k5), "k5: deterministic")
cu <- cv_make_folds("test:1,3-5", ch, n)
ok(identical(cu[[1]]$test, c("1", "3", "4", "5")) && length(cu[[1]]$train) == length(ch) - 4, "test:<chroms> range expansion")
ok(length(cv_make_folds("none", ch, n)) == 0, "none: no folds")
ok(fails(cv_make_folds("bogus", ch, n)), "unknown spec aborts")
ok(fails(cv_make_folds("test:1,99", ch, n)), "unknown chromosome in test: aborts")
ok(fails(cv_make_folds("k1", ch, n)) && fails(cv_make_folds("k99", ch, n)), "k out of range aborts")
ok(fails(cv_make_folds("evenodd", c("1", "2", "X"), c(1, 1, 1))), "evenodd with non-integer chromosome aborts")
ok(fails(cv_make_folds("test:1-22", ch, n)), "fold with no train chromosome aborts")
ok(fails(cv_check_folds_present(eo, c("1", "2", "3"))), "fold naming chromosomes with no input aborts")
ok(isTRUE(cv_check_folds_present(eo, ch)), "all chromosomes present passes")

## ---- get_sum_matrix on an empty chromosome filter (all-zero, not an error)
d0 <- tibble(sample1 = character(), pop_id1 = character(), pop_id2 = character(), chrom = character(), ibd = double())
pops <- c("A", "B", "C")
z <- try(align_matrix_cols(get_sum_matrix(d0, pops, sample1, pop_id2), c("t1", "t2")), silent = TRUE)
ok(!inherits(z, "try-error") && all(z == 0) && identical(dim(z), c(3L, 2L)), "get_sum_matrix + align on empty input gives a zero matrix")

## ---- synthetic: correct model ~ noise floor, over-parameterised model overfits
set.seed(42)
D <- 300; C <- 22; N_chr_t <- 120; N_chr_s <- 1500   # target IBD per chromosome, source IBD per chromosome
rdir <- function(k, a) { x <- rgamma(k, a); x / sum(x) }
K_true <- 5; K_extra <- 45
mu <- replicate(K_true + K_extra, rdir(D, 0.15))       # donor x source expected palettes
p0 <- c(0.5, 0.3, 0.2, 0, 0)
pi_t <- as.vector(mu[, 1:K_true] %*% p0)
n_tg <- 300   # the overfit effect is ~1% of the error, so it needs paired tests over many targets
tn <- paste0("t", seq_len(n_tg)); sn <- paste0("s", seq_len(K_true + K_extra))
by_t <- lapply(seq_len(C), function(i) sapply(tn, function(x) as.vector(rmultinom(1, N_chr_t, pi_t))) |> `dimnames<-`(list(paste0("d", 1:D), tn)))
by_s <- lapply(seq_len(C), function(i) sapply(seq_len(ncol(mu)), function(k) as.vector(rmultinom(1, N_chr_s, mu[, k]))) |> `dimnames<-`(list(paste0("d", 1:D), sn)))
sum_over <- function(m, ix) Reduce(`+`, m[ix])
fold <- cv_make_folds("evenodd", as.character(1:C), rep(1, C))[[1]]
tr <- as.integer(fold$train); te <- as.integer(fold$test)
self_row <- setNames(rep(1L, n_tg), tn); self_src <- setNames(rep(FALSE, n_tg), tn)
run <- function(cols) {
  cv_run_fold(fold, cv_fit_nnls,
    Ttr = normalize_matrix_cols(sum_over(by_t, tr)),
    Str = normalize_matrix_cols(sum_over(by_s, tr))[, cols, drop = FALSE],
    Tte_raw = sum_over(by_t, te),
    Ste_raw = sum_over(by_s, te)[, cols, drop = FALSE],
    self_row = self_row, self_src = self_src, n_chunks = 4)
}
good <- run(1:K_true)                 # the true source set
over <- run(seq_len(K_true + K_extra)) # true + 45 decoys
m <- function(x) median(x, na.rm = TRUE)
cat(sprintf("     correct: train %.5f test %.5f ll %.4f active %.1f\n", m(good$res_train_ex_self), m(good$res_test_ex_self), m(good$ll_test), m(good$n_active)))
cat(sprintf("     overfit: train %.5f test %.5f ll %.4f active %.1f\n", m(over$res_train_ex_self), m(over$res_test_ex_self), m(over$ll_test), m(over$n_active)))
keep <- seq_len(D) != 1
pt <- pi_t[keep] / sum(pi_t[keep])
floor_rmse <- sqrt(mean(pt * (1 - pt)) / (length(te) * N_chr_t * (1 - pi_t[1])))
cat(sprintf("     theoretical test-palette noise floor (RMSE): %.5f\n", floor_rmse))
ok(m(good$res_test_ex_self) > 0.9 * floor_rmse && m(good$res_test_ex_self) < 1.6 * floor_rmse, "correct model: held-out error within [0.9, 1.6] x the sampling noise floor")
paired <- function(x, label, expect) {
  tt <- t.test(x, alternative = expect)
  ok(tt$p.value < 1e-3, sprintf("paired %s: mean diff %+.2e, %.0f%% of targets, p=%.1e", label, mean(x), 100 * mean(x > 0), tt$p.value))
}
## overfit minus correct: train error should be LOWER, held-out error HIGHER, held-out ll LOWER
paired(over$res_train_ex_self - good$res_train_ex_self, "train RMSE (overfit fits better in-sample)", "less")
paired(over$res_test_ex_self - good$res_test_ex_self, "held-out RMSE (overfit generalises worse)", "greater")
paired(good$ll_test - over$ll_test, "held-out log-lik (correct beats overfit)", "greater")
ok(all(good$n_src_zero_test == 0) && all(is.finite(good$cm_test)), "no spurious zero-source flags in the healthy case")
ok(all(c("fold", "n_chr_train", "n_chr_test", "self_is_source", "p_lost_test") %in% colnames(good)) && nrow(good) == n_tg, "output shape")

## zero-source diagnostic: a source with no test-chromosome IBD must be flagged and its mass reported
Ste <- sum_over(by_s, te)[, 1:K_true, drop = FALSE]; Ste[, 1] <- 0
zz <- cv_run_fold(fold, cv_fit_nnls, normalize_matrix_cols(sum_over(by_t, tr)),
  normalize_matrix_cols(sum_over(by_s, tr))[, 1:K_true], sum_over(by_t, te), Ste, self_row, self_src, 2)
ok(all(zz$n_src_zero_test == 1) && m(zz$p_lost_test) > 0.3, sprintf("zero test column flagged (lost mass median %.2f)", m(zz$p_lost_test)))
cat("ALL PASSED\n")
