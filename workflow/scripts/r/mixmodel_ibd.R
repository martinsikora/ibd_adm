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
  library(purrr)
  library(tidyr)
  library(future)
  library(furrr)
  library(lsei)
})

options(future.globals.maxSize = 5e9)

## --------------------------------------------------
## functions

get_sum_matrix <- function(d, groups, g_x, g_r) {
  ## set up vars
  g_x1 <- enquo(g_x)
  g_r1 <- enquo(g_r)

  ## get grouped summary and reshape (raw sums)
  r1 <- d |>
    group_by(!!g_x1, !!g_r1) |>
    summarise(
      value = sum(ibd),
      .groups = "drop_last"
    ) |>
    ungroup() |>
    select(!!g_x1, !!g_r1, value) |>
    pivot_wider(
      values_from = "value",
      names_from = !!g_x1,
      values_fill = 0
    )

  m1 <- r1 |>
    select(-!!g_r1) |>
    as.matrix()

  rownames(m1) <- r1 |>
    select(!!g_r1) |>
    pull(!!g_r1)

  ## set up result matrix and return results
  m2 <- matrix(0,
    ncol = ncol(m1),
    nrow = length(groups)
  )
  rownames(m2) <- groups
  colnames(m2) <- colnames(m1)
  m2[rownames(m1), ] <- m1
  m2
}

normalize_matrix_cols <- function(m) {
  m <- as.matrix(m)
  cs <- colSums(m)
  nz <- cs > 0
  out <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  rownames(out) <- rownames(m)
  colnames(out) <- colnames(m)
  if (any(nz)) {
    out[, nz] <- sweep(m[, nz, drop = FALSE], 2, cs[nz], "/")
  }
  out
}

## Fit quality on the EXTERNAL donor palette.
##
## all_pops includes the target populations themselves, so a target's palette carries
## a column of IBD with its own cluster. When that cluster is not also a source, no
## mixture of sources can predict it, and the miss enters res_norm as a fixed penalty
## scaling with cohort size and endogamy rather than with ancestry fit. In a large,
## endogamous cohort most of the squared error can sit in that one column.
##
## Both vectors are renormalized over the retained donors, so this compares palette
## SHAPE outside the own cluster. Dropping the row without renormalizing leaves y
## deflated by (1 - self_share) while pred is not; the two variants agree to
## spearman 0.995 panel-wide and diverge only where self_share is large.
##
## Always an RMSE, for both estimators. Note that the existing res_norm is NOT
## comparable between them: bayesian res_norm is an RMSE (infer_sourcefind), while
## the NNLS one is pnnls()'s rnorm, a plain L2 norm, larger by sqrt(length(y)).
res_norm_ex_self <- function(y, pred, self_row) {
  if (length(self_row) != 1L || is.na(self_row)) {
    return(NA_real_)
  }
  keep <- seq_along(y) != self_row
  sy <- sum(y[keep])
  sp <- sum(pred[keep])
  if (sy <= 0 || sp <= 0) {
    return(NA_real_)
  }
  sqrt(mean((pred[keep] / sp - y[keep] / sy)^2))
}

align_matrix_cols <- function(m, cols) {
  m <- as.matrix(m)
  out <- matrix(0,
    nrow = nrow(m),
    ncol = length(cols)
  )
  rownames(out) <- rownames(m)
  colnames(out) <- cols
  present <- intersect(colnames(m), cols)
  if (length(present) > 0) {
    out[, present] <- m[, present, drop = FALSE]
  }
  out
}

get_palette_matrix <- function(d, groups, g_x, g_r) {
  get_sum_matrix(d, groups, {{ g_x }}, {{ g_r }}) |>
    normalize_matrix_cols()
}

rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

log_ddirichlet <- function(x, alpha) {
  sum((alpha - 1) * log(pmax(x, 1e-16))) +
    lgamma(sum(alpha)) -
    sum(lgamma(alpha))
}

ess_1d <- function(x, max_lag = 200L) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3) {
    return(NA_real_)
  }
  x <- x - mean(x)
  v <- sum(x^2) / n
  if (v == 0) {
    return(n)
  }

  max_lag <- min(max_lag, n - 1)
  acf_vals <- sapply(seq_len(max_lag), function(lag_i) {
    sum(x[seq_len(n - lag_i)] * x[seq.int(1 + lag_i, n)]) / ((n - lag_i) * v)
  })

  tau <- 1
  k <- 1
  while (k <= length(acf_vals)) {
    pair_sum <- acf_vals[k] + ifelse(k + 1 <= length(acf_vals), acf_vals[k + 1], 0)
    if (is.na(pair_sum) || pair_sum < 0) {
      break
    }
    tau <- tau + 2 * pair_sum
    k <- k + 2
  }

  ess <- n / tau
  max(1, min(n, ess))
}

rhat_1d <- function(chain_draws) {
  chain_draws <- as.matrix(chain_draws)
  n <- nrow(chain_draws)
  m <- ncol(chain_draws)
  if (m < 2 || n < 2) {
    return(NA_real_)
  }

  chain_means <- colMeans(chain_draws)
  chain_vars <- apply(chain_draws, 2, var)
  w <- mean(chain_vars)
  if (!is.finite(w) || w <= 0) {
    ## 2026-08-21: this branch used to return 1 unconditionally, which reports
    ## perfect convergence for a component that never moved. That is correct when all
    ## chains sit at the same constant, but it also hid the opposite case: every
    ## chain frozen at a DIFFERENT value, which is maximal non-convergence.
    ## Zero within-chain variance with non-zero between-chain variance is
    ## R-hat = Inf, so say so and let rhat_max carry it.
    b0 <- var(chain_means)
    if (!is.finite(b0) || b0 <= 0) {
      return(1)
    }
    return(Inf)
  }
  b <- n * var(chain_means)
  var_hat <- ((n - 1) / n) * w + (1 / n) * b
  sqrt(var_hat / w)
}

safe_pnnls <- function(source_mat, y, raw = FALSE) {
  k <- ncol(source_mat)
  fit <- tryCatch(
    if (raw) pnnls(source_mat, y) else pnnls(source_mat, y, sum = 1),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$x) || length(fit$x) != k) {
    return(rep(1 / k, k))
  }
  p <- pmax(as.numeric(fit$x), 0)
  s <- sum(p)
  if (!is.finite(s) || s <= 0) {
    return(rep(1 / k, k))
  }
  p / s
}

select_active_sources_slots <- function(
    y,
    source_mat,
    max_active_sources,
    expected_active_sources,
    num_slots = 200,
    n_iter = 1500,
    burnin = 500,
    thin = 10,
    jump_prob = 0.1,
    n_copies = 20000,
    raw_scale = FALSE
) {
  k <- ncol(source_mat)
  if (k <= 1 || max_active_sources >= k) {
    return(seq_len(k))
  }

  max_active_sources <- max(2L, min(k, as.integer(max_active_sources)))
  expected_active_sources <- max(1, min(max_active_sources, expected_active_sources))
  num_slots <- max(max_active_sources, as.integer(num_slots))
  n_iter <- max(10L, as.integer(n_iter))
  burnin <- max(0L, min(as.integer(burnin), n_iter - 2L))
  thin <- max(1L, as.integer(thin))

  ll_from_slots <- function(assign) {
    p <- tabulate(assign, nbins = k) / length(assign)
    q <- model_pred(source_mat, p, raw_scale)
    n_active <- sum(p > 0)
    ll <- n_copies * sum(y * log(pmax(q, 1e-16)))
    lp <- dpois(n_active, lambda = expected_active_sources, log = TRUE)
    ll + lp
  }

  active_init <- sample.int(k, size = max_active_sources, replace = FALSE)
  slots <- sample(active_init, size = num_slots, replace = TRUE)
  ll_prev <- ll_from_slots(slots)
  slots_prev <- slots
  kept <- list()
  keep_i <- 0L

  for (it in seq_len(n_iter)) {
    uniq <- sort(unique(slots))
    counts <- as.numeric(table(factor(slots, levels = uniq)))
    donor_out <- uniq[sample.int(length(uniq), size = 1, prob = counts / sum(counts))]
    donor_slots <- which(slots == donor_out)
    n_replace <- sample.int(length(donor_slots), size = 1)
    if (runif(1) < jump_prob) {
      n_replace <- length(donor_slots)
    }
    replace_slots <- donor_slots[sample.int(length(donor_slots), size = n_replace, replace = FALSE)]

    keep_slots <- setdiff(seq_len(num_slots), replace_slots)
    active_keep <- unique(slots[keep_slots])
    candidate_pool <- setdiff(active_keep, donor_out)
    if (length(candidate_pool) < max_active_sources) {
      add_pool <- setdiff(seq_len(k), c(candidate_pool, donor_out))
      n_add <- min(length(add_pool), max_active_sources - length(candidate_pool))
      if (n_add > 0) {
        candidate_pool <- c(candidate_pool, add_pool[sample.int(length(add_pool), size = n_add, replace = FALSE)])
      }
    }
    if (runif(1) < jump_prob) {
      candidate_pool <- setdiff(seq_len(k), c(active_keep, donor_out))
    }
    if (length(candidate_pool) == 0) {
      candidate_pool <- setdiff(seq_len(k), donor_out)
    }

    slots_prop <- slots
    donor_in <- candidate_pool[sample.int(length(candidate_pool), size = 1)]
    slots_prop[replace_slots] <- donor_in
    ll_prop <- ll_from_slots(slots_prop)
    if (log(runif(1)) < (ll_prop - ll_prev)) {
      slots <- slots_prop
      ll_prev <- ll_prop
      slots_prev <- slots_prop
    } else {
      slots <- slots_prev
    }

    if (it > burnin && ((it - burnin) %% thin == 0)) {
      keep_i <- keep_i + 1L
      kept[[keep_i]] <- tabulate(slots, nbins = k) / num_slots
    }
  }

  if (length(kept) == 0) {
    p_last <- tabulate(slots, nbins = k) / num_slots
    kept <- list(p_last)
  }
  p_mean <- colMeans(do.call(rbind, kept))
  ord <- order(p_mean, decreasing = TRUE)
  selected <- ord[seq_len(max_active_sources)]
  selected[order(selected)]
}

infer_sourcefind <- function(
    y,
    source_mat,
    n_iter = 4000,
    burnin = 1000,
    thin = 10,
    alpha_prior = 0.5,
    proposal_scale = 200,
    n_copies = 20000,
    raw_scale = FALSE,
    n_chains = 4,
    adapt_burnin_frac = 0.5,
    adapt_interval = 200,
    adapt_target_accept = 0.01,
    adapt_rate = 1.5,
    proposal_scale_min = 0.1,
    proposal_scale_max = 1e6,
    local_move_prob = 0.8,
    mean_active_sources = 0,
    active_eps = 1e-4,
    hybrid_active_search = TRUE,
    max_active_sources = 0,
    active_search_slots = 200,
    active_search_iter = 1500,
    active_search_burnin = 500,
    active_search_thin = 10,
    active_search_jump_prob = 0.1,
    se_n_copies = 0,
    selected_idx_override = NULL
) {
  ## two-stage fit: the weights come from the fit with n_copies observations, the SE from a second fit on the same
  ## sources with se_n_copies observations (a posterior whose width is calibrated; see docs/DIAGNOSTICS.md)
  if (se_n_copies > 0 && se_n_copies != n_copies) {
    call_args <- as.list(environment())
    fit_p <- do.call(infer_sourcefind, modifyList(call_args, list(se_n_copies = 0)))
    fit_se <- do.call(infer_sourcefind, modifyList(call_args, list(se_n_copies = 0, n_copies = se_n_copies,
                                                                   selected_idx_override = fit_p$selected_idx)))
    fit_p$se <- fit_se$se
    return(fit_p)
  }
  k <- ncol(source_mat)
  max_active_sources <- as.integer(max_active_sources)
  hybrid_active_search <- as.logical(hybrid_active_search)

  if (max_active_sources <= 0) {
    if (mean_active_sources > 0) {
      max_active_sources <- max(2L, min(k, as.integer(ceiling(mean_active_sources * 1.5))))
    } else {
      max_active_sources <- k
    }
  } else {
    max_active_sources <- max(2L, min(k, max_active_sources))
  }

  selected_idx <- seq_len(k)
  if (!is.null(selected_idx_override)) {
    selected_idx <- selected_idx_override
  } else if (hybrid_active_search && max_active_sources < k) {
    expected_active <- ifelse(mean_active_sources > 0, mean_active_sources, max_active_sources / 2)
    selected_idx <- select_active_sources_slots(
      y = y,
      source_mat = source_mat,
      max_active_sources = max_active_sources,
      expected_active_sources = expected_active,
      num_slots = active_search_slots,
      n_iter = active_search_iter,
      burnin = active_search_burnin,
      thin = active_search_thin,
      jump_prob = active_search_jump_prob,
      n_copies = n_copies,
      raw_scale = raw_scale
    )
  }

  source_mat_fit <- source_mat[, selected_idx, drop = FALSE]
  k_fit <- ncol(source_mat_fit)
  alpha_vec <- rep(alpha_prior, k_fit)
  use_active_prior <- is.finite(mean_active_sources) && mean_active_sources > 0

  log_active_prior <- function(p) {
    if (!use_active_prior) {
      return(0)
    }
    n_active <- sum(p > active_eps)
    dpois(n_active, lambda = mean_active_sources, log = TRUE)
  }

  log_post <- function(p) {
    q <- model_pred(source_mat_fit, p, raw_scale)
    ll <- n_copies * sum(y * log(pmax(q, 1e-16)))
    lp <- log_ddirichlet(p, alpha_vec)
    la <- log_active_prior(p)
    ll + lp + la
  }

  p_curr <- rep(1 / k_fit, k_fit)
  lp_curr <- log_post(p_curr)
  n_keep <- floor((n_iter - burnin) / thin)
  if (n_keep < 2) {
    stop("Need at least 2 kept posterior samples. Increase --mcmc_iter or reduce --burnin/--thin.")
  }
  n_chains <- as.integer(n_chains)
  if (n_chains < 1) {
    stop("n_chains must be >= 1")
  }
  adapt_phase <- as.integer(max(0, min(burnin, floor(n_iter * adapt_burnin_frac))))

  init_p <- safe_pnnls(source_mat_fit, y, raw = raw_scale)
  chain_samples <- vector("list", n_chains)
  chain_accept <- numeric(n_chains)
  chain_ess <- matrix(NA_real_, nrow = n_chains, ncol = k_fit)
  chain_prop_final <- numeric(n_chains)

  for (chain_idx in seq_len(n_chains)) {
    p_curr <- if (chain_idx == 1) {
      init_p
    } else {
      rdirichlet(init_p * 200 + 1e-3)
    }
    lp_curr <- log_post(p_curr)

    samples <- matrix(NA_real_, nrow = n_keep, ncol = k_fit)
    keep_i <- 0
    accepted <- 0L
    accepted_block <- 0L
    block_len <- 0L
    proposal_scale_chain <- proposal_scale

    for (i in seq_len(n_iter)) {
      do_local_move <- (k_fit > 1) && (runif(1) < local_move_prob)
      if (do_local_move) {
        pair <- sample.int(k_fit, size = 2, replace = FALSE)
        i_from <- pair[1]
        i_to <- pair[2]
        lo <- -p_curr[i_to]
        hi <- p_curr[i_from]
        delta <- runif(1, lo, hi)
        p_prop <- p_curr
        p_prop[i_from] <- p_prop[i_from] - delta
        p_prop[i_to] <- p_prop[i_to] + delta
        lp_prop <- log_post(p_prop)
        log_a <- lp_prop - lp_curr
      } else {
        prop_alpha <- p_curr * proposal_scale_chain + 1e-3
        p_prop <- rdirichlet(prop_alpha)
        lp_prop <- log_post(p_prop)

        log_q_forward <- log_ddirichlet(p_prop, prop_alpha)
        log_q_reverse <- log_ddirichlet(p_curr, p_prop * proposal_scale_chain + 1e-3)
        log_a <- lp_prop - lp_curr + log_q_reverse - log_q_forward
      }

      accepted_step <- FALSE
      if (log(runif(1)) < log_a) {
        p_curr <- p_prop
        lp_curr <- lp_prop
        accepted <- accepted + 1L
        accepted_step <- TRUE
      }

      if (i <= adapt_phase) {
        block_len <- block_len + 1L
        if (accepted_step) {
          accepted_block <- accepted_block + 1L
        }
        if (block_len >= adapt_interval) {
          block_accept <- accepted_block / block_len
          step_scale <- exp((adapt_target_accept - block_accept) * adapt_rate)
          proposal_scale_chain <- proposal_scale_chain * step_scale
          proposal_scale_chain <- min(proposal_scale_max, max(proposal_scale_min, proposal_scale_chain))
          block_len <- 0L
          accepted_block <- 0L
        }
      }

      if (i > burnin && ((i - burnin) %% thin == 0)) {
        keep_i <- keep_i + 1
        samples[keep_i, ] <- p_curr
      }
    }

    chain_samples[[chain_idx]] <- samples
    chain_accept[chain_idx] <- accepted / n_iter
    chain_ess[chain_idx, ] <- apply(samples, 2, ess_1d)
    chain_prop_final[chain_idx] <- proposal_scale_chain
  }

  pooled_samples <- do.call(rbind, chain_samples)
  p_mean_fit <- colMeans(pooled_samples)
  p_sd_fit <- apply(pooled_samples, 2, sd)
  ess_vec <- colSums(chain_ess, na.rm = TRUE)
  active_counts <- apply(pooled_samples, 1, function(p) sum(p > active_eps))
  rhat_vec <- rep(NA_real_, k_fit)
  if (n_chains > 1) {
    for (j in seq_len(k_fit)) {
      rhat_vec[j] <- rhat_1d(sapply(chain_samples, function(s) s[, j]))
    }
  }
  p_mean <- rep(0, k)
  p_sd <- rep(0, k)
  p_mean[selected_idx] <- p_mean_fit
  p_sd[selected_idx] <- p_sd_fit
  p_med <- rep(0, k)
  p_med[selected_idx] <- apply(pooled_samples, 2, median)
  p_med <- p_med / sum(p_med)                     # posterior median per source, renormalised to sum to 1
  pred <- model_pred(source_mat, p_mean, raw_scale)
  res_norm <- sqrt(mean((pred - y)^2))

  list(
    p = p_mean,
    p_median = p_med,
    se = p_sd,
    selected_idx = selected_idx,
    res_norm = res_norm,
    accept_rate = mean(chain_accept),
    accept_rate_min = min(chain_accept),
    accept_rate_max = max(chain_accept),
    ess_min = min(ess_vec, na.rm = TRUE),
    ess_median = median(ess_vec, na.rm = TRUE),
    rhat_max = ifelse(all(is.na(rhat_vec)), NA_real_, max(rhat_vec, na.rm = TRUE)),
    rhat_median = ifelse(all(is.na(rhat_vec)), NA_real_, median(rhat_vec, na.rm = TRUE)),
    active_sources_median = median(active_counts, na.rm = TRUE),
    selected_sources_n = length(selected_idx),
    n_keep = n_keep * n_chains,
    n_chains = n_chains,
    proposal_scale_final = median(chain_prop_final)
  )
}

se_jk <- function(theta_hat, theta_j, mi) {
  ## jackknife SE given resampled estimates and block sizes (Busing 1999)
  g <- length(mi)
  if (g < 2) {
    return(NA_real_)
  }
  n <- sum(mi)
  hi <- n / mi
  t2 <- sum((1 - mi / n) * theta_j)
  t1 <- (hi * theta_hat - (hi - 1) * theta_j - g * theta_hat + t2)^2 / (hi - 1)
  sqrt(sum(t1) / g)
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
  nargs = "+",
  help = "Files with ibd sharing data for each chromosome"
)

parser$add_argument("-s", "--sample_file",
  action = "store",
  dest = "sample_file",
  help = "File with sample to group mapping for model"
)

parser$add_argument("-g", "--group_file",
  action = "store",
  dest = "group_file",
  help = "File with sample to group mapping"
)

parser$add_argument("-i", "--individuals",
  action = "store",
  dest = "individuals_file",
  help = "File with sample labels"
)

parser$add_argument("-l", "--length_file",
  action = "store",
  dest = "length_file",
  help = "File with number of markers for each chromosome"
)

parser$add_argument("-o", "--out",
  action = "store",
  dest = "out_file",
  help = "Output filename"
)

parser$add_argument("-t", "--threads",
  action = "store",
  dest = "threads",
  type = "integer",
  default = 1L,
  help = "Number of threads [default %(default)s]"
)

parser$add_argument("--mcmc_iter",
  action = "store",
  dest = "mcmc_iter",
  type = "integer",
  default = 4000L,
  help = "MCMC iterations [default %(default)s]"
)

parser$add_argument("--burnin",
  action = "store",
  dest = "burnin",
  type = "integer",
  default = 1000L,
  help = "MCMC burnin [default %(default)s]"
)

parser$add_argument("--thin",
  action = "store",
  dest = "thin",
  type = "integer",
  default = 10L,
  help = "MCMC thinning [default %(default)s]"
)

parser$add_argument("--proposal_scale",
  action = "store",
  dest = "proposal_scale",
  type = "double",
  default = 20,
  help = "Dirichlet proposal concentration scale for Bayesian mode [default %(default)s]"
)

parser$add_argument("--mcmc_chains",
  action = "store",
  dest = "mcmc_chains",
  type = "integer",
  default = 4L,
  help = "Number of Bayesian MCMC chains [default %(default)s]"
)

parser$add_argument("--adapt_burnin_frac",
  action = "store",
  dest = "adapt_burnin_frac",
  type = "double",
  default = 0.5,
  help = "Fraction of iterations used for proposal adaptation [default %(default)s]"
)

parser$add_argument("--adapt_interval",
  action = "store",
  dest = "adapt_interval",
  type = "integer",
  default = 200L,
  help = "Iterations per adaptation update [default %(default)s]"
)

parser$add_argument("--adapt_target_accept",
  action = "store",
  dest = "adapt_target_accept",
  type = "double",
  default = 0.01,
  help = "Target acceptance rate for adaptive proposal updates [default %(default)s]"
)

parser$add_argument("--local_move_prob",
  action = "store",
  dest = "local_move_prob",
  type = "double",
  default = 0.8,
  help = "Probability of local pairwise mass-transfer proposal [default %(default)s]"
)

parser$add_argument("--mean_active_sources",
  action = "store",
  dest = "mean_active_sources",
  type = "double",
  default = 0,
  help = "Poisson prior mean on number of active sources; <=0 disables prior [default %(default)s]"
)

parser$add_argument("--active_eps",
  action = "store",
  dest = "active_eps",
  type = "double",
  default = 1e-4,
  help = "Threshold for counting active source contributions [default %(default)s]"
)

parser$add_argument("--hybrid_active_search",
  action = "store",
  dest = "hybrid_active_search",
  type = "integer",
  default = 1L,
  help = "Enable slot-based active-source search before continuous refinement (0/1) [default %(default)s]"
)

parser$add_argument("--palette_scale",
  action = "store",
  dest = "palette_scale",
  type = "character",
  default = "normalized",
  help = "normalized: every palette sums to 1, so a source is weighted by its ancestry share times its total IBD per individual. raw: sources are mean per-individual palettes in cM, the target palette is fitted up to a free scale and the weights are normalised afterwards; it removes that dependence on total IBD, but a source whose palette is far smaller than the others (a single low-sharing genome) can take any weight, so do not use it with such sources [default %(default)s]"
)

parser$add_argument("--se_genome_length_cm",
  action = "store",
  dest = "se_genome_length_cm",
  type = "double",
  default = 0,
  help = "Bayesian: if > 0, fit a second time on the same sources with this many observations in the likelihood and report its posterior SD as se (the weights stay those of the fixed 20000 fit). Set it to the length of the genome covered by the IBD data in cM (about 3500 for human autosomes), the number of trials in the SOURCEFIND likelihood. 0 reports the SD of the single fit [default %(default)s]"
)

parser$add_argument("--max_active_sources",
  action = "store",
  dest = "max_active_sources",
  type = "integer",
  default = 0L,
  help = "Maximum active sources in hybrid search; <=0 auto from mean_active_sources [default %(default)s]"
)

parser$add_argument("--active_search_slots",
  action = "store",
  dest = "active_search_slots",
  type = "integer",
  default = 200L,
  help = "Number of slots in hybrid active-source search [default %(default)s]"
)

parser$add_argument("--active_search_iter",
  action = "store",
  dest = "active_search_iter",
  type = "integer",
  default = 1500L,
  help = "Iterations for hybrid active-source search [default %(default)s]"
)

parser$add_argument("--active_search_burnin",
  action = "store",
  dest = "active_search_burnin",
  type = "integer",
  default = 500L,
  help = "Burnin for hybrid active-source search [default %(default)s]"
)

parser$add_argument("--active_search_thin",
  action = "store",
  dest = "active_search_thin",
  type = "integer",
  default = 10L,
  help = "Thinning for hybrid active-source search [default %(default)s]"
)

parser$add_argument("--active_search_jump_prob",
  action = "store",
  dest = "active_search_jump_prob",
  type = "double",
  default = 0.1,
  help = "Global-jump probability in hybrid active-source search [default %(default)s]"
)

parser$add_argument("--method",
  action = "store",
  dest = "method",
  default = "nnls",
  help = "Inference method: nnls or bayesian [default %(default)s]"
)

parser$add_argument("--seed",
  action = "store",
  dest = "seed",
  type = "integer",
  default = -1L,
  help = paste(
    "RNG seed. Negative means unseeded (the historical behaviour).",
    "Without a seed, furrr_options(seed = TRUE) derives its per-worker",
    "L'Ecuyer streams from ambient RNG state, so no run is reproducible --",
    "two runs of the same config agreed on only 9 of 14 selected sources",
    "in the hybrid slot search. [default %(default)s]"
  )
)

parser$add_argument("--read_prefilter_frac",
  action = "store",
  type = "double",
  default = 0.5,
  help = paste(
    "Row-prefilter the IBD tables to the samples this panel actually uses",
    "(targets + sources) when they are at most this FRACTION of the panel.",
    "0 disables. Only rows whose sample1 is a target or a source are ever",
    "read downstream, so the filter is lossless. [default %(default)s]"
  )
)

parser$add_argument("--cv",
  action = "store",
  dest = "cv",
  default = "none",
  help = paste(
    "Chromosome hold-out evaluation (see mixmodel_cv.R): none, evenodd, loco,",
    "k<K> (K marker-balanced blocks) or test:<chroms> e.g. test:1,3-5.",
    "Each fold fits on the train chromosomes and scores on the held-out ones;",
    "results go to a separate per-target table. [default %(default)s]"
  )
)

parser$add_argument("--cv_only",
  action = "store",
  dest = "cv_only",
  type = "integer",
  default = 0L,
  help = "Run only the --cv evaluation and skip the genome-wide fit and main table (0/1) [default %(default)s]"
)

parser$add_argument("--cv_out",
  action = "store",
  dest = "cv_out",
  help = "CV table path [default: --out with .tsv replaced by .cv.tsv]"
)

args <- parser$parse_args()
## number of independent observations in the Bayesian likelihood
if (args$se_genome_length_cm < 0) stop("--se_genome_length_cm must be >= 0")

if (!(args$method %in% c("nnls", "bayesian"))) {
  stop("--method must be one of: nnls, bayesian")
}
if (args$threads < 1) {
  stop("--threads must be >= 1")
}
if (args$method == "bayesian" && (args$mcmc_iter <= args$burnin || args$thin < 1 || args$proposal_scale <= 0)) {
  stop("For Bayesian mode: require mcmc_iter > burnin, thin >= 1, proposal_scale > 0")
}
if (args$method == "bayesian" && (args$mcmc_chains < 1 || args$adapt_burnin_frac < 0 || args$adapt_burnin_frac > 1 || args$adapt_interval < 1 || args$adapt_target_accept <= 0 || args$adapt_target_accept >= 1)) {
  stop("For Bayesian mode: require mcmc_chains >= 1, 0 <= adapt_burnin_frac <= 1, adapt_interval >= 1, and 0 < adapt_target_accept < 1")
}
if (args$method == "bayesian" && (args$local_move_prob < 0 || args$local_move_prob > 1 || args$active_eps <= 0 || args$mean_active_sources < 0)) {
  stop("For Bayesian mode: require 0 <= local_move_prob <= 1, active_eps > 0, and mean_active_sources >= 0")
}
if (args$method == "bayesian" && (!(args$hybrid_active_search %in% c(0, 1)) || args$max_active_sources < 0 || args$active_search_slots < 2 || args$active_search_iter < 10 || args$active_search_burnin < 0 || args$active_search_thin < 1 || args$active_search_jump_prob < 0 || args$active_search_jump_prob > 1)) {
  stop("For Bayesian mode: require hybrid_active_search in {0,1}, max_active_sources >= 0, active_search_slots >= 2, active_search_iter >= 10, active_search_burnin >= 0, active_search_thin >= 1, and 0 <= active_search_jump_prob <= 1")
}


if (!(args$cv_only %in% c(0L, 1L))) {
  stop("--cv_only must be 0 or 1")
}
if (args$cv_only == 1L && args$cv == "none") {
  stop("--cv_only needs a --cv spec")
}
if (args$cv == "loco" && args$method == "bayesian") {
  stop("--cv loco with the Bayesian estimator would run 22 MCMC refits per target; use evenodd, k<K> or test:<chroms>")
}
.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]))) else "."
source(file.path(.script_dir, "mixmodel_palette.R"))
if (args$cv != "none") {
  source(file.path(.script_dir, "mixmodel_cv.R"))
}


## --------------------------------------------------
## read input data

## The IBD tables are read in parallel, and only for the samples the fit needs.
##
##   1. Row prefilter. ibd_pop is consumed solely via
##      `filter(sample1 %in% target_samples)` / `%in% source_samples`, and
##      all_pops plus the matrix row space come from sample_map, not from
##      ibd_pop -- so dropping rows for excluded samples is LOSSLESS. It is applied
##      only when the panel is narrow enough to benefit.
##   2. Parallel read across the per-chromosome files, which only pays once (1) has
##      cut what the workers hand back: multisession must serialize the result to
##      the parent, and unfiltered that transfer eats the gain.
##
## readr is used rather than data.table::fread on purpose. The two parsers can
## differ by 1 ULP on a few values, which is numerically irrelevant but makes the
## MCMC redraw its whole trajectory, so a parser change would invalidate seeded
## results. With readr the prefiltered read is bit-identical to the unfiltered one,
## so seeded runs reproduce exactly. Wide panels (most of the samples) take the
## unchanged sequential path, where the read is a small part of the run time.
##
## Metadata is read BEFORE the IBD tables (needed to build the keep set), which
## also puts the whole read ahead of set.seed() so it can never perturb the
## seeded RNG state the model depends on.

read_ibd <- function(files, keep_ids, n_panel, threads, prefilter_frac) {
  frac <- if (is.null(keep_ids) || is.na(n_panel) || n_panel <= 0) {
    1
  } else {
    length(keep_ids) / n_panel
  }
  if (!(prefilter_frac > 0 && length(keep_ids) > 0 && frac <= prefilter_frac)) {
    cat(sprintf(
      "__ reading IBD data (sequential; panel uses %.0f%% of samples) __\n",
      100 * frac
    ))
    return(map_dfr(files, ~ read_tsv(.x, col_types = "ccccdi", show_col_types = FALSE)))
  }

  kf <- tempfile(fileext = ".ids")
  writeLines(keep_ids, kf)
  on.exit(unlink(kf), add = TRUE)
  nw <- max(1L, min(as.integer(threads), length(files)))
  cat(sprintf(
    "__ reading IBD data (prefilter to %d/%d samples = %.1f%%, %d workers) __\n",
    length(keep_ids), n_panel, 100 * frac, nw
  ))
  awk <- "NR==FNR{k[$1];next} FNR==1||($2 in k)"
  cmds <- vapply(files, function(f) sprintf(
    "zcat %s | awk -F'\t' %s %s -", shQuote(f), shQuote(awk), shQuote(kf)
  ), character(1), USE.NAMES = FALSE)

  old_plan <- future::plan()
  future::plan(future::multisession, workers = nw)
  on.exit(future::plan(old_plan), add = TRUE)
  dplyr::bind_rows(furrr::future_map(cmds, function(cm) {
    readr::read_tsv(pipe(cm), col_types = "ccccdi", show_col_types = FALSE)
  }))
}

cat("__ reading metadata __\n")
sample_map <- read_tsv(args$sample_file,
  show_col_types = FALSE
)

group_map <- read_tsv(args$group_file, show_col_types = FALSE)

## keep set = every sample this panel can use as a target or a source
.keep_ids <- character(0)
if (all(c("sample_id", "group") %in% colnames(group_map))) {
  .keep_ids <- unique(group_map$sample_id[group_map$group %in% c("target", "source")])
}

ibd_pop <- read_ibd(
  files = args$files,
  keep_ids = .keep_ids,
  n_panel = nrow(sample_map),
  threads = as.integer(args$threads),
  prefilter_frac = args$read_prefilter_frac
)

## Seed before plan(): furrr derives its parallel streams from the RNG state
## current at the future_map_dfr call, so seeding here makes both the hybrid
## source search and the MCMC reproducible across runs.
if (!is.null(args$seed) && !is.na(args$seed) && args$seed >= 0) {
  set.seed(as.integer(args$seed), kind = "L'Ecuyer-CMRG")
  cat(sprintf("__ seeded with %d __\n", as.integer(args$seed)))
}

plan(multisession, workers = as.integer(args$threads))
on.exit(plan(sequential), add = TRUE)


individuals <- read_tsv(args$individuals_file,
  show_col_types = FALSE
)

n_markers <- read_tsv(args$length_file,
  show_col_types = FALSE
)

# Keep chromosome keys type-stable across joins/bind_rows.
n_markers <- n_markers |>
  mutate(chrom = as.character(chrom))

req_ibd <- c("sample1", "pop_id1", "pop_id2", "chrom", "ibd")
if (!all(req_ibd %in% colnames(ibd_pop))) {
  stop("IBD input files must contain columns: ", paste(req_ibd, collapse = ", "))
}
if (!all(c("sample_id", "pop_id") %in% colnames(sample_map))) {
  stop("sample_file must contain columns: sample_id, pop_id")
}
if (!all(c("sample_id", "group") %in% colnames(group_map))) {
  stop("group_file must contain columns: sample_id, group")
}
if (!"sample_id" %in% colnames(individuals)) {
  stop("individuals file must contain column: sample_id")
}
if (!all(c("chrom", "n") %in% colnames(n_markers))) {
  stop("length_file must contain columns: chrom, n")
}
if (args$cv != "none") {
  # fail fast: a bad --cv spec must not surface only after the fit has finished
  cv_check_folds_present(
    cv_make_folds(args$cv, as.character(n_markers$chrom), n_markers$n),
    unique(ibd_pop$chrom)
  )
}
if (!"label" %in% colnames(individuals)) {
  individuals <- individuals |>
    mutate(label = sample_id)
}

group_map <- group_map |>
  select(sample_id, group) |>
  distinct()


## --------------------------------------------------
## set up helpers for models

sample_info <- sample_map |>
  left_join(group_map,
    by = c("sample_id" = "sample_id")
  )

label_map <- individuals |>
  select(sample_id, label)

## target samples to model
target_samples <- sample_info |>
  filter(
    group == "target",
    pop_id != "exclude"
  ) |>
  distinct(sample_id) |>
  pull(sample_id)

## source individuals
source_samples <- sample_info |>
  filter(group == "source") |>
  distinct(sample_id) |>
  pull(sample_id)

## Recipient (cluster_min_dist) samples carry a "_r" suffix on their pop_id so
## they can be separated in the output/plots. That suffix is a labeling device
## only: the donor palette and source populations must use the base cluster
## label, matching the aggregated IBD data (pop_id1/pop_id2 never carry "_r").
## Stripping it here keeps the model matrices on the base population set; the
## "_r" label is re-attached to the output via the sample_map join below.
source_pops <- sample_info |>
  filter(group == "source", pop_id != "exclude") |>
  mutate(pop_id = sub("_r$", "", pop_id)) |>
  distinct(pop_id) |>
  pull(pop_id)

## full set of populations
all_pops <- sample_info |>
  filter(pop_id != "exclude") |>
  mutate(pop_id = sub("_r$", "", pop_id)) |>
  distinct(pop_id) |>
  pull(pop_id)

if (length(target_samples) == 0) {
  stop("No target samples found in group_file")
}
if (length(source_samples) == 0) {
  stop("No source samples found in group_file")
}

## Row of each target's own cluster in the donor palette, for res_norm_ex_self().
## Both estimators build their palettes with rownames = all_pops (see
## get_sum_matrix), so one lookup serves both branches. The "_r" suffix is a
## labeling device only, as for source_pops above.
self_cl <- sample_info |>
  filter(sample_id %in% target_samples) |>
  mutate(cl = sub("_r$", "", pop_id)) |>
  distinct(sample_id, cl)
target_self_row <- setNames(match(self_cl$cl, all_pops), self_cl$sample_id)
## whether that cluster is itself a source: if so the model CAN fit the self
## column, and res_norm_ex_self means "fit away from home" rather than "fit on
## the part no source can reach" -- a different claim, so it is flagged, not
## silently folded in
target_self_src <- setNames(self_cl$cl %in% source_pops, self_cl$sample_id)

## --------------------------------------------------
## chromosome hold-out evaluation (--cv); folds, scoring and the per-fold driver
## live in mixmodel_cv.R. Everything here is inert unless --cv is set, and it runs
## after the main table is written (or instead of it with --cv_only), so it cannot
## change the genome-wide fit.

## top level on purpose: see cv_run_fold on why fit functions must not be closures
## over large frames
cv_bayes_fit <- function(y, S) {
  infer_sourcefind(
    y = y,
    source_mat = S,
    n_iter = args$mcmc_iter,
    burnin = args$burnin,
    thin = args$thin,
    proposal_scale = args$proposal_scale,
    n_chains = args$mcmc_chains,
    adapt_burnin_frac = args$adapt_burnin_frac,
    adapt_interval = args$adapt_interval,
    adapt_target_accept = args$adapt_target_accept,
    local_move_prob = args$local_move_prob,
    mean_active_sources = args$mean_active_sources,
    active_eps = args$active_eps,
    hybrid_active_search = as.logical(args$hybrid_active_search),
    max_active_sources = args$max_active_sources,
    active_search_slots = args$active_search_slots,
    active_search_iter = args$active_search_iter,
    active_search_burnin = args$active_search_burnin,
    active_search_thin = args$active_search_thin,
    active_search_jump_prob = args$active_search_jump_prob
  )$p
}

cv_out_path <- if (!is.null(args$cv_out)) {
  args$cv_out
} else {
  sub("(\\.tsv)?(\\.gz)?$", ".cv.tsv", args$out_file)
}

run_cv <- function() {
  cat(sprintf("__ chromosome hold-out evaluation: --cv %s (%s) __\n", args$cv, args$method))
  chroms <- as.character(n_markers$chrom)
  folds <- cv_make_folds(args$cv, chroms, n_markers$n)
  cv_check_folds_present(folds, unique(ibd_pop$chrom))
  for (fd in folds) {
    cat(sprintf(
      "   fold %-12s fit on %2d chr, score on %2d chr (%s)\n",
      fd$name, length(fd$train), length(fd$test), paste(fd$test, collapse = ",")
    ))
  }

  ## per-chromosome raw sums, same construction as the NNLS jackknife
  by_chrom_t <- map(setNames(chroms, chroms), function(i) {
    ibd_pop |>
      filter(sample1 %in% target_samples, chrom == i) |>
      get_sum_matrix(all_pops, sample1, pop_id2) |>
      align_matrix_cols(target_samples)
  })
  by_chrom_s <- map(setNames(chroms, chroms), function(i) {
    ibd_pop |>
      filter(sample1 %in% source_samples, chrom == i) |>
      get_sum_matrix(all_pops, pop_id1, pop_id2) |>
      align_matrix_cols(source_pops)
  })
  sum_over <- function(m, ch) reduce(m[ch], `+`)

  fit_fun <- if (args$method == "bayesian") cv_bayes_fit else cv_fit_nnls

  cv_tab <- map_dfr(folds, function(fd) {
    t0 <- Sys.time()
    res <- cv_run_fold(
      fold = fd,
      fit_fun = fit_fun,
      Ttr = normalize_matrix_cols(sum_over(by_chrom_t, fd$train)),
      Str = normalize_matrix_cols(sum_over(by_chrom_s, fd$train)),
      Tte_raw = sum_over(by_chrom_t, fd$test),
      Ste_raw = sum_over(by_chrom_s, fd$test),
      self_row = target_self_row,
      self_src = target_self_src,
      n_chunks = 4L * args$threads
    )
    cat(sprintf(
      "   fold %-12s done in %.1f min\n", fd$name,
      as.numeric(difftime(Sys.time(), t0, units = "mins"))
    ))
    res
  })

  cv_tab <- cv_tab |>
    left_join(distinct(select(sample_map, sample_id, pop_id), sample_id, .keep_all = TRUE),
      by = "sample_id"
    ) |>
    relocate(pop_id, .after = sample_id)

  cat("__ writing CV table: ", cv_out_path, " __\n", sep = "")
  dir.create(dirname(cv_out_path), showWarnings = FALSE, recursive = TRUE)
  write_tsv(cv_tab, file = cv_out_path)

  cv_tab |>
    group_by(fold) |>
    summarise(
      n = n(),
      n_na = sum(is.na(res_test_ex_self)),
      med_res_train = median(res_train_ex_self, na.rm = TRUE),
      med_res_test = median(res_test_ex_self, na.rm = TRUE),
      med_ll_test = median(ll_test, na.rm = TRUE),
      med_active = median(n_active),
      .groups = "drop"
    ) |>
    as.data.frame() |>
    print()
  invisible(cv_tab)
}

if (args$cv_only == 1L) {
  run_cv()
  cat("__ done! (cv only) __\n")
  quit(save = "no", status = 0)
}

## target sharing matrix (populations x target samples)
ibd_pop_target <- ibd_pop |>
  filter(sample1 %in% target_samples) |>
  get_palette_matrix(all_pops, sample1, pop_id2)

## source sharing matrix (populations x source populations)
ibd_pop_source <- ibd_pop |>
  filter(sample1 %in% source_samples) |>
  get_palette_matrix(all_pops, pop_id1, pop_id2)

if (!args$palette_scale %in% c("normalized", "raw")) stop("--palette_scale must be normalized or raw")
raw_scale <- identical(args$palette_scale, "raw")
if (raw_scale && args$cv != "none") stop("--palette_scale raw is not supported together with --cv; use --palette_scale normalized")
n_src_ind <- table(sub("_r$", "", sample_info$pop_id[sample_info$group %in% "source" & sample_info$pop_id != "exclude"]))
ibd_pop_source_fit <- if (raw_scale) {
  ibd_pop |>
    filter(sample1 %in% source_samples) |>
    get_sum_matrix(all_pops, pop_id1, pop_id2) |>
    raw_source_matrix(n_src_ind)
} else {
  ibd_pop_source
}

stopifnot(identical(colnames(ibd_pop_source_fit), colnames(ibd_pop_source)))
if (args$method == "bayesian") {
  cat("__ estimating Bayesian coefficients __\n")
  p_full <- future_map_dfr(target_samples, function(x) {
    fit <- infer_sourcefind(
      y = ibd_pop_target[, x],
      source_mat = ibd_pop_source_fit,
      se_n_copies = args$se_genome_length_cm,
      raw_scale = raw_scale,
      n_iter = args$mcmc_iter,
      burnin = args$burnin,
      thin = args$thin,
      proposal_scale = args$proposal_scale,
      n_chains = args$mcmc_chains,
      adapt_burnin_frac = args$adapt_burnin_frac,
      adapt_interval = args$adapt_interval,
      adapt_target_accept = args$adapt_target_accept,
      local_move_prob = args$local_move_prob,
      mean_active_sources = args$mean_active_sources,
      active_eps = args$active_eps,
      hybrid_active_search = as.logical(args$hybrid_active_search),
      max_active_sources = args$max_active_sources,
      active_search_slots = args$active_search_slots,
      active_search_iter = args$active_search_iter,
      active_search_burnin = args$active_search_burnin,
      active_search_thin = args$active_search_thin,
      active_search_jump_prob = args$active_search_jump_prob
    )

    tibble(
      sample_id = x,
      source_pop = colnames(ibd_pop_source),
      p = fit$p,
      p_median = fit$p_median,
      se = fit$se,
      res_norm = fit$res_norm,
      ## same quantity as res_norm here; carried so that res_norm_rmse means one
      ## thing across both estimators (see res_norm_ex_self header)
      res_norm_rmse = fit$res_norm,
      res_norm_ex_self = res_norm_ex_self(
        ibd_pop_target[, x],
        model_pred(ibd_pop_source_fit, fit$p, raw_scale),
        target_self_row[[x]]
      ),
      self_share = ifelse(is.na(target_self_row[[x]]),
        NA_real_,
        ibd_pop_target[target_self_row[[x]], x]
      ),
      self_is_source = target_self_src[[x]],
      accept_rate = fit$accept_rate,
      accept_rate_min = fit$accept_rate_min,
      accept_rate_max = fit$accept_rate_max,
      ess_min = fit$ess_min,
      ess_median = fit$ess_median,
      rhat_max = fit$rhat_max,
      rhat_median = fit$rhat_median,
      active_sources_median = fit$active_sources_median,
      selected_sources_n = fit$selected_sources_n,
      n_keep = fit$n_keep,
      n_chains = fit$n_chains,
      proposal_scale_final = fit$proposal_scale_final
    )
  }, .options = furrr_options(seed = TRUE))
} else {
  cat("__ estimating NNLS coefficients __\n")
  chroms <- as.character(n_markers$chrom)

  ## Precompute per-chromosome raw sharing matrices and derive leave-one-chrom matrices
  chrom_target_sum <- map(chroms, function(i) {
    ibd_pop |>
      filter(sample1 %in% target_samples, chrom == i) |>
      get_sum_matrix(all_pops, sample1, pop_id2) |>
      align_matrix_cols(target_samples)
  })
  names(chrom_target_sum) <- chroms

  chrom_source_sum <- map(chroms, function(i) {
    ibd_pop |>
      filter(sample1 %in% source_samples, chrom == i) |>
      get_sum_matrix(all_pops, pop_id1, pop_id2) |>
      align_matrix_cols(source_pops)
  })
  names(chrom_source_sum) <- chroms

  total_target_sum <- reduce(chrom_target_sum, `+`)
  total_source_sum <- reduce(chrom_source_sum, `+`)
  ibd_pop_target <- normalize_matrix_cols(total_target_sum)
  ibd_pop_source <- normalize_matrix_cols(total_source_sum)
  ## the NNLS branch aligns the source columns to source_pops; rebuild the raw matrix in the same order
  if (raw_scale) {
    ibd_pop_source_fit <- raw_source_matrix(total_source_sum, n_src_ind)
    stopifnot(identical(colnames(ibd_pop_source_fit), colnames(ibd_pop_source)),
              identical(rownames(ibd_pop_source_fit), rownames(ibd_pop_source)))
  } else {
    ## the matrix defined above has the columns in a different order; the weights follow source_pops
    ibd_pop_source_fit <- ibd_pop_source
  }

  p_genome <- future_map_dfr(target_samples, function(x) {
    r <- if (raw_scale) {
      pnnls(ibd_pop_source_fit, total_target_sum[, x])
    } else {
      pnnls(ibd_pop_source, ibd_pop_target[, x], sum = 1)
    }
    ## raw scale: free overall scale, weights normalised afterwards, residual reported on the proportion scale
    rn <- r$rnorm
    if (raw_scale) {
      r$x <- r$x / sum(r$x)
      rn <- sqrt(sum((model_pred(ibd_pop_source_fit, r$x, TRUE) - ibd_pop_target[, x])^2))
    }
      tibble(
        sample_id = x,
        source_pop = colnames(ibd_pop_source),
        p = r$x,
        res_norm = rn,
        ## pnnls returns an L2 norm, sqrt(length(y)) larger than the bayesian
        ## res_norm; this is the same figure on the bayesian scale
        res_norm_rmse = rn / sqrt(nrow(ibd_pop_target)),
        res_norm_ex_self = res_norm_ex_self(
          ibd_pop_target[, x],
          model_pred(ibd_pop_source_fit, r$x, raw_scale),
          target_self_row[[x]]
        ),
        self_share = ifelse(is.na(target_self_row[[x]]),
          NA_real_,
          ibd_pop_target[target_self_row[[x]], x]
        ),
        self_is_source = target_self_src[[x]],
        chrom = "0",
        accept_rate = NA_real_,
        accept_rate_min = NA_real_,
        accept_rate_max = NA_real_,
        ess_min = NA_real_,
        ess_median = NA_real_,
        rhat_max = NA_real_,
        rhat_median = NA_real_,
        active_sources_median = NA_real_,
        selected_sources_n = NA_real_,
        n_keep = NA_real_,
        n_chains = NA_real_,
        proposal_scale_final = NA_real_
      )
    })

  p_jk <- future_map_dfr(chroms, function(i) {
    ibd_pop_target_i <- normalize_matrix_cols(total_target_sum - chrom_target_sum[[i]])
    ibd_pop_source_i <- normalize_matrix_cols(total_source_sum - chrom_source_sum[[i]])

    if (raw_scale) ibd_pop_source_raw_i <- raw_source_matrix(total_source_sum - chrom_source_sum[[i]], n_src_ind)

    map_dfr(target_samples, function(x) {
      r <- if (raw_scale) {
        pnnls(ibd_pop_source_raw_i, (total_target_sum - chrom_target_sum[[i]])[, x])
      } else {
        pnnls(ibd_pop_source_i, ibd_pop_target_i[, x], sum = 1)
      }
      if (raw_scale) r$x <- r$x / sum(r$x)
      tibble(
        sample_id = x,
        source_pop = colnames(ibd_pop_source_i),
        p = r$x,
        res_norm = r$rnorm,
        res_norm_rmse = NA_real_,
        ## leave-one-chromosome fits feed the jackknife SE only; the reported
        ## residuals are all taken from chrom "0" below
        res_norm_ex_self = NA_real_,
        self_share = NA_real_,
        self_is_source = NA,
        chrom = as.character(i),
        accept_rate = NA_real_,
        accept_rate_min = NA_real_,
        accept_rate_max = NA_real_,
        ess_min = NA_real_,
        ess_median = NA_real_,
        rhat_max = NA_real_,
        rhat_median = NA_real_,
        active_sources_median = NA_real_,
        selected_sources_n = NA_real_,
        n_keep = NA_real_,
        n_chains = NA_real_,
        proposal_scale_final = NA_real_
      )
    })
  })

  p_full <- bind_rows(p_genome, p_jk) |>
    arrange(chrom) |>
    left_join(n_markers,
      by = c("chrom" = "chrom")
    ) |>
    group_by(sample_id, source_pop) |>
    summarise(
      se = se_jk(
        p[chrom == "0"],
        p[chrom != "0"],
        n[chrom != "0"]
      ),
      p = p[chrom == "0"],
      res_norm = res_norm[chrom == "0"],
      res_norm_rmse = res_norm_rmse[chrom == "0"],
      res_norm_ex_self = res_norm_ex_self[chrom == "0"],
      self_share = self_share[chrom == "0"],
      self_is_source = self_is_source[chrom == "0"],
      accept_rate = accept_rate[chrom == "0"],
      accept_rate_min = accept_rate_min[chrom == "0"],
      accept_rate_max = accept_rate_max[chrom == "0"],
      ess_min = ess_min[chrom == "0"],
      ess_median = ess_median[chrom == "0"],
      rhat_max = rhat_max[chrom == "0"],
      rhat_median = rhat_median[chrom == "0"],
      active_sources_median = active_sources_median[chrom == "0"],
      selected_sources_n = selected_sources_n[chrom == "0"],
      n_keep = n_keep[chrom == "0"],
      n_chains = n_chains[chrom == "0"],
      proposal_scale_final = proposal_scale_final[chrom == "0"],
      .groups = "drop"
    )
}

p_full <- p_full |>
  left_join(sample_map,
    by = c("sample_id" = "sample_id")
  ) |>
  left_join(label_map,
    by = c("sample_id" = "sample_id")
  ) |>
  mutate(group = "target") |>
  select(sample_id, label, pop_id, group, source_pop, p, any_of("p_median"), se, res_norm, res_norm_rmse, res_norm_ex_self, self_share, self_is_source, accept_rate, accept_rate_min, accept_rate_max, ess_min, ess_median, rhat_max, rhat_median, active_sources_median, selected_sources_n, n_keep, n_chains, proposal_scale_final)

p_source <- group_map |>
  filter(
    group == "source"
  ) |>
  select(sample_id, group) |>
  left_join(sample_map, by = "sample_id") |>
  left_join(label_map,
    by = c("sample_id" = "sample_id")
  ) |>
  mutate(
    source_pop = pop_id,
    p = 1,
    se = 0,
    res_norm = 0,
    res_norm_rmse = 0,
    res_norm_ex_self = 0,
    self_share = NA_real_,
    self_is_source = NA,
    accept_rate = NA_real_,
    accept_rate_min = NA_real_,
    accept_rate_max = NA_real_,
    ess_min = NA_real_,
    ess_median = NA_real_,
    rhat_max = NA_real_,
    rhat_median = NA_real_,
    active_sources_median = NA_real_,
    selected_sources_n = NA_real_,
    n_keep = NA_real_,
    n_chains = NA_real_,
    proposal_scale_final = NA_real_
  )

o <- bind_rows(p_full, p_source)

cat("__ writing output __\n")

write_tsv(o,
  file = args$out_file
)

if (args$cv != "none") {
  run_cv()
}

cat("__ done! __\n")
