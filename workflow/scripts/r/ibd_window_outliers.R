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

## --------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(argparse)
})


## --------------------------------------------------
## functions

kl_div <- function(p, q, eps = 1e-12) {
  p1 <- p + eps
  q1 <- q + eps
  sum(p1 * log(p1 / q1))
}

js_div <- function(p, q, eps = 1e-12) {
  m <- (p + q) / 2
  0.5 * kl_div(p, m, eps = eps) + 0.5 * kl_div(q, m, eps = eps)
}

chi2_dist <- function(p, q, eps = 1e-12) {
  sum((p - q)^2 / (q + eps))
}

mad_z_score <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_raw <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(mad_raw) || mad_raw == 0) {
    return(rep(0, length(x)))
  }
  0.67448975 * (x - med) / mad_raw
}

score_windows <- function(m, method = "jsd", eps = 1e-12) {
  row_sums <- rowSums(m)
  valid <- row_sums > 0
  p_mat <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  p_mat[valid, ] <- m[valid, , drop = FALSE] / row_sums[valid]

  q <- colSums(m)
  q <- q / sum(q)

  if (method == "jsd") {
    score <- apply(p_mat, 1, function(p) js_div(p, q, eps = eps))
  } else if (method == "kl") {
    score <- apply(p_mat, 1, function(p) kl_div(p, q, eps = eps))
  } else if (method == "chi2") {
    score <- apply(p_mat, 1, function(p) chi2_dist(p, q, eps = eps))
  } else if (method == "max_mad_z") {
    z_mat <- apply(m, 2, mad_z_score)
    score <- apply(abs(z_mat), 1, max, na.rm = TRUE)
  } else {
    stop("Unsupported method: ", method)
  }

  tibble(
    score = score,
    window_total = row_sums,
    n_pop_nonzero = rowSums(m > 0)
  )
}

contribution_matrix <- function(m, method = "jsd", eps = 1e-12) {
  n_row <- nrow(m)
  n_col <- ncol(m)

  row_sums <- rowSums(m)
  valid <- row_sums > 0
  p_mat <- matrix(0, nrow = n_row, ncol = n_col)
  p_mat[valid, ] <- m[valid, , drop = FALSE] / row_sums[valid]

  q <- colSums(m)
  q <- q / sum(q)
  q_mat <- matrix(q, nrow = n_row, ncol = n_col, byrow = TRUE)

  if (method == "jsd") {
    mix <- (p_mat + q_mat) / 2
    return(0.5 * p_mat * log((p_mat + eps) / (mix + eps)) +
      0.5 * q_mat * log((q_mat + eps) / (mix + eps)))
  }
  if (method == "kl") {
    return(p_mat * log((p_mat + eps) / (q_mat + eps)))
  }
  if (method == "chi2") {
    return((p_mat - q_mat)^2 / (q_mat + eps))
  }
  if (method == "max_mad_z") {
    z_cols <- lapply(seq_len(n_col), function(j) mad_z_score(m[, j]))
    z_mat <- do.call(cbind, z_cols)
    return(abs(z_mat))
  }

  stop("Unsupported method: ", method)
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
  nargs = "*",
  help = "Optional list of input windowed coverage tables"
)

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "Input windowed coverage table (chromosome, pop_id, window_start, window_end, avg_n_seg)"
)

parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output file with outlier windows"
)

parser$add_argument("--score_file",
  action = "store",
  dest = "score_file",
  default = NULL,
  help = "Optional output file with scores for all windows"
)

parser$add_argument("--contrib_file",
  action = "store",
  dest = "contrib_file",
  default = NULL,
  help = "Optional output file with top contributing clusters per outlier window"
)

parser$add_argument("--method",
  action = "store",
  dest = "method",
  default = "jsd",
  help = "Outlier score: jsd, kl, chi2, max_mad_z [default %(default)s]"
)

parser$add_argument("--tail_prob",
  action = "store",
  dest = "tail_prob",
  type = "double",
  default = 0.01,
  help = "Upper-tail quantile used as outlier threshold [default %(default)s]"
)

parser$add_argument("--epsilon",
  action = "store",
  dest = "epsilon",
  type = "double",
  default = 1e-12,
  help = "Small positive value for stable divergence calculations [default %(default)s]"
)

parser$add_argument("--threshold_scope",
  action = "store",
  dest = "threshold_scope",
  default = "genome",
  help = "Threshold scope: genome or chrom [default %(default)s]"
)

parser$add_argument("--min_window_total",
  action = "store",
  dest = "min_window_total",
  type = "double",
  default = 0,
  help = "Minimum window_total required to call outliers [default %(default)s]"
)

parser$add_argument("--min_cluster_n",
  action = "store",
  dest = "min_cluster_n",
  type = "integer",
  default = 4L,
  help = "Minimum number of individuals required per cluster/pop_id [default %(default)s]"
)

parser$add_argument("--contrib_top_n",
  action = "store",
  dest = "contrib_top_n",
  type = "integer",
  default = 5L,
  help = "Number of top contributing clusters to keep per outlier window [default %(default)s]"
)

args <- parser$parse_args()

if (!(args$method %in% c("jsd", "kl", "chi2", "max_mad_z"))) {
  stop("--method must be one of: jsd, kl, chi2, max_mad_z")
}
if (args$tail_prob <= 0 || args$tail_prob >= 1) {
  stop("--tail_prob must be in the open interval (0, 1)")
}
if (!(args$threshold_scope %in% c("genome", "chrom"))) {
  stop("--threshold_scope must be one of: genome, chrom")
}
if (args$min_cluster_n < 1) {
  stop("--min_cluster_n must be >= 1")
}
if (args$contrib_top_n < 1) {
  stop("--contrib_top_n must be >= 1")
}


## --------------------------------------------------
## read and score windows

cat("__ reading windowed data __\n")
input_files <- args$files
if (length(input_files) == 0 && !is.null(args$in_file)) {
  input_files <- c(args$in_file)
}
if (length(input_files) == 0) {
  stop("Provide input with positional files and/or --in_file")
}

d <- map_dfr(input_files, ~ read_tsv(.x, show_col_types = FALSE))

req_cols <- c("chromosome", "pop_id", "window_start", "window_end", "avg_n_seg")
missing_cols <- setdiff(req_cols, colnames(d))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

d1 <- d |>
  transmute(
    chromosome = as.character(chromosome),
    pop_id = as.character(pop_id),
    window_start = as.integer(window_start),
    window_end = as.integer(window_end),
    avg_n_seg = as.double(avg_n_seg)
  ) |>
  filter(
    !is.na(chromosome),
    !is.na(pop_id),
    !is.na(window_start),
    !is.na(window_end),
    !is.na(avg_n_seg)
  )

if ("norm_factor" %in% colnames(d)) {
  if ("group_n" %in% colnames(d)) {
    pop_n <- d |>
      transmute(
        pop_id = as.character(pop_id),
        group_n = suppressWarnings(as.integer(group_n))
      ) |>
      filter(!is.na(pop_id), !is.na(group_n), group_n > 0) |>
      group_by(pop_id) |>
      summarise(cluster_n = max(group_n), .groups = "drop")
  } else if ("norm_unit" %in% colnames(d) && any(as.character(d$norm_unit) == "n_pairs", na.rm = TRUE)) {
    pop_n <- d |>
      transmute(
        pop_id = as.character(pop_id),
        norm_factor = suppressWarnings(as.double(norm_factor))
      ) |>
      filter(!is.na(pop_id), !is.na(norm_factor), norm_factor > 0) |>
      group_by(pop_id) |>
      summarise(
        cluster_n = as.integer(round((1 + sqrt(1 + 8 * max(norm_factor))) / 2)),
        .groups = "drop"
      )
  } else {
    pop_n <- d |>
      transmute(
        pop_id = as.character(pop_id),
        norm_factor = suppressWarnings(as.integer(norm_factor))
      ) |>
      filter(!is.na(pop_id), !is.na(norm_factor), norm_factor > 0) |>
      group_by(pop_id) |>
      summarise(cluster_n = max(norm_factor), .groups = "drop")
  }
} else {
  warning("Input has no norm_factor column; assuming cluster_n = 1 for all clusters")
  pop_n <- d1 |>
    distinct(pop_id) |>
    mutate(cluster_n = 1L)
}

eligible_pops <- pop_n |>
  filter(cluster_n >= args$min_cluster_n) |>
  pull(pop_id)

d1 <- d1 |>
  filter(pop_id %in% eligible_pops)

if (nrow(d1) == 0) {
  stop("No clusters pass --min_cluster_n filter")
}

n_clusters_used <- d1 |>
  pull(pop_id) |>
  unique() |>
  length()

cat("__ reshaping matrix __\n")
d_wide <- d1 |>
  pivot_wider(
    names_from = pop_id,
    values_from = avg_n_seg,
    values_fill = 0
  ) |>
  arrange(chromosome, window_start, window_end)

meta <- d_wide |>
  mutate(window_id = row_number()) |>
  select(window_id, chromosome, window_start, window_end)

m <- d_wide |>
  select(-chromosome, -window_start, -window_end) |>
  as.matrix()

cat("__ scoring windows __\n")
scores <- score_windows(m, method = args$method, eps = args$epsilon)
o_all <- bind_cols(meta, scores) |>
  mutate(
    method = args$method,
    tail_prob = args$tail_prob,
    threshold_scope = args$threshold_scope,
    min_window_total = args$min_window_total,
    min_cluster_n = args$min_cluster_n,
    n_clusters_used = n_clusters_used,
    is_eligible = window_total >= args$min_window_total
  )

if (!any(o_all$is_eligible)) {
  stop("No windows pass --min_window_total filter")
}

if (args$threshold_scope == "genome") {
  th <- quantile(
    o_all$score[o_all$is_eligible],
    probs = 1 - args$tail_prob,
    na.rm = TRUE
  )

  o_all <- o_all |>
    mutate(
      threshold = as.double(th),
      is_outlier = is_eligible & (score >= threshold)
    )
} else {
  th_by_chr <- o_all |>
    filter(is_eligible) |>
    group_by(chromosome) |>
    summarise(
      threshold = as.double(quantile(score, probs = 1 - args$tail_prob, na.rm = TRUE)),
      .groups = "drop"
    )

  o_all <- o_all |>
    left_join(th_by_chr, by = "chromosome") |>
    mutate(
      is_outlier = is_eligible & !is.na(threshold) & (score >= threshold)
    )
}

o_all <- o_all |>
  arrange(desc(score), chromosome, window_start)

o_outlier <- o_all |>
  filter(is_outlier)


## --------------------------------------------------
## write output

cat("__ writing outputs __\n")
write_tsv(o_outlier, args$out_file)

if (!is.null(args$score_file)) {
  write_tsv(o_all, args$score_file)
}

if (!is.null(args$contrib_file)) {
  if (nrow(o_outlier) == 0) {
    contrib_empty <- tibble(
      window_id = integer(),
      chromosome = character(),
      window_start = integer(),
      window_end = integer(),
      method = character(),
      score = double(),
      threshold = double(),
      pop_id = character(),
      cluster_n = integer(),
      contribution = double(),
      contribution_pos = double(),
      contribution_frac = double(),
      contribution_rank = integer()
    )
    write_tsv(contrib_empty, args$contrib_file)
  } else {
    contrib_mat <- contribution_matrix(m, method = args$method, eps = args$epsilon)
    colnames(contrib_mat) <- colnames(m)

    contrib_long <- as_tibble(contrib_mat) |>
      mutate(window_id = seq_len(nrow(contrib_mat))) |>
      pivot_longer(
        cols = -window_id,
        names_to = "pop_id",
        values_to = "contribution"
      ) |>
      left_join(pop_n, by = "pop_id") |>
      mutate(cluster_n = if_else(is.na(cluster_n), 1L, as.integer(cluster_n)))

    outlier_meta <- o_outlier |>
      select(window_id, chromosome, window_start, window_end, method, score, threshold)

    contrib_outlier <- contrib_long |>
      inner_join(outlier_meta, by = "window_id") |>
      mutate(contribution_pos = pmax(contribution, 0)) |>
      group_by(window_id) |>
      mutate(
        contrib_pos_sum = sum(contribution_pos, na.rm = TRUE),
        contribution_frac = if_else(
          contrib_pos_sum > 0,
          contribution_pos / contrib_pos_sum,
          NA_real_
        ),
        contribution_rank = row_number(desc(contribution_pos))
      ) |>
      ungroup() |>
      filter(contribution_rank <= args$contrib_top_n) |>
      arrange(window_id, contribution_rank) |>
      select(
        window_id, chromosome, window_start, window_end, method, score, threshold,
        pop_id, cluster_n, contribution, contribution_pos, contribution_frac, contribution_rank
      )

    write_tsv(contrib_outlier, args$contrib_file)
  }
}

cat("__ done! __\n")
