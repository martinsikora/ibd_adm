## Chromosome hold-out evaluation for mixmodel_ibd.R.
##
## Sourced by mixmodel_ibd.R ONLY when --cv is not "none", so a default run never
## executes any of this. Pure functions: no I/O, no globals of their own.
##
## Expects these to be in scope (defined in mixmodel_ibd.R): normalize_matrix_cols(),
## res_norm_ex_self(). Packages: dplyr, purrr, furrr, lsei.
##
## Scoring convention. A fold fits weights p on the TRAIN chromosomes, then predicts
## the target's palette on the held-out TEST chromosomes as S_test %*% p, where S_test
## is the SOURCE palettes rebuilt from the test chromosomes (never the train ones).
## Everything is scored on the retained donors, i.e. with the target's own-cluster row
## dropped and both vectors renormalised, exactly as res_norm_ex_self does.
##
## Two criteria, deliberately:
##   res_test_ex_self  RMSE on the normalised palette. Same scale as the in-sample
##                     figure, but it weights donors by squared mass, so the ~20
##                     nearest donors decide it.
##   ll_test           held-out multinomial log-likelihood per unit of IBD,
##                     sum(y_frac * log(pred)) -- the same objective as log_post in
##                     infer_sourcefind. Prefer this to rank models. pred is mixed
##                     with the mean TRAIN palette at weight CV_LL_EPS so a donor the
##                     model gives zero mass cannot send the score to -Inf.
##
## res_train_ex_self is reported for reference only. A one-chromosome (or half-genome)
## test palette is far sparser than the train one, so test minus train mostly measures
## fold size, not optimism. Compare models on the SAME fold, as a paired difference.
##
## Held-out chromosomes do NOT test what is missing from every chromosome: a source
## absent from the panel is absent from the test half too (use the residual
## diagnostic / sink test), and relatives among the sources share haplotypes on the
## test chromosomes, so targets with kin in the source set score optimistically.

CV_LL_EPS <- 1e-3
CV_ACTIVE_EPS <- 1e-4

cv_expand_chroms <- function(txt) {
  parts <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
  unlist(lapply(parts, function(z) {
    if (grepl("^[0-9]+-[0-9]+$", z)) {
      ab <- as.integer(strsplit(z, "-", fixed = TRUE)[[1]])
      if (ab[1] > ab[2]) stop("descending chromosome range in --cv: ", z)
      as.character(seq(ab[1], ab[2]))
    } else {
      z
    }
  }), use.names = FALSE)
}

## Fold specs (chroms and n must be aligned; n = markers per chromosome):
##   none       no folds
##   loco       one fold per chromosome (refits; the jackknife fits are not reused)
##   evenodd    train_even (fit even, score odd) and train_odd (fit odd, score even)
##   k<K>       K blocks balanced on marker count (greedy, deterministic); each block
##              held out once
##   test:1,3-5a single custom fold holding out the listed chromosomes
cv_make_folds <- function(spec, chroms, n = NULL) {
  chroms <- as.character(chroms)
  if (anyDuplicated(chroms)) stop("duplicate chromosomes in length file")
  mk <- function(name, test) {
    test <- chroms[chroms %in% test]
    train <- setdiff(chroms, test)
    if (length(test) < 1L || length(train) < 1L) {
      stop("fold '", name, "' needs at least one train and one test chromosome")
    }
    list(name = name, train = train, test = test)
  }
  if (identical(spec, "none")) {
    return(list())
  }
  if (identical(spec, "loco")) {
    return(lapply(chroms, function(i) mk(paste0("loco_", i), i)))
  }
  if (identical(spec, "evenodd")) {
    ci <- suppressWarnings(as.integer(chroms))
    if (anyNA(ci)) stop("--cv evenodd needs integer chromosome names, got: ", paste(chroms, collapse = ","))
    return(list(
      mk("train_even", chroms[ci %% 2L == 1L]),
      mk("train_odd", chroms[ci %% 2L == 0L])
    ))
  }
  if (grepl("^k[0-9]+$", spec)) {
    k <- as.integer(sub("^k", "", spec))
    if (k < 2L || k > length(chroms)) stop("--cv ", spec, ": K must be in 2..", length(chroms))
    if (is.null(n) || length(n) != length(chroms)) stop("--cv ", spec, " needs marker counts per chromosome")
    blocks <- vector("list", k)
    load <- numeric(k)
    for (i in order(-n, seq_along(n))) {
      j <- which.min(load)
      blocks[[j]] <- c(blocks[[j]], chroms[i])
      load[j] <- load[j] + n[i]
    }
    return(lapply(seq_len(k), function(j) mk(sprintf("%s_%d", spec, j), blocks[[j]])))
  }
  if (startsWith(spec, "test:")) {
    test <- cv_expand_chroms(sub("^test:", "", spec))
    bad <- setdiff(test, chroms)
    if (length(bad)) stop("--cv test: unknown chromosome(s): ", paste(bad, collapse = ","))
    return(list(mk("custom", test)))
  }
  stop("--cv must be none, loco, evenodd, k<K> or test:<chroms>, got: ", spec)
}

## a fold naming a chromosome with no data would silently become an all-zero matrix
cv_check_folds_present <- function(folds, present) {
  for (fd in folds) {
    miss <- setdiff(c(fd$train, fd$test), present)
    if (length(miss)) {
      stop(
        "fold '", fd$name, "' names chromosome(s) with no IBD input: ",
        paste(miss, collapse = ","), " (input files cover: ",
        paste(sort(present), collapse = ","), ")"
      )
    }
  }
  invisible(TRUE)
}

cv_ll_ex_self <- function(y_raw, pred, self_row, q, eps = CV_LL_EPS) {
  if (length(self_row) != 1L || is.na(self_row)) {
    return(NA_real_)
  }
  keep <- seq_along(y_raw) != self_row
  y <- y_raw[keep]
  pk <- pred[keep]
  qk <- q[keep]
  sy <- sum(y)
  sp <- sum(pk)
  sq <- sum(qk)
  if (sy <= 0 || sp <= 0 || sq <= 0) {
    return(NA_real_)
  }
  pk <- (1 - eps) * pk / sp + eps * qk / sq
  sum(y / sy * log(pmax(pk, 1e-12)))
}

## one target, one fold. p: weights over the columns of S_tr / S_te.
## y_tr: normalised train palette; y_te_raw: raw (unnormalised) test palette.
cv_score_target <- function(p, y_tr, y_te_raw, S_tr, S_te, src_zero, self_row, q) {
  pred_tr <- as.vector(S_tr %*% p)
  pred_te <- as.vector(S_te %*% p)
  sy <- sum(y_te_raw)
  y_te <- if (sy > 0) y_te_raw / sy else y_te_raw
  c(
    cm_test = if (is.na(self_row)) NA_real_ else sum(y_te_raw[-self_row]),
    res_train_ex_self = res_norm_ex_self(y_tr, pred_tr, self_row),
    res_test_ex_self = res_norm_ex_self(y_te, pred_te, self_row),
    ll_test = cv_ll_ex_self(y_te_raw, pred_te, self_row, q),
    n_active = sum(p > CV_ACTIVE_EPS),
    p_lost_test = sum(p[src_zero])
  )
}

cv_fit_nnls <- function(y, S) lsei::pnnls(S, y, sum = 1)$x

## Fit every target on the train matrices and score it on the test matrices.
##   Ttr      donors x targets, column-normalised, TRAIN chromosomes
##   Str      donors x sources, column-normalised, TRAIN chromosomes
##   Tte_raw  donors x targets, raw sums, TEST chromosomes
##   Ste_raw  donors x sources, raw sums, TEST chromosomes
##   fit_fun  function(y, S) -> weights (length ncol(S))
## self_row / self_src are named by target sample_id. Targets are sent to workers in
## contiguous chunks so each worker receives only its own columns.
cv_run_fold <- function(fold, fit_fun, Ttr, Str, Tte_raw, Ste_raw, self_row, self_src,
                        n_chunks = 1L) {
  ids <- colnames(Ttr)
  if (!identical(ids, colnames(Tte_raw))) stop("train/test target columns differ")
  if (!identical(colnames(Str), colnames(Ste_raw))) stop("train/test source columns differ")
  S_te <- normalize_matrix_cols(Ste_raw)
  src_zero <- colSums(Ste_raw) == 0
  q <- rowMeans(Ttr)
  n_chunks <- max(1L, min(as.integer(n_chunks), length(ids)))
  chunk_ix <- split(seq_along(ids), cut(seq_along(ids), n_chunks, labels = FALSE))
  chunks <- lapply(chunk_ix, function(ix) {
    list(
      ids = ids[ix],
      ytr = Ttr[, ix, drop = FALSE],
      yte = Tte_raw[, ix, drop = FALSE],
      self = unname(self_row[ids[ix]])
    )
  })

  ## The worker gets a minimal environment on purpose: a closure defined here would
  ## carry this whole frame (Ttr, Tte_raw, every chunk) to every worker. fit_fun must
  ## likewise be defined at top level, not inside a function holding large objects.
  worker <- function(ch) {
    out <- vapply(seq_along(ch$ids), function(j) {
      p <- fit_fun(ch$ytr[, j], Str)
      cv_score_target(p, ch$ytr[, j], ch$yte[, j], Str, S_te, src_zero, ch$self[j], q)
    }, numeric(6))
    dplyr::as_tibble(t(out)) |> dplyr::mutate(sample_id = ch$ids, .before = 1)
  }
  environment(worker) <- list2env(
    list(fit_fun = fit_fun, Str = Str, S_te = S_te, src_zero = src_zero, q = q),
    parent = globalenv()
  )
  rows <- furrr::future_map(chunks, worker,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )

  ## scalars first: inside mutate() a new column called `fold` would mask the argument
  fold_name <- fold$name
  n_tr <- length(fold$train)
  n_te <- length(fold$test)
  n_zero <- sum(src_zero)
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      fold = fold_name,
      n_chr_train = n_tr,
      n_chr_test = n_te,
      n_src_zero_test = n_zero,
      self_is_source = unname(self_src[sample_id]),
      .after = sample_id
    ) |>
    dplyr::select(
      sample_id, fold, n_chr_train, n_chr_test, cm_test, res_train_ex_self,
      res_test_ex_self, ll_test, n_active, p_lost_test, n_src_zero_test, self_is_source
    )
}
