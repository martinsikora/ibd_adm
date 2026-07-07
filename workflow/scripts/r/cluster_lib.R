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

## Shared helpers for the staged IBD hierarchical-clustering pipeline.
## Sourced by cluster_build_matrix.R / cluster_distance.R / cluster_hclust.R /
## cluster_cut.R / cluster_plot.R. Contains no top-level side effects.

## --------------------------------------------------
## stage 1: feature matrix

## Read the per-chromosome IBD tables and build the dense symmetric pairwise
## IBD matrix for `inds`, summing sharing across chromosomes. Uses data.table
## for a lower-memory, faster read than readr::read_tsv + pivot_wider.
## Equivalent to cluster_hc.R:176-198 (one undirected value per pair; diag 0).
read_ibd_matrix <- function(files, inds, threads = 1L) {
  inds_set <- unique(inds)
  dt <- data.table::rbindlist(lapply(files, function(f) {
    x <- data.table::fread(
      f,
      header = FALSE,
      col.names = c("sample1", "sample2", "chrom", "ibd"),
      colClasses = list(character = 1:3, numeric = 4L),
      nThread = threads
    )
    x <- x[sample1 %in% inds_set & sample2 %in% inds_set]
    x[, .(ibd = sum(ibd)), by = .(sample1, sample2)]
  }))
  ## sum across chromosomes
  dt <- dt[, .(ibd = sum(ibd)), by = .(sample1, sample2)]

  n <- length(inds)
  m <- matrix(0, nrow = n, ncol = n, dimnames = list(inds, inds))
  i <- match(dt$sample1, inds)
  j <- match(dt$sample2, inds)
  ## place value on both triangles so the matrix is symmetric regardless of
  ## which direction each pair was stored in
  m[cbind(i, j)] <- dt$ibd
  m[cbind(j, i)] <- dt$ibd
  diag(m) <- 0
  m
}

## --------------------------------------------------
## stage 2: feature transforms (cluster_hc.R:204-216)

## Column z-score (per partner sample) then per-sample row L2 normalization,
## each applied only if requested. Returns the transformed matrix.
apply_feature_transforms <- function(m, standardize_features = FALSE,
                                     normalize_ibd_vectors = FALSE) {
  if (isTRUE(standardize_features)) {
    col_means <- colMeans(m)
    col_sds <- apply(m, 2, sd)
    col_sds[col_sds == 0] <- 1
    m <- sweep(m, 2, col_means, "-")
    m <- sweep(m, 2, col_sds, "/")
  }
  if (isTRUE(normalize_ibd_vectors)) {
    row_norms <- sqrt(rowSums(m^2))
    row_norms[row_norms == 0] <- 1
    m <- sweep(m, 1, row_norms, "/")
  }
  m
}

## --------------------------------------------------
## stage 3: tree hierarchy strings (cluster_hc.R:256-288)

## Derive, from the hclust dendrogram, the per-sample cluster_id path string and
## the ancestral-prefix expansion. Depends only on the tree, so it is computed
## once in the hclust stage and reused by every cut. Returns list(cl_hier,
## cl_ids_expand).
build_tree_hier <- function(res_hc, inds) {
  res_hc_el <- res_hc |>
    as.dendrogram() |>
    data.tree::as.Node(name = "0") |>
    data.tree::ToDataFrameNetwork() |>
    tibble::as_tibble()

  cl_hier <- res_hc_el |>
    dplyr::mutate(sample_id = gsub(".*/", "", to)) |>
    dplyr::filter(sample_id %in% inds) |>
    dplyr::mutate(cluster_id = gsub("/", "_", from)) |>
    dplyr::select(sample_id, cluster_id) |>
    dplyr::arrange(cluster_id, sample_id)

  cl_ids <- cl_hier |>
    dplyr::distinct(cluster_id) |>
    dplyr::pull(cluster_id)

  cl_ids_expand <- purrr::map_dfr(cl_ids, function(k) {
    r <- strsplit(k, "_") |>
      unlist()
    r1 <- purrr::map_chr(seq_along(r), ~ paste(r[1:.x], collapse = "_"))
    tibble::tibble(
      cluster_id = k,
      cluster_id_anc = r1,
      cluster_level = seq_along(r1) - 1L
    )
  })

  list(cl_hier = cl_hier, cl_ids_expand = cl_ids_expand)
}

## --------------------------------------------------
## readable cluster labels (cluster_hc.R:352-394)

cluster_depth_from_id <- function(cluster_id) {
  tokens <- strsplit(cluster_id, "_") |>
    unlist()
  max(0L, length(tokens) - 1L)
}

cluster_bits_from_id <- function(cluster_id) {
  tokens <- strsplit(cluster_id, "_") |>
    unlist()
  if (length(tokens) <= 1) {
    return("")
  }
  branch_tokens <- tokens[-1]
  bits <- purrr::map_chr(branch_tokens, function(x) {
    if (x == "1") {
      return("0")
    }
    if (x == "2") {
      return("1")
    }
    ## fallback for unexpected labels
    if (suppressWarnings(!is.na(as.integer(x)))) {
      return(as.character(as.integer(x) %% 2L))
    }
    "0"
  })
  paste0(bits, collapse = "")
}

bits_to_octal_triplets <- function(bits) {
  if (nchar(bits) == 0) {
    return("0")
  }
  pad_n <- (3 - (nchar(bits) %% 3)) %% 3
  if (pad_n > 0) {
    bits <- paste0(bits, paste(rep("0", pad_n), collapse = ""))
  }
  starts <- seq(1, nchar(bits), by = 3)
  oct <- purrr::map_chr(starts, ~ {
    as.character(strtoi(substr(bits, .x, .x + 2), base = 2))
  })
  paste(oct, collapse = "_")
}

## Attach cluster_depth, cluster_label and cluster_parent_label to a cl_final
## table that already carries (sample_id, cut_height, cluster_id, cluster_level,
## group, ...). Encapsulates cluster_hc.R:396-463, including the explicit
## "unassigned" handling for samples cutreeDynamic left out (cluster_terminal 0).
add_cluster_labels <- function(cl_final) {
  cl_label_map <- cl_final |>
    dplyr::distinct(cut_height, cluster_id) |>
    dplyr::arrange(cut_height, cluster_id) |>
    dplyr::mutate(
      cluster_depth = dplyr::if_else(
        cluster_id == "unassigned",
        NA_integer_,
        purrr::map_int(cluster_id, cluster_depth_from_id)
      ),
      cluster_label = dplyr::if_else(
        cluster_id == "unassigned",
        "unassigned",
        purrr::map_chr(cluster_id, ~ {
          bits <- cluster_bits_from_id(.x)
          paste0("C", bits_to_octal_triplets(bits))
        })
      )
    ) |>
    dplyr::group_by(cut_height, cluster_label) |>
    dplyr::mutate(
      label_dup_i = dplyr::row_number(),
      label_dup_n = dplyr::n(),
      cluster_label = ifelse(
        label_dup_n > 1,
        paste0(cluster_label, "_", sprintf("%02d", label_dup_i)),
        cluster_label
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(cut_height, cluster_id, cluster_depth, cluster_label)

  height_levels <- cl_final |>
    dplyr::filter(cut_height >= 0) |>
    dplyr::distinct(cut_height) |>
    dplyr::arrange(cut_height) |>
    dplyr::pull(cut_height)

  if (length(height_levels) > 0) {
    cl_final <- cl_final |>
      dplyr::left_join(cl_label_map, by = c("cut_height", "cluster_id"))

    parent_cut_map <- tibble::tibble(cut_height = unique(cl_final$cut_height)) |>
      dplyr::mutate(
        parent_cut_height = purrr::map_dbl(cut_height, ~ {
          if (.x < 0) {
            return(min(height_levels))
          }
          h_next <- height_levels[height_levels > .x]
          if (length(h_next) == 0) {
            return(NA_real_)
          }
          min(h_next)
        })
      )

    parent_lookup <- cl_final |>
      dplyr::select(sample_id, cut_height, cluster_label) |>
      dplyr::rename(
        parent_cut_height = cut_height,
        cluster_parent_label = cluster_label
      )

    cl_final <- cl_final |>
      dplyr::left_join(parent_cut_map, by = "cut_height",
        relationship = "many-to-one") |>
      dplyr::left_join(parent_lookup,
        by = c("sample_id", "parent_cut_height")) |>
      dplyr::select(-parent_cut_height)
  } else {
    cl_final <- cl_final |>
      dplyr::left_join(cl_label_map, by = c("cut_height", "cluster_id")) |>
      dplyr::mutate(cluster_parent_label = NA_character_)
  }
  cl_final
}
