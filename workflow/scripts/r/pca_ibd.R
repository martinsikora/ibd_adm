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
  library(ggplot2)
  library(tibble)
})

clr_transform <- function(x, pseudocount = 1e-8) {
  x1 <- x + pseudocount
  log_x <- log(x1)
  sweep(log_x, 1, rowMeans(log_x), "-")
}

apply_transform <- function(x, method) {
  if (method == "hellinger") {
    return(sqrt(x))
  }
  if (method == "clr") {
    return(clr_transform(x))
  }
  if (method == "raw") {
    return(x)
  }
  stop("Unknown transform: ", method)
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
  nargs = "+",
  help = "Files with ibd sharing data for each chromosome"
)

parser$add_argument("-c", "--color_file",
  action = "store",
  dest = "color_file",
  help = "File with color mapping"
)

parser$add_argument("-s", "--sample_file",
  action = "store",
  dest = "sample_file",
  help = "File with sample to population mapping"
)

parser$add_argument("-g", "--group_file",
  action = "store",
  dest = "group_file",
  default = NULL,
  help = "File with sample to group mapping"
)

parser$add_argument("-o", "--out",
  action = "store",
  dest = "out_file",
  help = "output file"
)

parser$add_argument("-p", "--plot_file",
  action = "store",
  dest = "plot_file",
  help = "Plot file"
)

parser$add_argument("-n", "--number_pcs",
  action = "store",
  dest = "n_pcs",
  type = "integer",
  default = 20L,
  help = "Number of pcs to plot [default %(default)s]"
)

parser$add_argument("--project",
  action = "store_true",
  dest = "project",
  default = FALSE,
  help = "Project target samples on PCA from source groups"
)

parser$add_argument("--transform",
  action = "store",
  dest = "transform",
  default = "hellinger",
  help = "Feature transform: hellinger, clr, raw [default %(default)s]"
)

args <- parser$parse_args()

valid_transforms <- c("hellinger", "clr", "raw")
if (!(args$transform %in% valid_transforms)) {
  stop("--transform must be one of: ", paste(valid_transforms, collapse = ", "))
}

if (args$project && is.null(args$group_file)) {
  stop("--group_file is required when --project is used")
}


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

ibd_pop <- map_dfr(args$files, ~ {
  r <- read_tsv(.x,
    col_types = "ccccdi"
  )
  r
})

sample_map <- read_tsv(args$sample_file,
  col_types = "cc"
)

color_map <- read_tsv(args$color_file,
  col_types = "ccci"
)

group_map <- NULL
if (!is.null(args$group_file)) {
  group_map <- read_tsv(args$group_file,
    col_types = "cc"
  )
}


## --------------------------------------------------
## calculate PCs

cat("__ calculating PCA __\n")

## aggregate IBD sharing across chromosomes
d <- ibd_pop |>
  group_by(sample1, pop_id2) |>
  summarise(ibd = sum(ibd), .groups = "drop_last") |>
  mutate(p_ibd = ibd / sum(ibd)) |>
  ungroup()


## reshape into matrix
m <- d |>
  filter(sample1 %in% sample_map$sample_id) |>
  select(-ibd) |>
  pivot_wider(
    names_from = pop_id2,
    values_from = p_ibd,
    values_fill = 0
  ) |>
  column_to_rownames("sample1") |>
  as.matrix() |>
  t()

## samples x features matrix for PCA
x <- t(m)
x_t <- apply_transform(x, args$transform)

if (args$project) {
  ## split matrix in ref and project
  samples_ref <- group_map |>
    filter(group == "source") |>
    pull(sample_id)

  samples_proj <- group_map |>
    filter(group == "target") |>
    pull(sample_id)

  ## keep only available samples
  samples_ref <- intersect(samples_ref, rownames(x_t))
  samples_proj <- intersect(samples_proj, rownames(x_t))

  ## PCA on sources, then project targets
  pca_fit <- prcomp(x_t[samples_ref, , drop = FALSE],
    center = TRUE,
    scale. = FALSE
  )

  ref_scores <- pca_fit$x
  proj_scores <- predict(pca_fit, newdata = x_t[samples_proj, , drop = FALSE])
  pc_names <- paste("PC", seq_len(ncol(ref_scores)), sep = "")
  colnames(ref_scores) <- pc_names
  colnames(proj_scores) <- pc_names

  ## results table
  pca_ref <- ref_scores |>
    as_tibble() |>
    mutate(
      sample_id = samples_ref,
      group = "source"
    ) |>
    left_join(sample_map, by = "sample_id") |>
    select(sample_id, pop_id, group, everything())

  pca_proj <- proj_scores |>
    as_tibble() |>
    mutate(
      sample_id = samples_proj,
      group = "target"
    ) |>
    left_join(sample_map, by = "sample_id") |>
    select(sample_id, pop_id, group, everything())

  pca_res <- bind_rows(pca_ref, pca_proj)
  pca_sdev <- pca_fit$sdev
} else {
  ## Compute the PC basis on cluster_full (donor_recipient) samples only, then
  ## project cluster_min_dist (recipient) samples onto it, so lower-confidence
  ## min_dist samples do not shape the axes. Recipients are flagged by the "_r"
  ## suffix their pop_id carries in the mixmodel sample_map.
  is_recip <- grepl("_r$", sample_map$pop_id)
  ref_ids <- intersect(sample_map$sample_id[!is_recip], rownames(x_t))
  proj_ids <- intersect(sample_map$sample_id[is_recip], rownames(x_t))

  pca_fit <- prcomp(x_t[ref_ids, , drop = FALSE],
    center = TRUE,
    scale. = FALSE
  )

  ref_scores <- pca_fit$x
  colnames(ref_scores) <- paste("PC", seq_len(ncol(ref_scores)), sep = "")
  scores <- ref_scores |>
    as_tibble() |>
    mutate(sample_id = ref_ids)

  if (length(proj_ids) > 0) {
    proj_scores <- predict(pca_fit, newdata = x_t[proj_ids, , drop = FALSE])
    colnames(proj_scores) <- colnames(ref_scores)
    scores <- bind_rows(
      scores,
      proj_scores |>
        as_tibble() |>
        mutate(sample_id = proj_ids)
    )
  }

  ## results table
  pca_res <- scores |>
    left_join(sample_map, by = "sample_id") |>
    mutate(group = "target") |>
    select(sample_id, pop_id, group, everything())
  pca_sdev <- pca_fit$sdev
}

write_tsv(pca_res,
  file = args$out_file,
  col_names = TRUE
)


## --------------------------------------------------
## set up plot

cat("__ preparing plot data __\n")

## helpers
pal_c <- color_map$color
names(pal_c) <- color_map$pop_id

pal_f <- color_map$fill
names(pal_f) <- color_map$pop_id

pal_s <- color_map$shape
names(pal_s) <- color_map$pop_id

th <- theme_bw() +
  theme(
    panel.grid.major = element_line(
      linewidth = 0.25,
      linetype = "dotted"
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.key.size = unit(0.0015, "npc"),
    legend.text = element_text(size = 6)
  )

var_explained <- pca_sdev^2 / sum(pca_sdev^2) * 100

n_pcs <- min(args$n_pcs, ncol(pca_res) - 3)
if (n_pcs %% 2 != 0) {
  n_pcs <- n_pcs - 1L
}
n_pcs <- max(n_pcs, 2L)

labs <- matrix(
  paste("PC", 1:n_pcs,
    " (",
    formatC(var_explained[1:n_pcs], format = "f", digits = 1),
    "%)",
    sep = ""
  ),
  nrow = 2
)
pcs <- matrix(paste("PC", 1:n_pcs, sep = ""),
  nrow = 2
)

## Recipients keep their "_r"-suffixed pop_id in the output table; strip it here
## so cluster_min_dist points are colored by their parent cluster.
d <- pca_res |>
  mutate(
    pop_id = sub("_r$", "", pop_id),
    pop_id = factor(pop_id, levels = color_map$pop_id)
  )

d_source <- d |>
  filter(group == "source")

d_target <- d |>
  filter(group == "target")

## plot
cat("__ generating plot __\n")

pdf(args$plot_file,
  width = 6,
  height = 6
)

walk(seq_len(ncol(pcs)), ~ {
  p <- ggplot(d, aes(
    x = !!sym(pcs[1, .x]),
    y = !!sym(pcs[2, .x]),
    color = pop_id,
    fill = pop_id,
    shape = pop_id
  ))
  print(p +
    geom_point(
      size = 2.5,
      alpha = 1,
      stroke = 1,
      data = d_source,
      show.legend = FALSE
    ) +
    geom_point(
      size = 1,
      alpha = 0.8,
      data = d_target,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal_c) +
    scale_fill_manual(values = pal_f) +
    scale_shape_manual(values = pal_s) +
    xlab(labs[1, .x]) +
    ylab(labs[2, .x]) +
    guides(color = guide_legend(
      ncol = 4,
      override.aes = list(size = 2)
    )) +
    th)
})
dev.off()

cat("__ done! __\n")
