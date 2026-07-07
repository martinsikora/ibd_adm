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
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(ggtree)
  library(ape)
  library(phytools)
  library(heatmap3)
  library(scico)
})


## --------------------------------------------------
## functions

## build square, symmetric TVD matrix from long-format table
build_tvd_matrix <- function(tvd) {
  pops <- sort(unique(c(tvd$pop_id1, tvd$pop_id2)))
  m1 <- tvd |>
    mutate(
      pop_id1 = factor(pop_id1, levels = pops),
      pop_id2 = factor(pop_id2, levels = pops)
    ) |>
    pivot_wider(
      names_from = pop_id2,
      values_from = tvd,
      values_fill = 0
    ) |>
    column_to_rownames("pop_id1") |>
    as.matrix()
  m1 <- m1[pops, pops]
  (m1 + t(m1)) / 2
}

## build midpoint-rooted NJ tree from a TVD matrix
build_nj_tree <- function(m1) {
  njs(m1) |>
    midpoint.root()
}

plot_tvd_nj <- function(tvd, m1, tr, color_map, plot_file) {
  ## set up and plot
  pal_c <- color_map$color
  names(pal_c) <- color_map$pop_id

  pal_f <- color_map$fill
  names(pal_f) <- color_map$pop_id

  pal_s <- color_map$shape
  names(pal_s) <- color_map$pop_id

  w <- tvd |>
    pull(pop_id1) |>
    unique() |>
    length() %/% 20 + 3

  pdf(plot_file,
    width = 5,
    height = w
  )
  p <- ggtree(tr,
    aes(
      x = x,
      y = y
    ),
    linewidth = 0.1,
  )
  print(p %<+% color_map +
    geom_tippoint(
      aes(
        colour = label,
        fill = label,
        shape = label
      ),
      size = 1,
      alpha = 1
    ) +
    geom_tiplab(
      size = 1.5
    ) +
    scale_color_manual(values = pal_c) +
    scale_fill_manual(values = pal_f) +
    scale_shape_manual(values = pal_s) +
    xlim(c(0, 1)) +
    theme_tree() +
    theme(legend.position = "none"))
  dev.off()
}

## heatmap of pairwise TVD, ordered by the NJ tree (same style as cluster_hc.R)
plot_tvd_heatmap <- function(m1, tr, color_map, plot_file) {
  ## convert NJ tree -> dendrogram for heatmap3
  ## (zero negative edges + force ultrametric; only ordering/topology matter)
  tr2 <- tr
  tr2$edge.length[tr2$edge.length < 0] <- 0
  tr2 <- force.ultrametric(tr2, method = "extend")
  hc <- as.hclust(tr2)

  ## align matrix rows/cols with leaf order (heatmap3 reorders by index, not name)
  m1 <- m1[hc$labels, hc$labels]
  dendro <- as.dendrogram(hc)

  ## colors
  pal_c <- colorRampPalette(c("black", scico(20, palette = "davos")))(1000)
  pal_side <- color_map$color
  names(pal_side) <- color_map$pop_id
  side_colors <- pal_side[colnames(m1)]

  w <- nrow(m1) %/% 30 + 6

  pdf(plot_file,
    width = w,
    height = w
  )
  heatmap3(m1,
    scale = "none",
    col = pal_c,
    bg = pal_c,
    symm = TRUE,
    ColSideLabs = "TVD",
    ColSideColors = side_colors,
    RowSideLabs = "",
    useRaster = FALSE,
    cexRow = 0.25,
    cexCol = 0.25,
    Rowv = dendro,
    Colv = dendro,
    margins = c(10, 10)
  )
  dev.off()
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in",
  action = "store",
  dest = "in_file",
  help = "TVD table"
)

parser$add_argument("-c", "--color_file",
  action = "store",
  dest = "color_file",
  help = "File with color mapping"
)

parser$add_argument("-p", "--plot_file",
  action = "store",
  dest = "plot_file",
  help = "Plot file"
)

parser$add_argument("-m", "--heatmap_file",
  action = "store",
  dest = "heatmap_file",
  default = NULL,
  help = "Heatmap plot file (pairwise TVD, ordered by NJ tree)"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

tvd <- read_tsv(args$in_file,
  col_types = "ccd"
)

color_map <- read_tsv(args$color_file,
  col_types = cols(.default = "c")
)
if (!"fill" %in% colnames(color_map)) {
  color_map$fill <- color_map$color
}
if (!"shape" %in% colnames(color_map)) {
  color_map$shape <- "21"
}
color_map$shape <- as.integer(color_map$shape)


## --------------------------------------------------
## plot tree

cat("__ generating plot __\n")
m1 <- build_tvd_matrix(tvd)
tr <- build_nj_tree(m1)

plot_tvd_nj(tvd, m1, tr, color_map, args$plot_file)

if (!is.null(args$heatmap_file) && nzchar(args$heatmap_file)) {
  cat("__ generating heatmap __\n")
  plot_tvd_heatmap(m1, tr, color_map, args$heatmap_file)
}

cat("__ done! __\n")
