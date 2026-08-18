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

## Stage 5: plot the dendrogram (with cut-height colouring) and the IBD heatmap
## for a single cut height. Kept separate so re-cutting/re-labelling never
## regenerates the large heatmap PDF. By default only the log10 heatmap is
## rendered; pass --full_heatmap to also render the raw-scale panel.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(scico)
  library(heatmap3)
})

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(sub("^--file=", "", f[1])) else "."
}
source(file.path(get_script_dir(), "cluster_lib.R"))

parser <- ArgumentParser()
parser$add_argument("--hc", dest = "hc_file", help = "hclust .rds (res_hc)")
parser$add_argument("--matrix", dest = "matrix_file",
  help = "feature-matrix .rds (m_raw)")
parser$add_argument("--clusters", dest = "clusters_file",
  help = "clusters.tsv from the cut stage")
parser$add_argument("-s", "--sample_file", dest = "sample_file",
  help = "Sample-to-group mapping table")
parser$add_argument("--out_cl", dest = "out_file_cl",
  help = "Output cluster dendrogram plot filename")
parser$add_argument("--out_hm", dest = "out_file_hm",
  help = "Output heatmap filename")
parser$add_argument("--height", dest = "height", type = "character",
  help = paste(
    "Cut height being plotted. Character, not numeric, so a gated cut's tag",
    "(e.g. 0.5g0.2) can be plotted too -- it is matched numerically when both",
    "sides parse as numbers, and as a string otherwise."
  ))
parser$add_argument("--full_heatmap", dest = "full_heatmap",
  action = "store_true", default = FALSE,
  help = "Also render the raw-scale heatmap panel (default: log10 only)")
args <- parser$parse_args()

cat("__ reading inputs __\n")
sample_label <- read_tsv(args$sample_file, col_types = "ccc")
inds <- sample_label |>
  filter(group != "exclude") |>
  pull(sample_id)
inds_cl_full <- sample_label |>
  filter(group == "cluster_full") |>
  pull(sample_id)

res_hc <- readRDS(args$hc_file)
m_raw <- readRDS(args$matrix_file)
## cut_height as character: a gated cut stamps a tag like "0.5g0.2", which
## col_double() would silently turn into NA and drop every row.
cl_final <- read_tsv(args$clusters_file, col_types = cols(.default = col_character()))

## tolerant match, mirroring make_agg_panel_from_clusters.py: numeric when both
## sides are numbers (so "0.50" still matches 0.5), string otherwise.
height_match <- function(val, target) {
  vn <- suppressWarnings(as.numeric(val))
  tn <- suppressWarnings(as.numeric(target))
  ifelse(!is.na(vn) & !is.na(tn), abs(vn - tn) < 1e-9, as.character(val) == as.character(target))
}

heights <- args$height

## --------------------------------------------------
## plot cluster hierarchies
cat("__ plotting clusters __\n")

g_plots_all <- res_hc |>
  as_tbl_graph() |>
  activate(nodes) |>
  mutate(
    label1 = paste(
      label, sample_label$label[match(label, sample_label$sample_id)],
      sep = " / "
    )
  )

w <- length(inds) %/% 50

pdf(args$out_file_cl, width = w + 7, height = w + 7)
walk(unique(heights), ~ {
  cl <- cl_final |>
    filter(height_match(cut_height, .x))

  cl_ids <- cl |>
    count(cluster_id, sort = T) |>
    pull(cluster_id)

  pal_c1 <- scico(length(cl_ids), palette = "batlow")
  names(pal_c1) <- cl_ids

  d_title <- tibble(x = 0, y = 0, label = paste("cut height", .x))

  g1 <- g_plots_all |>
    left_join(cl, by = c("label" = "sample_id"))

  l1 <- create_layout(g1, "dendrogram", height = height, circular = TRUE)

  p <- ggraph(l1)
  print(p +
    geom_edge_elbow(color = "grey40", lineend = "round", edge_width = 0.25) +
    geom_node_point(aes(filter = leaf, color = cluster_id)) +
    geom_node_text(
      aes(filter = leaf, angle = node_angle(x, y), label = label1),
      hjust = "outward", size = 1.5
    ) +
    geom_node_text(aes(label = label), data = d_title) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = pal_c1) +
    theme_void() +
    theme(
      plot.margin = unit(rep(100, 4), "points"),
      legend.position = "none",
    ))
})
dev.off()

## --------------------------------------------------
## plot heatmap
cat("__ plotting heatmap __\n")

w <- length(inds) %/% 30 + 6
pal_c <- colorRampPalette(c(rev(scico(20, palette = "davos")), "black"))(1000)
cl_dendro <- res_hc |>
  as.dendrogram()

m1 <- m_raw[inds_cl_full, inds_cl_full]
lab <- paste(inds_cl_full,
  sample_label$label[match(inds_cl_full, sample_label$sample_id)],
  sep = " / "
)
colnames(m1) <- lab
rownames(m1) <- lab

m_l <- log10(m1)
m_l[m_l == -Inf] <- NA

pdf(file = args$out_file_hm, width = w, height = w)
walk(unique(heights), ~ {
  cl <- cl_final |>
    filter(height_match(cut_height, .x), group == "cluster_full")

  cl_ids <- cl |>
    count(cluster_id, sort = T) |>
    pull(cluster_id)

  pal_c1 <- scico(length(cl_ids), palette = "batlow")
  names(pal_c1) <- cl_ids

  pal_c2 <- pal_c1[cl$cluster_id]
  names(pal_c2) <- cl$sample_id

  ## raw-scale panel only when explicitly requested (halves default PDF size)
  if (isTRUE(args$full_heatmap)) {
    heatmap3(m1,
      scale = "none", col = pal_c, bg = pal_c, symm = TRUE,
      ColSideLabs = paste("cluster height ", .x),
      ColSideColors = pal_c2[inds_cl_full], RowSideLabs = "",
      useRaster = TRUE, cexRow = 0.25, cexCol = 0.25,
      Rowv = cl_dendro, Colv = cl_dendro, margins = c(10, 10)
    )
  }
  heatmap3(m_l,
    scale = "none", col = pal_c, bg = pal_c,
    ColSideLabs = paste("cluster height ", .x),
    ColSideColors = pal_c2[inds_cl_full], RowSideLabs = "",
    useRaster = TRUE, cexRow = 0.25, cexCol = 0.25,
    Rowv = cl_dendro, Colv = cl_dendro, margins = c(10, 10)
  )
})
dev.off()

cat("__ done! __\n")
