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

## Standalone colour legend for the mixture-model source clusters: one swatch
## per source population, coloured with the same palette as the barplots, and
## labelled with the cluster id plus its dominant sample label. Written as an
## extra plot file alongside the mixmodel barplot / grid.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

parser <- ArgumentParser()
parser$add_argument("-i", "--in", dest = "in_file", help = "mixmodel results tsv")
parser$add_argument("-c", "--color_file", dest = "color_file", help = "colour-map tsv")
parser$add_argument("-o", "--out", dest = "out_file", help = "output legend pdf")
args <- parser$parse_args()

cat("__ reading data __\n")
d <- read_tsv(args$in_file, col_types = cols(.default = "c"))
cm <- read_tsv(args$color_file, col_types = cols(.default = "c"))
if (!"fill" %in% colnames(cm)) {
  cm$fill <- cm$color
}
if (!"shape" %in% colnames(cm)) {
  cm$shape <- NA_character_
}

## source clusters = the populations that act as sources; dominant sample label
src <- d |> filter(group == "source")
if (nrow(src) == 0) {
  src <- d |>
    distinct(source_pop) |>
    transmute(pop_id = source_pop, label = NA_character_)
}
dom <- src |>
  filter(!is.na(pop_id), pop_id != "") |>
  count(pop_id, label) |>
  group_by(pop_id) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(pop_id, label)

## order like the barplots (colour-map pop_id order), attach colour + shape
pal_f <- setNames(cm$fill, cm$pop_id)
pal_c <- setNames(cm$color, cm$pop_id)
pal_s <- setNames(suppressWarnings(as.integer(cm$shape)), cm$pop_id)
dom <- dom |>
  mutate(ord = match(pop_id, cm$pop_id)) |>
  arrange(ord, pop_id) |>
  mutate(
    fill = ifelse(is.na(pal_f[pop_id]), "grey80", pal_f[pop_id]),
    color = ifelse(is.na(pal_c[pop_id]), "grey30", pal_c[pop_id]),
    shape = ifelse(is.na(pal_s[pop_id]), 22L, pal_s[pop_id]),
    lab = ifelse(is.na(label) | label == "", pop_id, paste0(pop_id, "  |  ", label))
  )

n <- nrow(dom)
maxchar <- max(nchar(dom$lab), 1L)
dom$y <- rev(seq_len(n))

cat("__ generating legend (", n, " sources) __\n", sep = "")

p <- ggplot(dom, aes(x = 0, y = y)) +
  ## colour swatch
  geom_tile(aes(fill = fill, color = color), width = 0.55, height = 0.8, linewidth = 0.3) +
  ## cluster shape marker (same shape codes as the PCA / auto legend)
  geom_point(aes(shape = shape, color = color, fill = fill),
    x = 0.55, size = 3.6, stroke = 0.7) +
  geom_text(aes(label = lab), x = 0.95, hjust = 0, size = 3) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_shape_identity() +
  coord_cartesian(
    xlim = c(-0.45, 0.95 + maxchar * 0.3),
    ylim = c(0.3, n + 0.9),
    clip = "off"
  ) +
  labs(title = paste0("Source clusters (n = ", n, ")")) +
  theme_void() +
  theme(
    plot.title = element_text(size = 10, hjust = 0, margin = margin(b = 4)),
    plot.margin = margin(6, 8, 6, 8)
  )

w <- 1.6 + maxchar * 0.085
h <- max(1.6, n * 0.26 + 0.8)
ggsave(args$out_file, p, width = w, height = h, limitsize = FALSE)

cat("__ done! __\n")
