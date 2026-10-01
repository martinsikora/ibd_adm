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
  library(ggplot2)
  library(stringr)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in",
  action = "store",
  dest = "in_file",
  help = "File with mixmodel results"
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

parser$add_argument("-o", "--out",
  action = "store",
  dest = "out_file",
  help = "output file"
)

parser$add_argument("-m", "--min_p",
  action = "store",
  default = 0.01,
  help = "Minimum ancestry proportion to plot SE  [default %(default)s]"
)

parser$add_argument("--source_grid",
  action = "store_true",
  default = FALSE,
  help = "Plot as grid with sources in separate rows"
)

args <- parser$parse_args()

## Bars are drawn with a hairline stroke: at linewidth 0.25 the outline is wider
## than a bar in a panel this size, so it covers most of the band and the plot
## reads far darker than the palette. Error bars and the 0/1 rules
## keep 0.25.
BAR_LW <- 0.05



## --------------------------------------------------
## read input data

cat("__ reading data __\n")

d_mixmodel <- read_tsv(args$in_file,
  col_types = cols(
    .default = "c",
    p = "d",
    se = "d",
    res_norm = "d"
  )
)

sample_map <- read_tsv(args$sample_file,
  col_types = "cc"
)


color_map <- read_tsv(args$color_file,
  col_types = cols(.default = "c")
)
if (!"fill" %in% colnames(color_map)) {
  color_map$fill <- color_map$color
}


## --------------------------------------------------
## set up helpers and plot data

## plotting helpers

cat("__ preparing plot data __\n")

pal_c <- color_map$color
names(pal_c) <- color_map$pop_id

pal_f <- color_map$fill
names(pal_f) <- color_map$pop_id

## Recipient (cluster_min_dist) samples carry a "_r" suffix on their pop_id (added
## in the mixmodel sample_map). Expand the pop_id levels so each "<cluster>_r"
## facet sorts immediately after its parent "<cluster>", separating cluster_full
## from cluster_min_dist samples of the same cluster. source_pop keeps base levels
## (sources are always cluster_full / donor_recipient).
pop_levels <- as.vector(rbind(color_map$pop_id, paste0(color_map$pop_id, "_r")))

th <- theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.1),
    axis.text.x = element_text(
      angle = 90,
      size = 5,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.x = element_text(
      angle = 90,
      size = 7,
      hjust = 0,
      vjust = 0.5
    ),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0.1, "lines")
  )


## plot data
d1 <- d_mixmodel |>
  mutate(
    sample_id = factor(sample_id, levels = sample_map$sample_id),
    pop_id = factor(pop_id, levels = pop_levels),
    source_pop = factor(source_pop, levels = color_map$pop_id)
  ) |>
  arrange(sample_id, pop_id) |>
  mutate( 
    sample_label = paste(sample_id, label, sep = " | "),
    sample_label = factor(sample_label, levels = unique(sample_label))
  ) |>
  mutate(
    sample_id = droplevels(sample_id),
    pop_id = droplevels(pop_id),
    source_pop = droplevels(source_pop),
    sample_label = droplevels(sample_label)
  )

d_s1 <- d1 |>
  filter(group == "source") |>
  arrange(sample_id, pop_id)

d_target <- d1 |>
  filter(group == "target") |>
  mutate(
    sample_id = droplevels(sample_id),
    pop_id = droplevels(pop_id),
    source_pop = droplevels(source_pop),
    sample_label = droplevels(sample_label)
  )

n_samples <- d1 |>
  pull(sample_label) |>
  unique() |>
  length()
w <- max(1, n_samples %/% 15 + 3)

h <- d1 |>
  pull(sample_label) |>
  str_length() |>
  max(na.rm = TRUE)
h <- ifelse(is.finite(h), h %/% 20, 1)


## --------------------------------------------------
## plot barplot or grid

cat("__ generating plot __\n")

if (args$source_grid) {
  h1 <- d_target |>
    pull(source_pop) |>
    unique() |>
    length() %/% 1.5 + 3

  pdf(args$out_file,
    width = w,
    height = h + h1
  )

  p <- ggplot(d_target, aes(
    x = sample_label,
    y = p
  ))
  print(p +
    geom_col(
      aes(
        color = source_pop,
        fill = source_pop
      ),
      linewidth = BAR_LW
    ) +
    geom_errorbar(
      aes(
        ymax = p + se,
        ymin = p - se,
        group = source_pop
      ),
      linewidth = 0.25,
      width = 0.25,
    ) +
    geom_col(
      aes(
        color = source_pop,
        fill = source_pop
      ),
      linewidth = BAR_LW,
      data = d_s1
    ) +
    geom_text(
      y = 0.5,
      label = "S",
      size = 2,
      data = d_s1,
    ) +
    geom_hline(
      yintercept = c(0, 1),
      linewidth = 0.25
    ) +
    facet_grid(source_pop ~ pop_id,
      space = "free_x",
      scales = "free_x"
    ) +
    scale_fill_manual(
      name = "Source population",
      values = pal_f
    ) +
    scale_color_manual(
      name = "Source population",
      values = pal_c
    ) +
    scale_size_manual(values = c(0.2, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    xlab("") +
    ylab("") +
    guides(fill = guide_legend(nrow = 2)) +
    th +
    theme(strip.text.y = element_text(
      angle = 0,
      size = 7,
      hjust = 0,
      vjust = 0.5
    ), ))

  dev.off()
} else {
  d1_m <- d1 |>
    filter(p >= args$min_p) |>
    select(sample_id, source_pop)

  d2 <- d1 |>
    group_by(sample_id) |>
    arrange(desc(source_pop)) |>
    mutate(
      p_sum = cumsum(p),
      p_low = p_sum - se
    ) |>
    ungroup() |>
    semi_join(d1_m, by = c("sample_id", "source_pop"))

  pdf(args$out_file,
    width = w,
    height = h + 4
  )

  p <- ggplot(d1, aes(
    x = sample_label,
    y = p
  ))
  print(p +
    geom_col(
      aes(
        color = source_pop,
        fill = source_pop
      ),
      linewidth = BAR_LW
    ) +
    geom_errorbar(
      aes(
        ymax = p_sum,
        ymin = p_low,
        group = source_pop
      ),
      linewidth = 0.25,
      width = 0.25,
      data = d2
    ) +
    geom_text(
      y = 0.5,
      label = "S",
      size = 2,
      data = d_s1,
    ) +
    geom_point(
      aes(
        size = res_norm,
        alpha = res_norm
      ),
      fill = "grey",
      color = "grey40",
      y = 1.05,
      shape = 22
    ) +
    geom_hline(
      yintercept = c(0, 1),
      linewidth = 0.25
    ) +
    facet_grid(. ~ pop_id,
      space = "free_x",
      scales = "free_x",
    ) +
    scale_fill_manual(
      name = "Source population",
      values = pal_f
    ) +
    scale_color_manual(
      name = "Source population",
      values = pal_c
    ) +
    scale_size_continuous(range = c(0.1, 2)) +
    scale_alpha_continuous(range = c(0.1, 1)) +
    coord_cartesian(ylim = c(0, 1.07)) +
    xlab("") +
    ylab("") +
    guides(fill = guide_legend(nrow = 1)) +
    th)

  dev.off()
}

cat("__ done! __\n")
