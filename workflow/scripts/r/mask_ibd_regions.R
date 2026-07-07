#!/usr/bin/env Rscript
## --------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(argparse)
  library(scales)
})

## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in_file",
  action = "store",
  dest = "in_file",
  help = "Input file genomecov"
)

parser$add_argument("-b", "--bed_file",
  action = "store",
  dest = "bed_file",
  help = "Output BED filename"
)

parser$add_argument("-p", "--plot_file",
  action = "store",
  dest = "plot_file",
  help = "Output plot filename"
)

parser$add_argument("--trim",
  action = "store",
  dest = "trim",
  type = "double",
  default = 0.05,
  help = "Percentile for trimmed mean for excess IBD region removal [default %(default)s]"
)

parser$add_argument("--sd",
  action = "store",
  dest = "sd",
  type = "double",
  default = 5,
  help = "Standard deviations from trimmed mean for excess IBD region removal [default %(default)s]"
)

args <- parser$parse_args()


## ---------------------------------------------------------
## helpers

th <- theme_bw() +
  theme(
    panel.grid.major = element_line(
      linetype = "dotted",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
  )


## ---------------------------------------------------------
## read and process genomecov data

d_cov <- read_tsv(
  file = args$in_file,
  col_names = c("chromosome", "pos_start", "pos_end", "n_seg")
)

d_cov_trim <- d_cov |>
  arrange(n_seg) |>
  slice(max(1L, round(n() * args$trim)):round(n() * (1 - args$trim))) |>
  summarise(
    ibd_mean_tr = mean(n_seg),
    ibd_sd_tr = if (n() >= 2) sd(n_seg) else 0
  )


## ---------------------------------------------------------
## plot and write bed file with annotations

ibd_cutoff <- d_cov_trim$ibd_mean_tr + args$sd * d_cov_trim$ibd_sd_tr
o <- d_cov |>
  mutate(flag = case_when(
    n_seg >= ibd_cutoff ~ paste("ibd_excess_", args$sd, "_sd", sep = ""),
    TRUE ~ "0"
  ))

w <- max(o$pos_end) %/% 2e7 + 5

pdf(
  file = args$plot_file,
  width = w,
  height = 4
)
p <- ggplot(o, aes(
  x = pos_start / 1e6,
  y = n_seg
))
p +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_hline(
    yintercept = ibd_cutoff,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_step(
    aes(
      color = flag,
      group = NA
    ),
    linewidth = 0.5
  ) +
  scale_colour_manual(values = setNames(
    c("black", "tomato"),
    c("0", paste0("ibd_excess_", args$sd, "_sd"))
  )) +
  scale_y_continuous(label = comma) +
  xlab("Position (Mb)") +
  ylab("IBD pairs") +
  th
dev.off()

o1 <- o |>
  mutate(idx = cumsum(flag != lag(flag, default = "0"))) |>
  group_by(idx, chromosome) |>
  summarise(
    pos_start = min(pos_start),
    pos_end = max(pos_end),
    flag = flag[1],
    .groups = "drop"
  ) |>
  select(-idx)

write_tsv(o1,
  file = args$bed_file,
  col_names = FALSE
)
