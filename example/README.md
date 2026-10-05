# Four-population example

A small simulated dataset that runs through the whole workflow in a few minutes on a laptop-sized machine:
48 individuals, 22 human-length chromosomes, 2.2 MB of IBD segments.

```
        O   S1  S2   X        O: outgroup, split 120 generations ago
        |    \  /    |        S1, S2: two sources, split 60 generations ago
        |     S12    |        X: formed 12 generations ago from S1 (0.6) and S2 (0.4), size 400 since then
        +-----+------+
```

| population | individuals | effective size |
|---|---|---|
| O, S1, S2 | 12 each | 2000 |
| X (admixed) | 12 | 400 (small, so X drifts and forms its own cluster) |

The true ancestry of the 12 X individuals, from the simulated genealogies, averages 0.604 S1 and 0.396 S2
(`expected/realized_ancestry.tsv`). `simulation/scenario.yaml` is the scenario file (generated with the simulator's
`ibdsim`, one segment file per chromosome, minimum length 1 cM).

## Run

```
snakemake --cores 8
```

This uses `config/config.yml` as shipped: masking at 2 cM, clustering with base height 0.5 gated to 0.2
(`min_sharing` 250000), and the default mixture settings. It runs:

1. **Clustering.** Ten clusters, each from a single true population: three for O, three for S2, two for S1 and two for X
   (`expected/clusters.tsv`). The clustering is deterministic, so you should get the same table.
2. **The mixture panel `four_pop`** (`config/panels/default/mixture_four_pop.tsv`): O, S1 and S2 individuals as
   sources, the 12 X individuals as targets. Sources are pooled by cluster.
3. **The automatic mixture panel `auto`** (see below).

Results appear under `results/panels/cluster_h0.5g0.2_dcosine_nscale_mward_D2/mixmodel/`.

## What to expect for `four_pop`

Mean estimated shares for the X individuals (`expected/mixture_four_pop_summary.tsv`):

| method | S1 | S2 | realized |
|---|---|---|---|
| Bayesian | 0.696 | 0.304 | 0.604 / 0.396 |
| NNLS | 0.617 | 0.383 | 0.604 / 0.396 |

The outgroup gets 0. Individual estimates vary by about 0.03 around these means. The seeds are fixed in
`config/config.yml` (`aggregation.seed`, `mixture.seed`), so a rerun should give the same numbers; small differences
can come from the software versions.

With the default normalized palettes the Bayesian estimate is 0.09 above the realized S1 share and NNLS 0.01. Normalized
palettes weight a source by its ancestry share times its total IBD per individual (see
[docs/DIAGNOSTICS.md](../docs/DIAGNOSTICS.md#palette-scale)). The offset is that known bias. Set `mixture.palette_scale: raw` to compare.

## About the `auto` panel

The automatic source picker is a convenience, not a guarantee. On this dataset it picks the outgroup, one S2 cluster and
one X cluster as sources, so X is not a target in the `auto` panel (the other individuals of X are fitted against the
sources it chose, and X itself gets no estimate). The picker takes one source per top-level clade and cannot tell that a
small, drifted admixed population should be a target. For a real analysis, check the picked sources, and prefer a
hand-written mixture file like `four_pop` when you know your sources.
