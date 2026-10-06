# Four-population example

A small simulated dataset that runs through the whole workflow in a few minutes on a laptop-sized machine:
48 individuals, 22 human-length chromosomes, 2.2 MB of IBD segments.

```
        O   S1  S2   X        O: outgroup, split 120 generations ago
        |    \  /    |        S1, S2: two sources, split 60 generations ago
        |     S12    |        X: formed 12 generations ago from S1 (0.6) and S2 (0.4), size 1200 since then
        +-----+------+
```

| population | individuals | effective size |
|---|---|---|
| O, S1, S2 | 12 each | 2000 |
| X (admixed) | 12 | 1200 (smaller than the others, so X drifts a little and forms its own cluster) |

The true ancestry of the 12 X individuals, from the simulated genealogies, averages 0.595 S1 and 0.405 S2
(`expected/realized_ancestry.tsv`). `simulation/scenario.yaml` is the scenario file (generated with the simulator's
`ibdsim`, one segment file per chromosome, minimum length 1 cM).

## Run

```
snakemake --cores 8
```

This uses `config/config.yml` as shipped: masking at 2 cM, clustering with a plain cut at height 1.0 and the default
`deep_split` of 1 (see below), and the default mixture settings. It runs:

1. **Clustering.** Four clusters, each exactly one of the true populations: `C0` is O, `C6` is S1, `C4` is S2 and `C7` is
   X (`expected/clusters.tsv`; the individuals of each also carry deeper sub-labels, which are not used here). The
   clustering is deterministic, so you should get the same table.
2. **The mixture panel `four_pop`** (`config/panels/default/mixture_four_pop.tsv`): O, S1 and S2 individuals as
   sources, the 12 X individuals as targets. Sources are pooled by cluster.
3. **The automatic mixture panel `auto`** (see below).

Results appear under `results/panels/cluster_h1.0_dcosine_nscale_mward_D2/`: `clustering/` (the cluster table, heatmap
and hierarchy), `aggregation/` and `mixmodel/four_pop/` and `mixmodel/auto/` (tables, plots and diagnostics).
`aggregation.color_map_embedding` is set to `mds3` because the default `tsne3` needs more clusters than its perplexity
allows.

## What the clustering settings do

The clusters come from a hierarchical tree of the individuals, cut with the adaptive `dynamicTreeCut` method. Two
settings decide how many clusters you get:

- **`base_height` / `gate_height`**: the height at which the tree is cut. Equal values give a plain cut (no gating).
  Lower heights give more, smaller clusters.
- **`deep_split`** (0-4; the `dynamicTreeCut` default is 1, used here): how readily a branch below the cut height is split
  into separate clusters. At 0 only a clearly bimodal branch is split, and higher values also split off cohesive
  sub-branches, so the same height gives more clusters.

In this example the tree has three large merges, at heights 0.40, 1.35 and 3.82, and everything else is below 0.22:
the four populations are four well separated branches, each with some structure of its own (sub-branches of 2-3
individuals). Number of clusters (the four true populations are four clusters):

| `deep_split` | height 0.2 | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|
| 0 | 9 | 4 | 4 | 3 |
| 1 (default, used) | 10 | 6 | **4** | 4 |
| 2 | 10 | 8 | 4 | 4 |
| 3 | 10 | 10 | 8 | 6 |

With the default, any height from 0.7 to 1.5 gives the four populations; the example uses 1.0, the middle. Larger
`deep_split` (3 is used for the large worldwide analyses) splits each population into several sub-clusters here and
needs a height of 2.0 or more to give four. Finer clusters are wanted when they correspond to real sub-populations.

## What to expect for `four_pop`

Mean estimated shares for the X individuals (`expected/mixture_four_pop_summary.tsv`):

| method | S1 | S2 | realized |
|---|---|---|---|
| Bayesian | 0.624 | 0.375 | 0.595 / 0.405 |
| NNLS | 0.612 | 0.387 | 0.595 / 0.405 |

The outgroup gets about 0 (0.001 or less), and individual estimates vary by about 0.025 around these means. The seeds are
fixed in `config/config.yml` (`aggregation.seed`, `mixture.seed`), so a rerun gives the same numbers; small differences
can come from software versions.

Both estimates are slightly above the realized S1 share (Bayesian +0.03, NNLS +0.02). Two things contribute, and for NNLS
(checked by refitting the population-level palettes) removing either one removes the offset:

- **Sharing inside X.** 39% of an X individual's IBD is shared with other X individuals (`self_share` in the output
  tables). No source can explain that part, and the fit spreads it over S1 and S2. The more drifted a target is, the larger
  this share and the offset. The workflow fits all targets together against the full palette; as a diagnostic only,
  refitting without the X donor column gives 0.591 for S1.
- **The palette scale.** Normalized palettes weight a source by its ancestry share times its total IBD per individual,
  and S1 has 9% more IBD per individual than S2 (see [docs/DIAGNOSTICS.md](../docs/DIAGNOSTICS.md#palette-scale)). With
  `mixture.palette_scale: raw` S1 is 0.605 (Bayesian) and 0.593 (NNLS).

## About the `auto` panel

The automatic source picker is a convenience, not a guarantee. On this dataset it picks O, S2 and X as sources, so the S1
individuals become targets and are explained as 99.9% X: X, which carries 60% S1 ancestry, is closer to S1 than anything else it was offered. The picker takes one source per top-level clade and cannot tell that a small, drifted
admixed population should be a target and not a source. For a real analysis, check the picked sources, and prefer a
hand-written mixture file like `four_pop` when you know your sources.
