# ibd_adm

A Snakemake workflow for **IBD-based ancestry and admixture modelling**.

`ibd_adm` is a *downstream* pipeline: it does **not** phase genotypes or call
identity-by-descent (IBD) itself. It starts from **precomputed pairwise IBD
segments** (one file per chromosome) plus a sample sheet, and turns them into

- an artefact mask over regions of excess IBD coverage,
- masked per-pair total IBD sharing,
- a hierarchical clustering of individuals into populations,
- population × population IBD-sharing profiles and a TVD (total variation
  distance) matrix with an automatically derived colour scheme,
- **admixture / mixture-model estimates** per target individual (NNLS and a
  Bayesian MCMC estimator), with residual diagnostics,
- PCA of the sharing profiles, and
- (optionally) a genome-window scan for population-specific IBD peaks.

IBD calling and phasing (e.g. IBDseq, hap-IBD, refined-IBD) happen **before**
this workflow; `ibd_adm` consumes their segment output.

---

## Pipeline overview

The workflow is assembled in [`workflow/Snakefile`](workflow/Snakefile) from
seven rule modules under [`workflow/rules/`](workflow/rules). Each stage is a
`config.yml` section and can be tuned or (for clustering / mixture / window
peaks) switched off independently.

```
Precomputed IBD segments ({chrom}...ibd.gz)  +  config/individuals.tsv
      │
 1. ibd_mask        get_genomecov → mask_ibd
      │             per-base IBD coverage → excess-coverage mask BED (+ QC PDF)
      │
 2. ibd_tot         get_excl_mask → get_ibd
      │             length/LOD filter, subtract mask, sum total IBD per pair
      │
 3. cluster_ibd     build_feature_matrix → distance_matrix → hclust
      │                 → cut_tree → plot_clusters
      │             sample×sample IBD matrix → distance → hierarchical clustering
      │             → adaptive tree cut (one panel per height) → cluster labels
      │
 4. aggregate_ibd   make_default_agg_panel → aggregate_ibd → tvd
      │                 → default_color_map → tvd_plot
      │             per-population IBD sharing, TVD matrix, colour/shape map
      │
      ├─ 5. mixmodel_ibd   make_mix_sample_map → make_default_mixture_auto
      │                        → run_models (nnls + bayesian)
      │                        → plot_mixmodel / grid / source_legend
      │                        → residual_profiles → residual_diagnostic
      │                    per-target ancestry proportions + SE/ESS/Rhat, barplots
      │
      ├─ 6. pca_ibd        pca (+ pca projected onto source populations)
      │
      └─ 7. ibd_window_peaks  coverage_by_pop → window_average → outliers → plots
                           (disabled by default)
```

Everything from stage 4 onward runs in **two parallel panel families**:

- **default** — the built-in hierarchical clustering, producing one panel per cut
  height (e.g. `cluster_h0.75_...`, `cluster_h1.0_...`); and
- **custom** — user-supplied panels listed under `aggregation.panels`, each with
  its own population definitions and colours (the shipped example is
  `example_panel`).

---

## Repository layout

```
config/
  config.yml                     # all workflow parameters (see docs/CONFIGURATION.md)
  chromosomes.txt                # chromosomes to process, one per line
  n_markers.tsv                  # per-chromosome marker counts (chrom, n)
  individuals.tsv                # sample sheet (sample_id, label, group)   [EXAMPLE]
  metadata/                      # inputs for build_example_panel.R          [EXAMPLE]
    sample_info.tsv
    cluster_info.tsv
  panels/
    example_panel/               # a custom panel                            [EXAMPLE]
      aggregate.tsv              #   sample_id, pop_id, group
      color_map.tsv             #   pop_id, color, fill, shape
      mixture_example.tsv       #   sample_id, group (target|source)
workflow/
  Snakefile
  rules/*.smk                    # the 7 pipeline stages
  scripts/{python,r,awk}/        # step implementations
docs/
  CONFIGURATION.md               # full per-knob reference
results/                         # all outputs (generated; git-ignored)
```

Files marked **[EXAMPLE]** are small synthetic placeholders that document the
required schema. Replace them with your real data before running (see
[Supplying real data](#supplying-real-data)).

---

## Inputs

### IBD segments (`input_data.ibd`)

One file per chromosome; the `{chrom}` wildcard is filled from
`config/chromosomes.txt`. Files are gzip-compressed, whitespace/TSV delimited,
with a one-line header (skipped). The workflow reads these columns:

| column | meaning                    |
|-------:|----------------------------|
| 1      | sample id (first haplotype owner)  |
| 2      | sample id (second haplotype owner) |
| 3      | chromosome                 |
| 4      | segment start (bp)         |
| 5      | segment end (bp)           |
| 6      | LOD / score                |
| 9      | segment length (cM)        |

Columns 7–8 are ignored. Both sample ids must appear in `individuals.tsv`.

### Sample sheet (`input_data.individuals` → `config/individuals.tsv`)

| column      | meaning |
|-------------|---------|
| `sample_id` | must match the ids in the IBD segment files |
| `label`     | free-text population/date label |
| `group`     | `cluster_full` (used to build the clustering tree), `cluster_min_dist` (assigned to a cluster afterwards by k-NN vote), or `exclude` |

### Reference (`ref`)

- `genome` — chromosome-length file (`chrom  length`) used by `bedtools genomecov`.
- `chromosomes` — the chromosome list (`config/chromosomes.txt`).
- `marker_file` — per-chromosome marker counts (`config/n_markers.tsv`), used by
  the mixture model to weight IBD by marker density.
- `fasta` — reserved; read from config but not currently used by any rule.

### Custom-panel files (`config/panels/<panel>/`)

| file                 | columns                          | purpose |
|----------------------|----------------------------------|---------|
| `aggregate.tsv`      | `sample_id, pop_id, group`       | maps samples to populations; `group` ∈ `{donor_recipient, recipient}`. A `pop_id` of `exclude` drops the sample. |
| `color_map.tsv`      | `pop_id, color, fill, shape`     | plotting colours/shapes; its `pop_id` set must match `aggregate.tsv`. Optional (panels without it skip TVD/PCA plots). |
| `mixture_<set>.tsv`  | `sample_id, group`               | defines one mixture experiment named `<set>`; `group` ∈ `{target, source}`. |

The default clustering panels generate their `aggregate.tsv` / `color_map.tsv`
automatically. For the default panels a mixture set named `auto` is always
available (sources auto-selected from the TVD tree); add
`config/panels/default/mixture_<set>.tsv` to define named sets by hand.

---

## Outputs

All outputs are written under `results/` (git-ignored):

```
results/
  masking/{tables,plots}/            # coverage tables, mask BEDs, QC PDFs
  ibd_tot/tables/                    # masked per-pair total IBD
  cluster_cache/<tag>/               # shared clustering intermediates (.rds)
  panels/<panel>/
    clustering/                      # cluster assignments + dendrogram/heatmap
    aggregation/{tables,plots,panels}/   # ibd_pop.tsv.gz, TVD matrix, colour map
    mixmodel/<set>/{tables,plots,diagnostics}/   # model estimates, barplots, flags
    pca/<set>/{tables,plots}/        # PCA tables + plots
    ibd_window_peaks/{tables,plots}/ # (when enabled)
```

`<panel>` is either a generated clustering panel (`cluster_h<height>_<tag>`) or a
custom panel name. `<tag>` encodes the distance/transform/agglomeration choices
(`d<dist>_n<transform>_m<clust>`), so changing those settings caches to a fresh
directory rather than overwriting.

---

## Requirements

- **Snakemake** (workflow engine).
- **Command-line tools** on `PATH`: `bedtools`, GNU `datamash`, `gawk`, `gzip`,
  `sort`.
- **Python 3** with `pandas` and `numpy`.
- **R** with:
  `argparse`, `readr`, `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`, `rlang`,
  `data.table`, `ggplot2`, `scales`, `scico`, `grid`, `grDevices`,
  `future`, `furrr`, `parallelDist`, `fastcluster`, `dynamicTreeCut`,
  `lsei`, `Rtsne`, `ape`, `phytools`, `data.tree`, `igraph`, `tidygraph`,
  `ggraph`, `ggtree`, `heatmap3`.

No conda environment is committed; install the above with your preferred manager
(conda/mamba, `renv`, system packages, …).

---

## Running

From the repository root (where `config/` lives):

```bash
# 1. Dry run: build the DAG and validate config without executing anything
snakemake -n

# 2. Run locally with N cores
snakemake --cores 16

# 3. Build a specific target
snakemake --cores 8 results/ibd_tot/tables/1.example_dataset.ibd_tot.tsv.gz
```

Several rules declare `mem_mb` and `runtime` resources (e.g. the clustering and
mixture-model steps need a lot of memory and CPU). These imply a cluster
executor. Supply a Snakemake **profile / SLURM executor** at invocation
(e.g. `snakemake --workflow-profile <profile>` or `--executor slurm`). No profile
is committed to this repo; provide one suited to your scheduler.

Two `config.yml` knobs throttle IO-heavy fan-out via Snakemake global resources:
`aggregation.max_concurrent_ibd_jobs` and
`ibd_window_peaks.max_concurrent_coverage_jobs`.

---

## Supplying real data

1. Point `input_data.ibd` at your per-chromosome IBD segment files (keep the
   `{chrom}` wildcard) and `input_data.individuals` at your sample sheet.
2. Set `ref.genome` (and `ref.chromosomes` / `ref.marker_file`) to your reference.
3. Set `prefix` to your dataset name (it appears in every output filename).
4. Optionally set `tmpdir` to a fast scratch location.
5. Replace or remove the `example_panel` entry under `aggregation.panels`, and
   provide the corresponding `config/panels/<name>/` files (or rely solely on the
   default clustering panels).

---

## Building a custom panel

A custom panel is the three files under `config/panels/<name>/` described
above, which you can write by hand. Alternatively,
[`workflow/scripts/r/build_example_panel.R`](workflow/scripts/r/build_example_panel.R)
builds `aggregate.tsv` and `color_map.tsv` together, so their `pop_id` sets stay
identical, from two metadata tables:

- `sample_info.tsv` — `sample_id, cluster_label, cluster_alias, cluster_assignment`
- `cluster_info.tsv` — `cluster_label, cluster_alias, color, fill, shape`

Run it against the shipped example metadata to regenerate `example_panel`:

```bash
Rscript workflow/scripts/r/build_example_panel.R \
  --sample_info  config/metadata/sample_info.tsv \
  --cluster_info config/metadata/cluster_info.tsv \
  --out_dir      config/panels/example_panel
```

---

## Mixture-model notes

For each **target** sample, its IBD-sharing vector is modelled as a non-negative
mixture over **source**-population profiles. Two estimators are available
(`mixture.method`):

- **`nnls`** — sum-to-one non-negative least squares (`lsei::pnnls`) with
  per-chromosome block-jackknife standard errors.
- **`bayesian`** — a SOURCEFIND-style MCMC with a Dirichlet proposal, adaptive
  proposal scaling, and an active-source search (spike-and-slab over the source
  palette); reports acceptance rate, ESS, and R-hat. Judge convergence on `rhat_median`; `rhat_max` becomes large for near-zero sources.

Source populations can be listed explicitly (a `mixture_<set>.tsv` with
`group == source`) or auto-selected from the TVD / neighbour-joining tree
(`mixture_auto`), controlled by the `mixture.auto_source_*` knobs. A post-hoc
**residual diagnostic** (when both estimators run) flags source populations that
behave as poor proxies and names the unused population they are standing in for.

Each fit reports `res_norm_ex_self` (the residual excluding the target's own cluster), which is the statistic to compare across targets. Sources are also screened for an R scale offset (`source_R_flags.tsv`), and NNLS fits can be evaluated by chromosome hold-out CV (`mixture.cv`). The output columns are described in [`docs/CONFIGURATION.md`](docs/CONFIGURATION.md#mixture-output-tables).

See [`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) for every knob.

---

## License

GPL-2.0-or-later. See [`LICENSE`](LICENSE). Source files carry
`Copyright Martin Sikora <martin.sikora@sund.ku.dk>`.
