#!/usr/bin/env python3
"""Curate the outputs of a finished example run and render the walkthrough pages.

  example/make_walkthrough.py curate <run dir> [--raw <run dir>] [--noS2 <run dir>] [--ds3 <run dir>]
  example/make_walkthrough.py pages
  example/make_walkthrough.py check <run dir>

curate  copies the tables and figures the walkthrough shows from a finished run of the example (`snakemake --cores 8` in
        <run dir>) into example/results/, renders the PDF plots to PNG (needs Ghostscript `gs`) and draws the dendrogram with the
        cut (needs Rscript). --raw, --noS2 and --ds3 are finished runs of the variants in example/variants/.
pages   writes docs/walkthrough/*.md from example/results/, so every number and excerpt in the text comes from the shipped
        outputs.
check   compares a finished run with example/results/ (clusters identical, mixture tables equal within 1e-6), to catch drift.
"""
import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = f"{REPO}/example/results"
PAGES = f"{REPO}/docs/walkthrough"
PFX = "example_dataset"


# ------------------------------------------------------------------ helpers
def panel_dir(run):
    d = sorted(glob.glob(f"{run}/results/panels/cluster_h*"))
    if not d:
        sys.exit(f"no results/panels/cluster_h* in {run}")
    return d[0]


def cut_height(pdir):
    return float(re.search(r"cluster_h([\d.]+?)_d", os.path.basename(pdir)).group(1))


def put(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy(src, dst)


def png(pdf, out, dpi=130):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", f"-r{dpi}", "-dFirstPage=1", "-dLastPage=1", "-o", out, pdf], check=True)


def md(df, fmt=None, index=False):
    """DataFrame -> GitHub markdown table."""
    fmt = fmt or {}
    cols = ([df.index.name or ""] if index else []) + list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for i in range(len(df)):
        cells = [str(df.index[i])] if index else []
        for c in df.columns:
            v = df[c].iloc[i]
            if c in fmt and pd.notna(v):
                v = fmt[c].format(v)
            cells.append("" if pd.isna(v) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def aligned(text, right=()):
    """Tab-separated text -> columns padded to a common width (the files themselves stay tab-separated).
    Columns whose name is in `right` are right-aligned."""
    rows = [l.split("\t") for l in text.rstrip("\n").split("\n")]
    w = [max(len(r[i]) for r in rows if i < len(r)) for i in range(max(len(r) for r in rows))]
    hdr = rows[0]
    out = []
    for r in rows:
        cells = [c.rjust(w[i]) if i < len(hdr) and hdr[i] in right else c.ljust(w[i]) for i, c in enumerate(r)]
        out.append("  ".join(cells).rstrip())
    return "\n".join(out)


def head_text(path, n):
    with open(path) as f:
        return "".join([next(f) for _ in range(n)])


# ------------------------------------------------------------------ curate
def curate_mix(P, setname, dst):
    m = f"{P}/mixmodel/{setname}"
    for meth in ("bayesian", "nnls"):
        put(f"{m}/tables/{PFX}.mixmodel_{meth}.tsv", f"{dst}/mixmodel_{meth}.tsv")
        png(f"{m}/plots/{PFX}.mixmodel_{meth}.pdf", f"{dst}/mixmodel_{meth}.png")
    for f in ("source_R_flags", "target_R_flags", "residual_diagnostic.cluster_residuals", "residual_diagnostic.source_flags",
              "residual_diagnostic.source_sink_by_stratum"):
        p = f"{m}/diagnostics/{PFX}.{f}.tsv"
        if os.path.exists(p):
            put(p, f"{dst}/diagnostics/{f.replace('residual_diagnostic.', '')}.tsv")


def curate(a):
    run = a.run
    P = panel_dir(run)
    h = cut_height(P)
    shutil.rmtree(RES, ignore_errors=True)
    # clustering
    c = pd.read_csv(f"{P}/clustering/default.{PFX}.clusters.tsv", sep="\t")
    c = c[c.cut_height.astype(float).round(4) == round(h, 4)].rename(columns={"label": "true_population"})
    os.makedirs(f"{RES}/clustering", exist_ok=True)
    c[["sample_id", "cluster_label", "true_population", "cluster_id", "cluster_depth"]].to_csv(f"{RES}/clustering/clusters.tsv", sep="\t", index=False)
    png(f"{P}/clustering/default.{PFX}.clusters_hierarchy.pdf", f"{RES}/clustering/hierarchy.png")
    png(f"{P}/clustering/default.{PFX}.clusters_heatmap.pdf", f"{RES}/clustering/heatmap.png")
    rds = glob.glob(f"{run}/results/cluster_cache/*/{PFX}.res_hc.rds")[0]
    r = f"""hc <- readRDS('{rds}'); h <- sort(hc$height, decreasing = TRUE)[1:8]
write.table(data.frame(rank = 1:8, height = round(h, 3)), '{RES}/clustering/tree_heights.tsv', sep = '\\t', quote = FALSE, row.names = FALSE)
cl <- read.delim('{RES}/clustering/clusters.tsv'); lab <- hc$labels[hc$order]
cc <- cl$cluster_label[match(lab, cl$sample_id)]; pp <- cl$true_population[match(lab, cl$sample_id)]
pal <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666')
colc <- setNames(pal[seq_along(unique(cc))], unique(cc)); colp <- setNames(pal[seq_along(unique(pp))], unique(pp))
png('{RES}/clustering/dendrogram_cut.png', width = 1500, height = 750, res = 130)
par(mar = c(5, 4, 1, 1)); plot(hc, labels = FALSE, hang = -1, main = '', xlab = '', sub = '', ylab = 'merge height')
abline(h = {h}, lty = 2, col = 'red'); par(xpd = NA); m <- max(hc$height); n <- length(lab); x <- seq_len(n)
rect(x - 0.5, -0.10 * m, x + 0.5, -0.05 * m, col = colc[cc], border = NA)
rect(x - 0.5, -0.17 * m, x + 0.5, -0.12 * m, col = colp[pp], border = NA)
for (k in unique(cc)) text(mean(x[cc == k]), -0.215 * m, k, cex = 0.9)
text(0.2, -0.075 * m, 'cluster', adj = c(1, 0.5), cex = 0.8); text(0.2, -0.145 * m, 'population', adj = c(1, 0.5), cex = 0.8)
for (k in unique(pp)) text(mean(x[pp == k]), -0.28 * m, k, cex = 0.9)
dev.off()"""
    subprocess.run(["Rscript", "-e", r], check=True)
    # palettes and TVD
    d = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in glob.glob(f"{P}/aggregation/tables/*.{PFX}.ibd_pop.tsv.gz")])
    d["ibd"] = d.ibd.astype(float)
    pal = d.groupby(["sample1", "pop_id2"]).ibd.sum().unstack(fill_value=0).round(1)
    pal.index.name = "sample_id"
    os.makedirs(f"{RES}/palettes", exist_ok=True)
    pal.to_csv(f"{RES}/palettes/palettes_cM.tsv", sep="\t")
    put(f"{P}/aggregation/tables/{PFX}.ibd_pop_tvd.tsv", f"{RES}/palettes/tvd.tsv")
    png(f"{P}/aggregation/plots/{PFX}.ibd_pop_tvd.pdf", f"{RES}/palettes/tvd_tree.png", 110)
    # mixtures, PCA
    curate_mix(P, "four_pop", f"{RES}/mixture/four_pop")
    curate_mix(P, "auto", f"{RES}/mixture/auto")
    put(f"{P}/mixmodel/panels/mixture_auto.tsv", f"{RES}/mixture/auto/mixture_auto.tsv")
    png(f"{P}/pca/default/plots/{PFX}.pca.pdf", f"{RES}/pca/pca.png", 110)
    png(f"{P}/pca/four_pop/plots/{PFX}.pca_proj.pdf", f"{RES}/pca/pca_proj_four_pop.png", 110)
    # variants
    if a.raw:
        curate_mix(panel_dir(a.raw), "four_pop", f"{RES}/variants/raw")
    if a.noS2:
        curate_mix(panel_dir(a.noS2), "noS2", f"{RES}/variants/noS2")
    if a.ds3:
        P3 = panel_dir(a.ds3)
        c3 = pd.read_csv(f"{P3}/clustering/default.{PFX}.clusters.tsv", sep="\t")
        c3 = c3[c3.cut_height.astype(float).round(4) == round(cut_height(P3), 4)].rename(columns={"label": "true_population"})
        os.makedirs(f"{RES}/variants/deep_split3", exist_ok=True)
        c3[["sample_id", "cluster_label", "true_population"]].to_csv(f"{RES}/variants/deep_split3/clusters.tsv", sep="\t", index=False)
    with open(f"{RES}/README.md", "w") as f:
        f.write("Curated outputs of one run of the example (`snakemake --cores 8`), rendered by `example/make_walkthrough.py`. "
                "They are shown and explained in [docs/walkthrough](../../docs/walkthrough/README.md).\n")
    print("curated into", RES)


# ------------------------------------------------------------------ facts used by the pages
class Facts:
    def __init__(self):
        r = lambda p, **k: pd.read_csv(f"{RES}/{p}", sep="\t", **k)
        self.cl = r("clustering/clusters.tsv")
        self.pop_of = self.cl.groupby("cluster_label").true_population.agg(lambda s: s.iloc[0]).to_dict()
        self.cl_of = {v: k for k, v in self.pop_of.items()}
        self.heights = r("clustering/tree_heights.tsv")
        self.pal = r("palettes/palettes_cM.tsv", index_col=0)
        t = r("palettes/tvd.tsv")
        self.tvd = t.pivot(index="pop_id1", columns="pop_id2", values="tvd")
        self.mix = {m: {meth: r(f"mixture/{m}/mixmodel_{meth}.tsv") for meth in ("bayesian", "nnls")} for m in ("four_pop", "auto")}
        for v in ("raw", "noS2"):
            if os.path.isdir(f"{RES}/variants/{v}"):
                self.mix[v] = {meth: r(f"variants/{v}/mixmodel_{meth}.tsv") for meth in ("bayesian", "nnls")}
        self.exp = pd.read_csv(f"{REPO}/example/expected/mixture_four_pop_summary.tsv", sep="\t").set_index(["method", "source"])

    def mean_p(self, name, meth, targets="X"):
        t = self.mix[name][meth]
        t = t[t.label == targets] if "label" in t else t
        g = t.groupby(["sample_id", "source_pop"]).p.sum().unstack().fillna(0)
        g = g.div(g.sum(axis=1), axis=0)
        return g.mean(), g.std()


def pct(x):
    return f"{x:.3f}"


# ------------------------------------------------------------------ pages
FOOT = "\n---\n[Walkthrough home](README.md) · generated by `example/make_walkthrough.py` from `example/results/`\n"


def write(name, text):
    os.makedirs(PAGES, exist_ok=True)
    with open(f"{PAGES}/{name}", "w") as f:
        f.write(text.strip() + "\n" + FOOT)


def pages(a):
    F = Facts()
    S1, S2, O, X = (F.cl_of[p] for p in ("S1", "S2", "O", "X"))
    R = "../../example/results"
    h = F.heights.height.tolist()
    cross = pd.crosstab(F.cl.cluster_label, F.cl.true_population)[["O", "S1", "S2", "X"]]
    cross.index.name = "cluster"
    bay, bsd = F.mean_p("four_pop", "bayesian")
    nn, nsd = F.mean_p("four_pop", "nnls")
    real = pd.read_csv(f"{REPO}/example/expected/realized_ancestry.tsv", sep="\t")
    real_s1 = real[(real.population == "X") & (real.source == "S1")].share.mean()
    xt = F.mix["four_pop"]["bayesian"]
    xt = xt[xt.label == "X"]
    self_share = xt.groupby("sample_id").self_share.first().mean()

    # ---------------- home
    write("README.md", f"""
# Walkthrough: a four-population example

This walkthrough follows the example in [`example/`](../../example) through the workflow, stage by stage. Every table and
figure is a real output of that example, shipped in [`example/results/`](../../example/results); the text is generated from
those files.

## The data

48 individuals on 22 chromosomes, simulated with a known history:

```
        O   S1  S2   X        O: outgroup, split 120 generations ago
        |    \\  /    |        S1, S2: two sources, split 60 generations ago
        |     S12    |        X: formed 12 generations ago from S1 (0.6) and S2 (0.4)
        +-----+------+
```

12 individuals from each of O, S1, S2 and X, all with effective size 2000. The true ancestry of X, from the simulated genealogies,
is {real_s1:.3f} S1 and {1 - real_s1:.3f} S2 on average ([`example/expected/realized_ancestry.tsv`](../../example/expected/realized_ancestry.tsv)),
so every estimate below can be compared with the truth.

## Run it

```bash
snakemake --cores 8     # a few minutes; results appear under results/
```

## Pages

| page | stage | what you see |
|---|---|---|
| [1. Inputs](1-inputs.md) | | the files the workflow reads |
| [2. Clustering](2-clustering.md) | `cluster_ibd` | the tree, the cut, four clusters |
| [3. Palettes and TVD](3-palettes-and-tvd.md) | `aggregate_ibd` | what an IBD palette is, distances between populations |
| [4. Mixture modelling](4-mixture.md) | `mixmodel_ibd` | ancestry proportions of X, against the truth |
| [5. Reading the diagnostics](5-diagnostics.md) | `mixmodel_ibd` | residuals, source R flags, `self_share` |
| [6. Source selection](6-source-selection.md) | `mixmodel_ibd` | the automatic panel and how to check its picks |
| [7. Variants to try](7-variants.md) | | raw palettes; a missing source |

For every knob, see [`docs/CONFIGURATION.md`](../CONFIGURATION.md); for how to read the outputs in general,
[`docs/DIAGNOSTICS.md`](../DIAGNOSTICS.md).
""")

    # ---------------- inputs
    ind_head = aligned(head_text(f"{REPO}/config/individuals.tsv", 4))
    seg_head = aligned(subprocess.run(f"zcat {REPO}/example/ibd_segments/22.example_dataset.ibdseq.ibd.gz | head -4", shell=True, capture_output=True, text=True).stdout,
                       right=("chrom", "start", "end", "lod", "hap1", "hap2", "cm"))
    write("1-inputs.md", f"""
# 1. Inputs

## IBD segments

One file per chromosome (`input_data.ibd`), here `example/ibd_segments/{{chrom}}.example_dataset.ibdseq.ibd.gz`. The workflow
reads the sample ids (columns 1-2), chromosome (3), start and end in bp (4-5), LOD (6) and length in cM (9) (the files are tab-separated; the excerpts on this page are padded into aligned columns for reading):

```
{seg_head.rstrip()}
```

## Sample sheet

`config/individuals.tsv` (`input_data.individuals`) lists the individuals and their group:

```
{ind_head.rstrip()}
```

`label` is the population name used for plots and for the true population here. `group` is `cluster_full` (used to build the
clusters), `cluster_min_dist` (assigned to a cluster afterwards by k-NN) or `exclude`.

## Other files

| file | key | what it holds |
|---|---|---|
| `config/chromosomes.txt` | `ref.chromosomes` | chromosomes to process, one per line; fills the `{{chrom}}` wildcard |
| `config/genome.txt` | `ref.genome` | chromosome lengths in bp, used for the coverage mask |
| `config/n_markers.tsv` | `ref.marker_file` | marker counts per chromosome; block sizes for the jackknife SEs |
| `config/panels/default/mixture_four_pop.tsv` | | which individuals are sources and which are targets (page 4) |

See [`docs/CONFIGURATION.md`](../CONFIGURATION.md) for `ref`, `input_data` and `masking`.
""")

    # ---------------- clustering
    top = ", ".join(f"{x:.2f}" for x in h[:3])
    T_h = md(F.heights.head(5).rename(columns={"rank": "merge rank"}))
    T_cross = md(cross, index=True)
    TRY = "**Try this:** `example/variants/deep_split3.yml` cuts the same tree at the same height with `deep_split` 3 (instructions in the file)."
    if os.path.exists(f"{RES}/variants/deep_split3/clusters.tsv"):
        c3 = pd.read_csv(f"{RES}/variants/deep_split3/clusters.tsv", sep="\t")
        n3 = c3.groupby("true_population").cluster_label.nunique()
        TRY += (f" It gives {c3.cluster_label.nunique()} clusters instead of 4 (O {n3['O']}, S1 {n3['S1']}, S2 {n3['S2']}, X {n3['X']}): every cluster is still "
                "one population, but the populations are split into sub-clusters ([`example/results/variants/deep_split3/clusters.tsv`](../../example/results/variants/deep_split3/clusters.tsv)).")
    write("2-clustering.md", f"""
# 2. Clustering

**In:** the per-pair IBD totals (after masking regions of excess IBD coverage). **Out:** `clustering/default.example_dataset.clusters.tsv`
and two plots. The individuals are clustered hierarchically (Ward's method on a cosine distance of their IBD vectors) and the
tree is cut with `dynamicTreeCut` at a fixed height.

![dendrogram](../../example/results/clustering/dendrogram_cut.png)

The tree has three large merges, at heights {top}; everything else is below {math.ceil(h[3] * 100) / 100:.2f}. The dashed line is the cut height
(`clustering.base_height` = `gate_height` = 1.0). The coloured strips under the leaves are the clusters the workflow
assigns and the true populations. A plain cut of the tree at 1.0 would give three clusters, because S1 and X join at {h[2]:.2f}.
The adaptive cut (`dynamicTreeCut`) works on the shape of the branches below the height and keeps those two cohesive groups
apart; how readily it splits is `deep_split`.

{T_h}

Compared with the true populations (`example/results/clustering/clusters.tsv`), each cluster is exactly one population:

{T_cross}

![heatmap](../../example/results/clustering/heatmap.png)

## Knobs

| key | here | effect |
|---|---|---|
| `clustering.base_height`, `gate_height` | 1.0, 1.0 | cut height; equal = plain cut, lower `gate_height` = sharing-gated finer cut |
| `clustering.deep_split` | 1 | how readily a branch under the cut is split further (0-4) |
| `clustering.cl_size` | 3 | smallest cluster the cut will emit |
| `clustering.knn` | 7 | k for assigning `cluster_min_dist` individuals |

With `deep_split` 1 these four clusters appear over a range of cut heights around 1.0; with `deep_split` 3 the same tree is cut
into more clusters at that height and needs a higher cut. The ranges are in
[`example/README.md`](../../example/README.md#what-the-clustering-settings-do).

{TRY}
""")

    # ---------------- palettes and TVD
    pal = F.pal.copy()
    pal.insert(0, "population", [F.cl.set_index("sample_id").true_population[i] for i in pal.index])
    ex = pal.loc[["O_001", "S1_001", "S2_001", "X_001"]]
    sh = ex.drop(columns="population")
    shares = sh.div(sh.sum(axis=1), axis=0).round(3)
    shares.insert(0, "population", ex.population)
    tv = F.tvd.copy().round(3)
    tv.index.name = "cluster"
    pn = lambda c: F.pop_of[c]
    far = F.tvd.drop(index=O, columns=O).max().max()
    xs = shares.loc["X_001"]
    T_ex = md(ex.reset_index().rename(columns={"index": "sample_id"}), fmt={c: "{:.0f}" for c in sh.columns})
    T_sh = md(shares.reset_index().rename(columns={"index": "sample_id"}), fmt={c: "{:.3f}" for c in sh.columns})
    T_tv = md(tv, index=True)
    write("3-palettes-and-tvd.md", f"""
# 3. Palettes and TVD

**In:** the clusters. **Out:** `aggregation/tables/*.ibd_pop.tsv.gz` (one per chromosome), `ibd_pop_tvd.tsv` and plots.

## The palette

Every individual has a **palette**: the total IBD (cM) it shares with each donor population (here, each of the four clusters).
The mixture model reproduces a target's palette as a weighted sum of source palettes. Palettes of one individual from each
population, summed over chromosomes (`example/results/palettes/palettes_cM.tsv`):

{T_ex}

As a share of the individual's total:

{T_sh}

Individuals share most IBD with their own cluster. X_001 has {xs[S1]:.2f} of its palette with S1 ({S1}), {xs[S2]:.2f} with S2
({S2}) and {xs[X]:.2f} with other X individuals ({X}). Its S1:S2 ratio is the signal the mixture model uses; the X-to-X sharing is
the part no source can explain (`self_share`, page 5).

To keep the sums over the same number of donors, `aggregate_ibd` sums each individual's within-cluster entry over the n-1 other
members and its between-cluster entries over a random n-1 of the n donors. `aggregation.seed` fixes that draw.

## TVD

The total variation distance between the palettes of the populations (`example/results/palettes/tvd.tsv`) is the distance
used for the population tree and for choosing sources automatically:

{T_tv}

![tvd tree](../../example/results/palettes/tvd_tree.png)

The outgroup ({O}) is far from everything ({F.tvd.loc[O].drop(O).min():.2f} or more). Among the others, X ({X}) is closest to S1 ({S1}, {F.tvd.loc[X, S1]:.2f}), then to S2 ({S2}, {F.tvd.loc[X, S2]:.2f}), and S1 and S2 are {F.tvd.loc[S1, S2]:.2f} apart: X lies between its two sources, nearer the one it takes more ancestry from.

## Knobs

`aggregation.ibd_params` (segment length and LOD filters for the palette), `aggregation.seed`, the colour-map keys
(`aggregation.color_*`; this example uses `color_map_embedding: mds3` because `tsne3` needs more clusters). See
[`docs/CONFIGURATION.md`](../CONFIGURATION.md#aggregation).
""")

    # ---------------- mixture
    mix_head = aligned(head_text(f"{REPO}/config/panels/default/mixture_four_pop.tsv", 3))
    x1 = F.mix["four_pop"]["bayesian"]
    x1 = x1[x1.sample_id == "X_001"].set_index("source_pop")
    n1 = F.mix["four_pop"]["nnls"]
    n1 = n1[n1.sample_id == "X_001"].set_index("source_pop")
    one = pd.DataFrame({"Bayesian p": x1.p, "Bayesian se": x1.se, "NNLS p": n1.p, "NNLS se": n1.se}).reset_index()
    res = pd.DataFrame({"S1": [real_s1, bay[S1], nn[S1]], "S2": [1 - real_s1, bay[S2], nn[S2]], "O": [0.0, bay[O], nn[O]]},
                       index=["realized", "Bayesian", "NNLS"]).round(3)
    res.index.name = "X individuals, mean"
    tot = F.pal.sum(axis=1)
    pop_s = F.cl.set_index("sample_id").true_population
    s1_more = 100 * (tot[pop_s[tot.index] == "S1"].mean() / tot[pop_s[tot.index] == "S2"].mean() - 1)
    T_res = md(res, fmt={"S1": "{:.3f}", "S2": "{:.3f}", "O": "{:.3f}"}, index=True)
    T_one = md(one, fmt={c: "{:.3f}" for c in one.columns if c != "source_pop"})
    xr = x1.iloc[0]
    write("4-mixture.md", f"""
# 4. Mixture modelling

**In:** the palettes and a mixture file, `config/panels/default/mixture_four_pop.tsv`, which says which individuals are sources and
which are targets:

```
{mix_head}
...
```

O, S1 and S2 individuals are the sources (individuals of one cluster are pooled into one source) and the 12 X individuals are the
targets. **Out:** `mixmodel/four_pop/tables/example_dataset.mixmodel_{{nnls,bayesian}}.tsv`, one row per target and source, and plots.

Two estimators run on the same palettes:
- `nnls`: sum-to-one non-negative least squares, with a jackknife standard error over chromosomes;
- `bayesian`: a MCMC with a Dirichlet prior and a search over which sources are active.

## Result

Mean estimated ancestry of the 12 X individuals ({S1} = S1, {S2} = S2, {O} = O), against the truth from the simulation:

{T_res}

Individual estimates vary by about {bsd[S1]:.3f} around the Bayesian mean.

![bayesian](../../example/results/mixture/four_pop/mixmodel_bayesian.png)

The panels are the O, S1, S2 and X individuals and each bar is one individual, coloured by the sources it is assigned to. The O, S1
and S2 individuals are the sources, so each is shown as its own source; the X individuals are the targets and show the mixture.

## One individual

The rows of X_001 in the two tables (`mixmodel_bayesian.tsv`, `mixmodel_nnls.tsv`):

{T_one}

`p` is the weight of the source and `se` its standard error (NNLS: jackknife over chromosomes; Bayesian: posterior). Weights sum
to 1 for a target. The same rows carry the residual of the fit ({xr.res_norm_rmse:.3f} on the RMSE scale for both estimators) and
`self_share` ({xr.self_share:.2f}); they are explained on [page 5](5-diagnostics.md).

## Why the estimates are above the truth

Both estimators give S1 more than {real_s1:.3f}. The main reason is that {self_share:.0%} of an X individual's palette is sharing with
other X individuals, which no source explains; the fit spreads it over S1 and S2. A smaller reason is the palette scale: with the
default `palette_scale: normalized` a source is weighted by its ancestry share times its total IBD per individual, so a source
with more IBD per individual (S1 here, {s1_more:.0f}% more than S2) is slightly over-credited. [Page 7](7-variants.md) refits with raw palettes.

The own-cluster effect is exaggerated in this example: X is a small, recently founded and strongly drifted population, and the panel has only four populations, so a large share of its sharing is within X. In larger datasets with many populations most targets share a much smaller part of their IBD with their own cluster, and the bias from this is usually small. Strongly endogamous cohorts are the exception and can still have a large `self_share`.

## Knobs

`mixture.method` (`nnls`, `bayesian`, `both`), `mixture.seed`, `mixture.palette_scale`, `mixture.mean_active_sources`,
`mixture.mcmc_iter` and the other MCMC keys, `mixture.two_stage_se`. See [`docs/CONFIGURATION.md`](../CONFIGURATION.md#mixture).
""")

    # ---------------- diagnostics
    srf = pd.read_csv(f"{RES}/mixture/four_pop/diagnostics/source_R_flags.tsv", sep="\t", comment="#")
    trf = pd.read_csv(f"{RES}/mixture/four_pop/diagnostics/target_R_flags.tsv", sep="\t")
    cres = pd.read_csv(f"{RES}/mixture/four_pop/diagnostics/cluster_residuals.tsv", sep="\t")
    sf = pd.read_csv(f"{RES}/mixture/four_pop/diagnostics/source_flags.tsv", sep="\t")
    srf = srf.rename(columns={"component": "source"})
    o = srf[srf.source == O].iloc[0]
    rep = {k: trf.risk.value_counts().get(k, 0) for k in ("none", "LOW", "MODERATE", "HIGH")}
    T_srf = md(srf[["source", "n", "R_excl_self", "fold_vs_median", "direction", "flag"]], fmt={"R_excl_self": "{:.0f}", "fold_vs_median": "{:.2f}"})
    trs = trf.sort_values("q_flagged", ascending=False)
    T_trf = md(trs.head(4)[["sample_id", "p_warn", "q_warn", "risk", "res_ratio", "poor_fit"]], fmt={"p_warn": "{:.4f}", "q_warn": "{:.4f}"})
    tt = trs.iloc[0]
    T_cres = md(cres[cres.pop_id == X][["pop_id", "res_norm", "top_source", "top_p", "miss_pop", "miss_score"]], fmt={"res_norm": "{:.4f}", "top_p": "{:.3f}", "miss_score": "{:.2f}"})
    res_ex = xt.groupby("sample_id").res_norm_ex_self.first().mean()
    write("5-diagnostics.md", f"""
# 5. Reading the diagnostics

The mixture stage writes its checks to `mixmodel/four_pop/diagnostics/`. This page reads them on the example. The general
reference is [`docs/DIAGNOSTICS.md`](../DIAGNOSTICS.md).

## Residuals and `self_share`

The residual measures how well the weighted sources reproduce the target's palette. For the X individuals
`res_norm_ex_self` (the RMSE after dropping the target's own-cluster row) is {res_ex:.4f}
on average, so the fit is good. `self_share` is the part of the palette in the target's own cluster row, {self_share:.2f} on
average: X shares that much IBD with its own members, which no source can reproduce. A large `res_norm_ex_self` for a group of
targets would point to a missing source (see the [missing-source variant](7-variants.md#a-missing-source)).

## Source R flags

`source_R_flags.tsv` compares the total IBD each source emits per individual (R) with the panel median:

{T_srf}

The outgroup {O} (O) emits {o.fold_vs_median:.1f} times less IBD than the median source and gets a `{o.flag}`. That is a statement
about scale, not quality: the weights of a low-R source are deflated. It is computed from the palettes before any fit and is
never a reason to drop a source. Here O is correctly given weight about 0.

## Target R risk

`target_R_flags.tsv` shows, for each target, how much of its estimate rests on flagged sources. For the 12 X individuals
the risk is `none` for {rep["none"]}, `LOW` for {rep["LOW"]}, `MODERATE` for {rep["MODERATE"]} and `HIGH` for {rep["HIGH"]}. The only flagged
source is O, which carries almost no weight, but a tiny raw weight on a source with a low R grows when it is divided by R: the
highest-risk target below ({tt.sample_id}) has a raw weight of {tt.p_warn:.4f} on O, which becomes a corrected share of {tt.q_warn:.3f}. That is an
upper bound on what O could contribute, not a proportion, and is why the tier uses the corrected share.

{T_trf}

`res_ratio` is the target's residual over the median of the targets; `poor_fit` is `yes` at 3 or more.

## Residual diagnostic

`cluster_residuals.tsv` asks, for each target cluster, whether the leftover looks like a population that is not a source. The row of the X cluster ({X}):

{T_cres}

The leftover is small and `miss_pop` is empty: every other population is already a source, so there is no unused population to blame. (The source clusters, which are fitted by themselves, have a residual of 0 and are listed in the file too.) `source_flags.tsv`, `source_sink_by_stratum.tsv` and the
other diagnostic tables are described in [`docs/DIAGNOSTICS.md`](../DIAGNOSTICS.md#residual-diagnostic).

## Knobs

`mixture.r_flag_warn`, `r_flag_severe`, `r_target_*`, `r_fit_gate`, `diag_*`, `mixture.cv`. See
[`docs/CONFIGURATION.md`](../CONFIGURATION.md#mixture).
""")

    # ---------------- source selection
    ma = pd.read_csv(f"{RES}/mixture/auto/mixture_auto.tsv", sep="\t")
    ma["population"] = ma.sample_id.str.split("_").str[0]
    roles = pd.crosstab(ma.population, ma.group)
    roles.index.name = "population"
    srcs = set(ma[ma.group == "source"].population)
    tgts = set(ma[ma.group == "target"].population)
    ab, _ = F.mean_p("auto", "bayesian", "X") if "X" in tgts else F.mean_p("auto", "bayesian", "S1")
    an, _ = F.mean_p("auto", "nnls", "X") if "X" in tgts else F.mean_p("auto", "nnls", "S1")
    if srcs == {"O", "S1", "S2"}:
        outcome = f"""It picked O, S1 and S2 as sources and left X as the target, which is the same split as the `four_pop` panel, and gives X
{ab.get(S1, 0):.3f} S1 (Bayesian) and {an.get(S1, 0):.3f} S1 (NNLS), identical to page 4. The picker takes one source per top-level clade, and
a screen removes populations that look like a mixture of others; neither knows which populations you want as targets."""
    else:
        outcome = f"""It picked {", ".join(sorted(srcs))} as sources and {", ".join(sorted(tgts))} as targets, which differs from the `four_pop` panel.
The picker takes one source per top-level clade and cannot know which populations are admixed targets."""
    write("6-source-selection.md", f"""
# 6. Source selection

Besides the hand-written `four_pop` panel, the workflow builds an automatic panel, `auto`
(`mixmodel/panels/mixture_auto.tsv`). It picks sources from the tree of population distances (page 3): it splits the tree
into clades, screens out populations that look like a mixture of others, and takes the most drifted population per clade.
The sources and targets it chose here:

{md(roles.reset_index())}

{outcome}

Always check the picked sources of an automatic panel, and prefer a hand-written mixture file when you know your sources.

## Knobs

`mixture.auto_k_min`, `auto_k_max`, `auto_source_pick_method` (`tree_spread`, `farthest`, `cluster_medoids`, `differentiated`,
`differentiated_spread`, `differentiated_unadmixed`) and the `auto_source_*` keys. See
[`docs/CONFIGURATION.md`](../CONFIGURATION.md#auto-source-selection).
""")

    # ---------------- variants
    rows = []
    if "raw" in F.mix:
        rb, _ = F.mean_p("raw", "bayesian")
        rn, _ = F.mean_p("raw", "nnls")
    ns = F.mix.get("noS2")
    var = [f"""
# 7. Variants to try

Each variant changes one setting. The outputs are shipped in `example/results/variants/`. The raw-palette variant overwrites the
mixture results, so run it in a copy of the repository; the missing-source variant only adds an output directory. The `deep_split` variant is on [page 2](2-clustering.md).
"""]
    if "raw" in F.mix:
        t = pd.DataFrame({"S1": [real_s1, bay[S1], rb[S1], nn[S1], rn[S1]], "S2": [1 - real_s1, bay[S2], rb[S2], nn[S2], rn[S2]]},
                         index=["realized", "Bayesian, normalized", "Bayesian, raw", "NNLS, normalized", "NNLS, raw"]).round(3)
        t.index.name = "X individuals, mean"
        T_raw = md(t, fmt={"S1": "{:.3f}", "S2": "{:.3f}"}, index=True)
        var.append(f"""
## Raw palettes

```bash
cp -r . ../ibd_adm_raw && cd ../ibd_adm_raw && snakemake --cores 8 --configfile example/variants/raw.yml
```

`mixture.palette_scale: raw` fits source palettes in cM (the target up to a free scale) instead of normalizing each palette to sum 1,
which reduces the over-credit of the source with more IBD per individual:

{T_raw}

Use `raw` only with sources of comparable total sharing ([`docs/DIAGNOSTICS.md`](../DIAGNOSTICS.md#palette-scale)).
""")
    if ns:
        nb, nbsd = F.mean_p("noS2", "bayesian")
        nnn, _ = F.mean_p("noS2", "nnls")
        cr = pd.read_csv(f"{RES}/variants/noS2/diagnostics/cluster_residuals.tsv", sep="\t")
        cx = cr[cr.pop_id == X].iloc[0]
        e_all = F.mix["four_pop"]["bayesian"]
        e_all = e_all[e_all.label == "X"].groupby("sample_id").res_norm_ex_self.first().mean()
        e_no = ns["bayesian"]
        e_no = e_no[e_no.label == "X"].groupby("sample_id").res_norm_ex_self.first().mean()
        T_cr = md(cr[cr.pop_id == X][["pop_id", "res_norm", "top_source", "top_p", "miss_pop", "miss_score"]], fmt={"res_norm": "{:.4f}", "top_p": "{:.3f}", "miss_score": "{:.2f}"})
        var.append(f"""
## A missing source

```bash
cp example/variants/noS2/mixture_noS2.tsv config/panels/default/ && snakemake --cores 8
```

The mixture file `mixture_noS2.tsv` is `mixture_four_pop.tsv` without the S2 individuals: O and S1 are the only sources, but
the S2 cluster stays in the donor palette. This is what a real analysis looks like when an ancestry is missing from the
source set. Outputs go to `mixmodel/noS2/`.

- **The weight lands on the nearest source.** X gets {nb.get(S1, 0):.3f} on S1 (Bayesian; NNLS {nnn.get(S1, 0):.3f}), against {real_s1:.3f} with S2 present.
  The proportion looks confident and is wrong.
- **The residual shows it.** `res_norm_ex_self` of the X individuals rises from {e_all:.4f} to {e_no:.4f}, {e_no / e_all:.0f} times.
- **The residual diagnostic names the missing source.** In `cluster_residuals.tsv` the X cluster ({X}) has `miss_pop` = {cx.miss_pop}
  (the S2 cluster) with `miss_score` {cx.miss_score:.2f}: the leftover correlates with the S2 palette.

{T_cr}

Judge a missing source on the residual, not on the proportion; then add the suspected population to the sources and check that the
residual falls.
""")
    write("7-variants.md", "\n".join(var))
    print("pages written to", PAGES)


# ------------------------------------------------------------------ check
def check(a):
    P = panel_dir(a.run)
    h = cut_height(P)
    ok = True
    c = pd.read_csv(f"{P}/clustering/default.{PFX}.clusters.tsv", sep="\t")
    c = c[c.cut_height.astype(float).round(4) == round(h, 4)].set_index("sample_id").cluster_label
    e = pd.read_csv(f"{RES}/clustering/clusters.tsv", sep="\t").set_index("sample_id").cluster_label
    same = c.sort_index().equals(e.sort_index())
    print("clusters identical:", same)
    ok &= same
    for meth in ("bayesian", "nnls"):
        n = pd.read_csv(f"{P}/mixmodel/four_pop/tables/{PFX}.mixmodel_{meth}.tsv", sep="\t").set_index(["sample_id", "source_pop"]).p
        o = pd.read_csv(f"{RES}/mixture/four_pop/mixmodel_{meth}.tsv", sep="\t").set_index(["sample_id", "source_pop"]).p
        d = (n.sort_index() - o.sort_index()).abs().max()
        print(f"{meth}: max |p - shipped p| = {d:.2e}")
        ok &= d < 1e-6
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    c = sp.add_parser("curate"); c.add_argument("run"); c.add_argument("--raw"); c.add_argument("--noS2"); c.add_argument("--ds3"); c.set_defaults(f=curate)
    sp.add_parser("pages").set_defaults(f=pages)
    k = sp.add_parser("check"); k.add_argument("run"); k.set_defaults(f=check)
    a = ap.parse_args()
    a.f(a)
