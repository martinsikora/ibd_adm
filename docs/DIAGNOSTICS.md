# Fit diagnostics

How to judge a mixture fit, from the numbers the workflow writes. Column and key names
are given as they appear in the output files; the config keys are described in
[`CONFIGURATION.md`](CONFIGURATION.md).

Contents: [What is fitted](#what-is-fitted) · [Estimates and uncertainty](#estimates-and-uncertainty) ·
[Residuals](#residuals) · [Sampler diagnostics](#sampler-diagnostics-bayesian-only) ·
[Source R flags](#source-r-flags) · [Target R risk](#target-r-risk) ·
[Residual diagnostic](#residual-diagnostic) · [Chromosome hold-out CV](#chromosome-hold-out-cv) ·
[Common problems](#common-problems) · [Checklist](#checklist-for-a-run)

## What is fitted

Every individual has a palette: the total IBD (cM) it shares with each donor population.
For a target, the model finds non-negative weights `p` that sum to 1 so that the weighted
sum of the source palettes reproduces the target palette. How the palettes are scaled is
set by `mixture.palette_scale` (see [Palette scale](#palette-scale)): with the default
`raw`, the source palettes are mean per-individual palettes in cM and the target palette
is fitted up to a free overall scale; with `normalized`, every palette is first divided
by its own total. Two estimators are available:

- `nnls`: sum-to-one non-negative least squares. One point estimate per target, with a
  jackknife standard error.
- `bayesian`: a SOURCEFIND-style MCMC with a Dirichlet prior on the weights, an
  optional Poisson prior on the number of active sources (`mixture.mean_active_sources`),
  and a search over which sources are active.

Two things follow from this and matter for every number below.

1. With `palette_scale: normalized`, `p` is a share of **IBD**, not of the genome. A source
   that emits more IBD into the panel per unit of ancestry (a higher R, see
   [Source R flags](#source-r-flags)) gets more weight than its ancestry alone justifies,
   and the reverse for a low-R source. With the default `raw` this dependence on the
   total IBD of the sources is removed, but not every offset is: a source that is only a
   relative of the true source still shares less IBD with the target than with itself.
   Compare `p` between targets fitted with the same sources, and be careful comparing it
   across sources.
2. A palette can only be reproduced from the sources you supply. If an ancestry is
   missing from the source set, its signal is absorbed by the nearest available source.
   Most of the diagnostics below exist to make that visible.

Which file holds what:

| file | one row per | estimators |
|---|---|---|
| `tables/*.mixmodel_<method>.tsv` | target and source | either |
| `tables/*.mixmodel_<method>.cv.tsv` | target and fold | NNLS in the workflow; Bayesian by hand |
| `diagnostics/*.residual_diagnostic.cluster_residuals.tsv` | target cluster | Bayesian |
| `diagnostics/*.residual_diagnostic.source_sink_by_stratum.tsv` | stratum and source | Bayesian |
| `diagnostics/*.residual_diagnostic.source_flags.tsv` | source | both |
| `diagnostics/*.source_R_flags.tsv` | source | none (uses profiles only) |
| `diagnostics/*.target_R_flags.tsv` | target | Bayesian |

## Palette scale

`mixture.palette_scale` (command line `--palette_scale`) has two settings.

- `raw` (default): each source is the mean per-individual palette of its source
  individuals, in cM. No donor-count correction is applied: the aggregation sums every
  individual's within-cluster entry over the n - 1 other cluster members and its
  between-cluster entries over a random subset of n - 1 of the n donors, so source
  individuals and targets are compared on the same number of donors. The target palette
  is fitted up to a free overall scale and the weights are normalised afterwards, so they
  are ancestry fractions. In the Bayesian model the prediction is rescaled to sum to 1
  before the likelihood.
- `normalized`: every palette is divided by its own total before fitting, which was the
  only behaviour before this option existed. A mixture of normalised palettes weights a
  source by its ancestry share times its total IBD per individual, so a source that
  carries more total IBD is over-credited. The size of this effect depends on the ratio
  of the total IBD of the sources over the donor panel, and changes with the donor set.

In simulations with known ancestry (two sources, one admixed target, 22 chromosomes) the
`normalized` fit overestimated the share of the source with more total IBD by 0.03 to 0.04
at mixed targets, and `raw` changed this to a small bias of 0.01 to 0.02 in the other
direction, from sources that are relatives of the true sources and not the sources
themselves. Results from earlier runs used `normalized`; set `palette_scale: normalized`
to reproduce them. `raw` does not support `--cv`, and the workflow stops with a message if
`mixture.cv` is set together with it.

## Estimates and uncertainty

| column | meaning |
|---|---|
| `p` | Weight of the source in the target. Sums to 1 over the sources of a target. |
| `se` | NNLS: leave-one-chromosome-out jackknife, weighted by chromosome size. Bayesian: posterior standard deviation. By default it comes from a likelihood with a fixed 20000 observations and is about 3 times narrower than the error between individuals in simulations. With `mixture.two_stage_se: true` (and `mixture.genome_length_cm` set to the genome length), a second fit on the same sources gives a wider posterior (about 1.3 times too narrow in simulations) and the weights are unchanged. In scenarios with many close-relative candidate sources that wide posterior, used for the weights as well, shifted them toward the middle of the simplex, which is why only the SE is taken from it. |
| `active_sources_median` | Bayesian: median number of sources with weight above `mixture.active_eps` (default 1e-4) per posterior draw. |
| `selected_sources_n` | Bayesian: number of sources passed to the continuous sampler after the active-source search. Equal to all sources when the search is off or `max_active_sources` covers them all. |

How to read them:

- A weight is only meaningful relative to its `se`. Treat `p / se` below about 3 as
  indistinguishable from zero for reporting purposes. In one panel, sources the Bayesian
  fit did not select sat at `p / se` around 0.6 to 0.8 (the prior leaves them slightly
  above zero), while selected sources were at 3 or more. That gap is what you look for.
- Bayesian weights are sparse. A source with a tiny `p` and a tiny `se` is not a small
  contribution measured precisely; it is a source the sampler mostly switched off.
- NNLS weights are exactly 0 for sources the fit does not need, so its `se` for those is
  uninformative. It is more useful to look at the sources with non-zero `p` and check that
  their jackknife `se` is small relative to `p`.
- The two estimators can disagree when sources are similar to each other, because
  collinear sources trade weight. Compare the sum over a group of related sources
  instead of each member. See [Collinear sources](#common-problems).
- The standard errors reflect sampling noise in the IBD, not error from a missing or
  mis-specified source.

## Residuals

The residual measures how well the fitted mixture reproduces the target palette. All
variants are root mean squared differences between the observed and predicted palette
fractions, so smaller is better and values are only comparable within a panel.

| column | meaning |
|---|---|
| `res_norm` | Bayesian: RMSE over all donors. NNLS: the plain L2 norm from the solver, larger than the RMSE by the square root of the number of donors. Not comparable between estimators. |
| `res_norm_rmse` | The same figure on the RMSE scale for both estimators. Use this to compare NNLS and Bayesian. |
| `res_norm_ex_self` | RMSE after dropping the target's own-cluster donor row and renormalising both vectors. `NA` if the target's cluster is not a donor. |
| `self_share` | Fraction of the target's palette that lies in its own cluster row. |
| `self_is_source` | Whether the target's own cluster is also a source. |

Reading them:

- `res_norm` is dominated by the own-cluster row for large endogamous cohorts, which no
  source can reproduce. This makes big, tightly related groups look badly fitted.
  `res_norm_ex_self` removes that row and is the figure to compare across targets.
- `self_share` tells you how large that effect is for a given target. A high value with a
  low `res_norm_ex_self` means the fit is good and the target is inbred.
- If `self_is_source` is true, the model can fit the own-cluster column, so
  `res_norm_ex_self` measures the fit away from home. It does not measure the part of
  the palette no source can reach, which is what it means otherwise.
- A large `res_norm_ex_self` for a group of targets, with a structured leftover, points to a
  missing source. The [residual diagnostic](#residual-diagnostic) looks for it.
- Adding a source can only lower the in-sample NNLS residual, so a small improvement alone
  does not justify keeping a source. Use the hold-out CV to check whether it generalises.

## Sampler diagnostics (Bayesian only)

These columns are `NA` for NNLS rows.

| column | meaning |
|---|---|
| `accept_rate`, `_min`, `_max` | Proposal acceptance rate, averaged over chains, and the lowest and highest chain. |
| `proposal_scale_final` | Dirichlet proposal concentration at the end (median over chains). Adapted during burn-in when `adapt_burnin_frac > 0`. |
| `ess_min`, `ess_median` | Effective sample size per source, summed over chains, then the minimum and median over sources. |
| `rhat_median`, `rhat_max` | Gelman-Rubin statistic per source across chains, then the median and maximum. Needs `mcmc_chains` of 2 or more; `NA` with one chain. |
| `n_keep`, `n_chains` | Posterior draws kept in total, and number of chains. |

Reading them:

- **Acceptance.** The workflow's default `adapt_target_accept` is low (0.01), because
  Dirichlet moves in a sparse simplex are mostly rejected. A low acceptance rate is
  expected and is not a problem in itself. A rate of exactly 0 in a chain, or a large gap
  between `accept_rate_min` and `accept_rate_max`, suggests a chain got stuck.
- **R-hat.** Judge convergence on `rhat_median`. Values close to 1 (below 1.1 is the usual
  working limit) mean the chains agree. `rhat_max` is dominated by sources sitting at the
  simplex corner, where the weight is near zero in every chain and the statistic
  becomes unstable without any real problem. It is `Inf` when chains are frozen at different
  values, which is a real failure (each chain stayed at its own starting solution).
- **ESS.** `ess_min` can be set by a near-zero source and is then pessimistic.
  Look at `ess_median` for the sources that matter, and at the sources with non-negligible
  `p`. As a general rule of thumb, a few hundred effective draws are enough for a stable
  posterior mean; this is a convention, not a threshold tuned for this pipeline.
- **What to do about poor sampling.** More iterations (`mcmc_iter`), more chains, a longer
  burn-in, or fewer sources. Persistent disagreement between chains with a well-fitted
  residual can mean several nearly equivalent source combinations.
- Convergence is separate from correctness. A well-converged chain can still describe a
  mis-specified model.

## Source R flags

`diagnostics/<prefix>.source_R_flags.tsv`, one row per source.

R is the total IBD an individual emits into the donor panel, averaged over the members of
the source. Because `p` is an IBD share (roughly ancestry times R), a source with a much
lower R than the panel median has its weight deflated, and a much higher R inflated. The
flag is computed from the panel alone, before looking at any fit.

| column | meaning |
|---|---|
| `component`, `n` | Source and number of individuals in it. |
| `R_excl_self` | Mean IBD emitted per individual, excluding sharing inside the source itself (the reported statistic). |
| `R_genomewide` | The same including within-source sharing, for reference. |
| `fold_vs_median` | R divided by the panel median, or its inverse if lower, so always 1 or more. |
| `direction` | `underestimated` (low R, weight deflated) or `overestimated`; `-` if the source is `ok`. |
| `flag` | `ok`, `WARN` at `mixture.r_flag_warn` (2.5×) or more, `SEVERE` at `mixture.r_flag_severe` (10×) or more. |

Reading it:

- The tiers are relative to the panel. Adding donors of one ancestry changes every other
  source's tier without any change in those sources. Do not compare flags between panels.
- The flag gives the direction of a possible scale offset, not its size. Validation
  against published qpAdm results found offsets for one and the same flagged source ranging
  from 1.3× to 6.8× across cohorts, and no per-source correction that predicted which one
  applies. Do not rescale weights by R.
- A flag is not a statement about source quality. A source with the lowest R in a panel
  can still be the only proxy for its ancestry and receive well-supported weights. Its
  ranking of targets can be excellent while its scale is off. Never drop a source only
  because of its R.
- The main use is to know which comparisons are safe. Comparing one target's weight on a
  flagged source with another target's weight on the same source is reasonable. Comparing
  it with the weight on an unflagged source is not.

## Target R risk

`diagnostics/<prefix>.target_R_flags.tsv`, one row per target.

Shows how much of each target's estimate rests on flagged sources. It uses the R-corrected
share `q = (p / R)`, renormalised over the sources, because tiering on the raw `p` would
miss the worst cases (a badly deflated source has a small `p`). Dividing by R over-corrects
as an estimator, so `q_flagged` is an upper bound on what the flagged sources could
contribute. It is a screen, not a corrected proportion.

| column | meaning |
|---|---|
| `p_severe`, `p_warn` | Raw weight on `SEVERE` and `WARN` sources. |
| `q_severe`, `q_warn`, `q_flagged` | The same after the R correction; `q_flagged` is their sum. |
| `top_flagged_source` | Flagged source with the largest corrected share, or `-`. |
| `risk` | `HIGH` if `q_severe >= r_target_high` (0.02); else `MODERATE` if `q_warn >= r_target_moderate` (0.05); else `LOW` if `q_flagged >= r_target_low` (0.01); else `none`. |

A source whose raw weight is below `r_target_pmin` (0.002) cannot raise the tier.

Reading it: `HIGH` means a `SEVERE` source could plausibly carry a real contribution that
the raw weights hide, so treat that target's proportions on those sources as a direction,
not a value. `none` means the flagged sources do not matter for this target. Most targets
in a panel without deep African or other extreme-R sources come out as `none`.

## Residual diagnostic

Post-hoc, no refit. It reconstructs the target and source palettes, works on target
**clusters** (using the mean Bayesian weights of each cluster), and asks whether the
unexplained part of the palette looks like a population that is not among the sources.
Its parameters are the `mixture.diag_*` keys.

### `cluster_residuals.tsv`

| column | meaning |
|---|---|
| `pop_id` | Target cluster. |
| `res_norm` | RMSE of the leftover (observed minus predicted) over all donors, including the cluster's own row. Not the same as `res_norm_ex_self` in the model table. |
| `top_source`, `top_p` | Source with the largest mean weight, and that weight. |
| `miss_pop`, `miss_score` | The unused population whose palette correlates best (Pearson r) with the leftover, and that correlation. |
| `stratum` | Stratum the cluster was assigned to (see below). |
| `n_samples`, `endogamy` | Cluster size, and the excess of within-cluster IBD relative to its size (1 is the panel average, above 1 is more inbred). |

`miss_score` is the useful column. A high value (roughly 0.3 or more is where the flags
start to react) means the leftover is structured and resembles a specific unused
population: a candidate missing source. A low value with a large `res_norm` means noise or
drift without an obvious candidate. The `miss_pop` is a lead to check, not proof, because
populations related to the true one score similarly.

### `source_flags.tsv`

Written only when both estimators are run. One row per source.

| column | meaning |
|---|---|
| `distality`, `is_distal` | Mean TVD between the source and the target clusters (weighted by cluster size), and whether it reaches the `mixture.diag_distal_quantile` quantile (default 0.5) of all sources. |
| `n_loaded`, `mean_p_loaded` | Number of target clusters where the source has weight of 0.15 or more, and their mean weight. |
| `mean_res_loaded`, `mean_miss_score` | Mean residual and mean `miss_score` over the loaded clusters. |
| `coupling_p_res` | Correlation, across target clusters, between the source's weight and the cluster residual. |
| `modal_miss_pop` | The `miss_pop` that occurs most often among the loaded clusters. |
| `disc` | Mean absolute difference between the NNLS and Bayesian cluster weights for this source. |
| `sink_strata`, `sink_stratum`, `max_r_res`, `max_r_endog` | From the stratum test below. |
| `n_eff_strata`, `top_stratum`, `top_stratum_share` | Clade confinement, below. |
| `absorber`, `poor_fit`, `sink`, `flag_type` | The flags; `flag_type` is the first that applies in the order absorber, poor_fit, sink. |

`absorber` and `poor_fit` apply to distal sources loaded in at least 3 clusters:

- `absorber`: the residual is not worse where the source is loaded (`coupling_p_res <= 0.05`),
  and either the leftover there is structured (`mean_miss_score >= 0.30`) or the
  estimators disagree (`disc >= 0.08`). The source fits, but it is redundant or standing in
  for something else.
- `poor_fit`: the residual is worse where the source is loaded (`coupling_p_res > 0.05`),
  the leftover is structured (`mean_miss_score >= 0.30`), and the mean residual is above
  the typical cluster. A deep source forced onto ancestry it cannot span.
- `sink`: see the stratum test.

A flag is a prompt to inspect the source, not a verdict. Look at `modal_miss_pop` for what
it may be standing in for, and try adding that population as a source.

### Stratum test: `source_sink_by_stratum.tsv`

A source can be a real ancestry component or an unfitted intercept that absorbs drift,
and the two are indistinguishable across the whole panel: a sink is quiet wherever it
has competition and only acts where the panel has no close source. So the target clusters
are cut into strata (`mixture.diag_sink_strata`, default 12, from the palette geometry),
and within each stratum the source's weight is correlated with the cluster residual
(`r_res`) and with endogamy (`r_endog`).

A source is a `sink` in a stratum when its weight rises with the residual (`r_res` at
`mixture.diag_sink_min_r`, default 0.4), it has non-trivial weight there, the stratum has
enough targets, it is the only source over that bar, and it leads the runner-up by
`diag_sink_min_gap` (0.15). If several sources rise together the stratum contains a badly
fitted sub-block that pulls every distant source at once; that is reported as
`shared_gradient` for the stratum and is not a sink. Correlation with endogamy alone is
not enough (a close source legitimately takes more weight in the more inbred members of
its stratum); it only corroborates (`sink_via` says whether it was misfit alone or
misfit plus drift).

### Clade confinement

`n_eff_strata` is the effective number of strata a source's palette spreads over, with
`top_stratum` and `top_stratum_share`. A deep source that is genuinely ancestral to one
clade shares broadly inside that clade and little outside, giving a low `n_eff_strata`. An
unanchored profile is smeared across many. Low distality with few strata is a deep source
doing real work; many strata suggests an intercept. In one panel this separated two
equally old, equally distal lineages: one tracked a real ancestry axis and the other
absorbed drift in the most inbred clusters.

## Chromosome hold-out CV

`<out>.cv.tsv`, written when `mixture.cv` is set. Fits use the training chromosomes and are
scored on the held-out ones, with source palettes rebuilt from the held-out chromosomes.
This guards against overfitting the weights.

| column | meaning |
|---|---|
| `fold`, `n_chr_train`, `n_chr_test` | Fold and chromosome counts. |
| `cm_test` | Held-out IBD (cM, excluding own cluster) the score rests on. |
| `ll_test` | Held-out multinomial log-likelihood per unit IBD, `sum(y * log(pred))`. Higher is better. The prediction is mixed with the mean training palette at weight 1e-3 so a zero-mass donor cannot give minus infinity. |
| `res_test_ex_self`, `res_train_ex_self` | RMSE on the held-out and training chromosomes. |
| `n_active`, `p_lost_test`, `n_src_zero_test` | Sources with weight above 1e-4; weight on sources with no held-out IBD; number of such sources. |

How to use it:

- Rank models on `ll_test`, as the paired per-target difference between two models on the
  **same folds**. The mean difference is in nats per unit IBD. A positive value favours the
  first model; look at the fraction of targets that improve as well as the mean.
- Do not compare `res_test_ex_self` with `res_train_ex_self` to measure overfitting. A
  held-out block is a small part of the genome and its palette is sparser and noisier, so
  the gap mostly reflects fold size.
- `p_lost_test` above zero means the model puts weight on sources that have no IBD in the
  held-out chromosomes, which happens with very small sources.
- Held-out chromosomes cannot reveal a source that is missing from every chromosome.
  Relatives among the sources share haplotypes across chromosomes, so targets with
  relatives in the source set score optimistically.

## Common problems

| symptom | check | likely cause | what to do |
|---|---|---|---|
| Large weight on the most basal or oldest source, with a structured residual | `cluster_residuals` `miss_pop` and `miss_score` | A source is missing; its ancestry lands on the nearest, often the deepest, one | Add the suspected source and see if the residual falls. Judge on the residual, not on the proportion |
| High `res_norm` for big groups, low `res_norm_ex_self` | `self_share`, `endogamy` | Own-cluster row of an endogamous cohort | Compare on `res_norm_ex_self` |
| `rhat_max` very large, `rhat_median` near 1 | `p` of the offending sources | Sources at the simplex corner | Ignore; use `rhat_median` and `ess_median` |
| `rhat_max` is `Inf` | Trace of the chains | Chains frozen at different values | Rerun with more iterations or check the proposal scale |
| Weights differ between NNLS and Bayesian | `disc`, group sums | Collinear sources trading weight | Compare the total over the related group; remove or merge near-duplicate sources |
| A minor component follows genome quality | Correlation of weight with `gp_avg_target`; z per individual | Absorbed noise from low-quality samples | Do not report it as ancestry; exclude low-quality sources; never filter on depth of imputed data |
| Weight rises with residual inside a stratum | `source_sink_by_stratum` | Source acting as a sink for drift | Add a closer source for that stratum, or drop the sink |
| Deep source gets near-zero weight | Its palette against the targets | Gap too deep for IBD to convert into a proportion | Expect this beyond a few thousand years; use closer proxies and treat the deep member as a hypothesis |
| `HIGH` in `target_R_flags` | `top_flagged_source` | Target relies on an extreme-R source | Report direction only |
| Proportion changes when a source is added or removed | Sum over related sources | Collinear alternatives | Report the summed component, not each member |

## Checklist for a run

1. Sampler: `rhat_median` near 1 and `ess_median` reasonable for the sources that matter
   (Bayesian).
2. Residuals: compare `res_norm_ex_self` within the panel; look for groups of targets with
   large, structured residuals.
3. Residual diagnostic: look at `miss_pop` for those groups, then at flagged sources and
   any `sink` strata.
4. Source R flags: note which sources are `WARN` or `SEVERE`, and check `target_R_flags` for
   the targets that rely on them.
5. Weights: read them with their `se`, and as sums over groups of related sources.
6. Optional: run the hold-out CV to rank alternative source sets.
7. Change the source set only for a reason you can name, then rerun and check that the
   residual and the flags moved as expected.
