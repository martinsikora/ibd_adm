# Configuration reference

Every knob lives in [`config/config.yml`](../config/config.yml) and is read by
[`workflow/Snakefile`](../workflow/Snakefile) (line numbers below refer to it).
Values are validated at load time — an out-of-range enum or range raises a
`ValueError` before any job runs. Defaults shown are the values shipped in this
repo; where a key is optional, the code fallback is given.

Contents: [Top level](#top-level) · [`ref`](#ref) · [`input_data`](#input_data) ·
[`masking`](#masking) · [`clustering`](#clustering) ·
[`aggregation`](#aggregation) · [colour-map knobs](#colour-map-knobs) ·
[`mixture`](#mixture) · [auto source selection](#auto-source-selection) ·
[`ibd_window_peaks`](#ibd_window_peaks) · [Enable/disable & derived tags](#enabledisable-flags-and-derived-tags)

---

> **Defaults.** The *default* column is the value the workflow uses when the key is omitted (the fallback in `workflow/Snakefile`). The shipped `config/config.yml` sets several keys differently (e.g. `mixture.threads`, `mixture.auto_k_min/max`), so a copy of it does not behave like an empty config.

## Top level

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `prefix` | `example_dataset` | string | Dataset/run name used verbatim as `PREFIX` in **every** output filename. Keep it filesystem-safe. |
| `tmpdir` | `/path/to/scratch/ibd_adm` | path or unset | If set, every shell command runs with `TMPDIR` exported here (dir created first). Omit to use the system temp dir. |

## `ref`

| key | default | type | controls |
|-----|---------|------|----------|
| `ref.fasta` | `/path/to/reference/genome.fa` | path | **Reserved** — read into the config dict but not consumed by any rule. Safe to leave as a placeholder. |
| `ref.genome` | `/path/to/reference/genome.genome` | path | Chromosome-length file (`chrom  length`) passed to `bedtools genomecov` during masking. Must cover every chromosome in `ref.chromosomes`. |
| `ref.chromosomes` | `config/chromosomes.txt` | path | Plain list of chromosomes (one per line) → `CHROMS`; fills the `{chrom}` wildcard everywhere. Empty file → error. |
| `ref.marker_file` | `config/n_markers.tsv` | path | Per-chromosome marker counts (`chrom  n`); the mixture model weights IBD by marker density. Fallback for `mixture.marker_file`. |

## `input_data`

| key | default | type | controls |
|-----|---------|------|----------|
| `input_data.ibd` | `resources/ibd_segments/{chrom}.example_dataset.ibdseq.ibd.gz` | path template | Per-chromosome precomputed IBD segments; must keep the `{chrom}` wildcard. Column layout consumed: `$1,$2`=sample ids, `$3`=chrom, `$4,$5`=start,end (bp), `$6`=LOD, `$9`=length (cM). |
| `input_data.individuals` | `config/individuals.tsv` | path | Sample sheet (`sample_id, label, group`); also the sample set for the default clustering panel. |

## `masking`

Stage 1 (`ibd_mask.smk`): build a mask over regions of excess IBD coverage.

| key | default | type | controls |
|-----|---------|------|----------|
| `masking.max_concurrent_ibd_jobs` | `8` | int | Cap on simultaneous masking coverage jobs (IO-heavy `bedtools genomecov`). |
| `masking.ibd_params.min_l_cm` | `2` | number | Min segment length (cM) counted toward coverage. |
| `masking.ibd_params.max_l_cm` | `16` | number | Max segment length (cM). |
| `masking.ibd_params.min_lod` | `3` | number | Min LOD/score. |
| `masking.ibd_params.ibd_trim` | `0.05` | fraction | Fractional trim of the coverage distribution tails before computing the threshold. |
| `masking.ibd_params.ibd_sd` | `5` | number | Windows above `trimmed_mean + ibd_sd · SD` are flagged as masked. |

## `clustering`

Stage 3 (`cluster_ibd.smk`): hierarchical clustering of individuals. Disable with
`enabled: false` (prunes the whole default panel family).

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `clustering.enabled` | `true` | bool | Master switch for the default (auto-clustering) pipeline. |
| `clustering.clust_method` | `ward.D2` | hclust method | Agglomeration method; part of the output cache tag. |
| `clustering.dist_method` | `euclidean` | distance metric | Distance on the IBD feature vectors; part of the cache tag. |
| `clustering.normalize_ibd_vectors` | `false` | bool | L2-normalize each sample's row vector before the distance. |
| `clustering.standardize_features` | `false` | bool | Z-score each feature (centre **and** scale). **Keep off under cosine** — centring turns a low row total into a systematic negative offset in every coordinate, so all low-sharing samples point the same way and cluster together regardless of ancestry. Mutually exclusive with `scale_features`. |
| `clustering.scale_features` | `false` | bool | Divide each feature by its SD **without** centring. Keeps the useful half of the z-score (up-weighting low-variance donor columns) without the offset. Sets the transform tag together with the two above (`raw`/`norm`/`scale`/`scalenorm`/`zscore`/`zscorenorm`). |
| `clustering.cl_size` | `2` | int | `dynamicTreeCut` minimum cluster size. At `2` the cut emits 1–2 member clusters, which are not resolved sub-populations and account for most of the apparent shredding of endogamous groups at finer cuts. |
| `clustering.deep_split` | `3` | int (0–4) | `dynamicTreeCut` `deepSplit` sensitivity. |
| `clustering.knn` | `1` | int | k for the k-NN majority vote that assigns `cluster_min_dist` samples to a cluster (`1` = single nearest neighbour). |
| `clustering.threads` | `24` | int | Threads for the matrix/distance/clustering rules. |
| `clustering.default_panel` | `default` (fallback) | string | Name of the default panel whose clusters seed aggregation. |
| `clustering.base_height` | *(required)* | number | The one default clustering panel's coarse/fallback cut height. Required whenever `clustering.enabled` is true — the workflow raises an error at load time if either this or `gate_height` is unset. |
| `clustering.gate_height` | *(required)* | number | Fine cut height, at or below `base_height`. **Equal to `base_height`**: no gating, the panel is a plain cut at that height. **Below `base_height`**: the panel is the sharing-gated combination of the two (`cluster_h{base}g{gate}_{tag}`) — a coarse cluster keeps its `base_height` label unless its median per-sample genome-wide IBD clears `min_sharing`, in which case its finer `gate_height` subdivision is used. There is exactly one resulting panel either way, not one per height. |
| `clustering.min_sharing` | `250000` | cM ≥ 0 | Median per-sample genome-wide IBD a coarse cluster needs before its finer subdivision is used. A data-sufficiency gate, not a split-quality test. Unused when `gate_height == base_height`. |

## `aggregation`

Stage 4 (`aggregate_ibd.smk`): per-population IBD sharing + TVD + colour map.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `aggregation.max_concurrent_ibd_jobs` | `8` | int | Global-resource cap on concurrent `aggregate_ibd` jobs. |
| `aggregation.ibd_params.min_l_cm` | `1` | number | Min segment length (cM) for the masked total-IBD pass (stage 2 `ibd_tot`). |
| `aggregation.ibd_params.max_l_cm` | `16` | number | Max segment length (cM). |
| `aggregation.ibd_params.min_lod` | `3` | number | Min LOD/score. |
| `aggregation.full_cluster_pop_overrides` | `{}` | mapping: panel name → options | **Scoped per panel** — keyed by a custom panel's directory name under `config/panels/`, or the literal `default` for the raw clustering panel. A panel not listed gets no override at all. Do not set this globally: a pop_id does not mean the same thing in every panel (e.g. `unassigned` can be a deliberate, wholly-recipient catch-all in one manually curated panel and an ordinary, mostly-`donor_recipient` `cut_tree` terminal cluster in another). Each panel's options: `include_recipient_only_pops` (bool, default `false`) treats every pop_id with no `donor_recipient` sample as a full cluster — kept in the TVD matrix (`tvd_matrix.py`) and with its `_r` suffix dropped in the mixmodel sample_map (`make_mix_sample_map.py`), so it appears under its own name in the TVD tree, PCA and mixture plots; `include_pops` (list of pop_ids, default `[]`) applies the same treatment by name; `exclude_pops` (list of pop_ids, default `[]`) is applied after the includes — in the TVD matrix the pops are dropped outright, in the sample_map they are only excluded from the include set (their samples are kept, with `_r`), so a catch-all bin is never silently dropped from the mixture model. |
| `aggregation.panels` | `[example_panel]` | list of names | Custom panels; each must be a `config/panels/<name>/` directory. Drives custom aggregation, mixture, PCA, and peaks. Omit/empty to use only the default clustering panels. |

### Colour-map knobs

Consumed by `make_color_map_mds.R` (`default_color_map` rule) to derive a
perceptually spread colour + shape map from a 3-D embedding of the TVD matrix.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `aggregation.color_map_embedding` | `tsne3` | `tsne3` \| `mds3` | 3-D embedding of the TVD matrix. |
| `aggregation.color_map_mapping` | `radial` | `radial` \| `pca_axes` | How the embedding maps to hue. |
| `aggregation.color_map_shapes` | `[0..19]` | non-empty int list | Plotting shape ids cycled across populations. |
| `aggregation.color_tsne_chroma_min` | `20` | number | Min HCL chroma. |
| `aggregation.color_tsne_chroma_max` | `130` | number | Max HCL chroma (**must exceed** the min). |
| `aggregation.color_tsne_lum_min` | `15` | number | Min HCL luminance. |
| `aggregation.color_tsne_lum_max` | `95` | number | Max HCL luminance (**must exceed** the min). |
| `aggregation.color_tsne_gamma_c` | `0.7` | number > 0 | Chroma gamma. |
| `aggregation.color_tsne_gamma_l` | `0.8` | number > 0 | Luminance gamma (tuned for `lc_spread: rank`; use `0.8` with `raw`). |
| `aggregation.color_tsne_hue_scale` | `1.15` | number | Hue scaling factor. |
| `aggregation.color_tsne_hue_rotate` | `25` | degrees | Hue rotation. |
| `aggregation.color_tsne_hue_spread` | `range` | `raw` \| `range` \| `rank` | Hue spreading mode. |
| `aggregation.color_tsne_lc_spread` | `raw` | `raw` \| `rank` | Chroma/luminance spreading mode. With `raw` (min–max), a skewed embedding axis bunches most clusters at the dark, desaturated end and the palette loses its reds and yellows to maroon and olive. |

## `mixture`

Stage 5 (`mixmodel_ibd.smk`): admixture / mixture modelling. Disable with
`enabled: false`. Requires a marker file (`mixture.marker_file` or `ref.marker_file`).

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `mixture.enabled` | `true` | bool | Master switch. |
| `mixture.method` | `nnls` | `nnls` \| `bayesian` \| `both` \| list | Which estimator(s) to run. `both` = NNLS + Bayesian. |
| `mixture.threads` | `12` | int | `future` workers per `run_models` job (parallel across targets). |
| `mixture.seed` | `-1` | int | RNG seed for both estimators (the hybrid slot search and the NNLS jackknife resampling are stochastic too). Negative = unseeded. |
| `mixture.r_flag_warn` | `2.5` | number > 1 | Source-level R scale QC: a source whose emitted IBD R differs from the panel median by at least this fold-change is flagged `WARN`. 2.5 is where validated cohorts start to show a real offset (Morocco 5.3x is a confirmed offset; the Steppe/WHG/EEF axis at <2.4x is clean). |
| `mixture.r_flag_severe` | `10.0` | number > `r_flag_warn` | Fold-change at which the flag becomes `SEVERE` (separates those cases from African sources at 20-500x). |
| `mixture.cv` | `none` | `none` \| `evenodd` \| `loco` \| `k<K>` \| `test:<chroms>` | Chromosome hold-out CV for the NNLS fits: weights are fitted on the training chromosomes and scored on the held-out ones, writing `<out>.cv.tsv` next to the main table (main table unchanged). `evenodd` = fit even/score odd and the reverse; `loco` = one fold per chromosome; `k<K>` = K marker-balanced blocks; `test:1,3-5` = one custom fold. Rank models on `ll_test` (held-out multinomial log-likelihood per unit IBD) using paired same-fold differences. Held-out chromosomes cannot reveal a source that is missing from every chromosome, and relatives among the sources make targets score optimistically. Not applied to Bayesian fits; run those by hand with `--cv evenodd --cv_only 1`. |
| `mixture.marker_file` | (falls back to `ref.marker_file`) | path | Per-chromosome marker counts for IBD length weighting. |

### Mixture output tables

`mixmodel/<set>/tables/*.tsv` has one row per target x source. Columns that are not obvious:

| column | meaning |
|---|---|
| `p`, `se` | Mixture weight and its standard error (NNLS: per-chromosome block jackknife; Bayesian: posterior). `p` is an **IBD share**, not a genome share; see the source R flag below. |
| `res_norm` | Fit residual. **Not comparable between estimators**: Bayesian is an RMSE, NNLS is a plain L2 norm (larger by sqrt of the number of donors). |
| `res_norm_rmse` | The same figure on the RMSE scale for both estimators. |
| `res_norm_ex_self` | RMSE after dropping the target's own-cluster donor row and renormalising both vectors. Use this to compare fits across targets: a large endogamous cohort has a dominant own-cluster row that inflates `res_norm`. `NA` if the target's cluster is not a donor. |
| `self_share` | Fraction of the target's palette in its own cluster row. |
| `self_is_source` | Whether that cluster is itself a source. If TRUE the model can fit the self column, so `res_norm_ex_self` means "fit away from home" rather than "fit on the part no source can reach". |
| `rhat_median`, `rhat_max` | Chain convergence. Judge convergence on `rhat_median` and `ess_*`; `rhat_max` blows up on sources sitting at the simplex corner (near-zero weight) without indicating a problem. `rhat_max` is `Inf` when chains are frozen at different values. NNLS rows carry `NA`. |
| `accept_rate*`, `ess_*`, `active_sources_median` | Bayesian sampler diagnostics (`NA` for NNLS). |

**Source R flags** (`diagnostics/<prefix>.source_R_flags.tsv`): flags sources whose total emitted IBD R sits far from the panel median, which deflates (low R) or inflates (high R) their proportions. Panel-relative, computed before any fit, and never a reason to drop a source: it flags the risk and direction of a scale offset, not a magnitude, and no per-source correction was found that generalises. See `workflow/scripts/awk/mixmodel_source_r_flags.awk`.

### Chromosome hold-out CV (`mixture.cv`)

With `cv` set, the NNLS run also writes `<out>.cv.tsv` next to the main table (`--cv_out` overrides the path). One row per target x fold: weights are fitted on the training chromosomes and the target's palette is predicted on the held-out chromosomes, with source palettes rebuilt from the held-out chromosomes.

| column | meaning |
|---|---|
| `fold`, `n_chr_train`, `n_chr_test` | Fold name and chromosome counts. |
| `cm_test` | Held-out IBD (cM, excluding own cluster) the score is based on. |
| `ll_test` | **Ranking statistic.** Held-out multinomial log-likelihood per unit IBD (`sum(y * log(pred))`), prediction mixed with the mean training palette at weight 1e-3 so a zero-mass donor cannot give -Inf. Higher is better. |
| `res_test_ex_self`, `res_train_ex_self` | RMSE on the held-out / training chromosomes. `res_test - res_train` measures fold size, not overfitting. |
| `n_active`, `p_lost_test`, `n_src_zero_test` | Active sources; weight on sources with no held-out IBD; number of such sources. |

Compare models on the **same folds** as a paired per-target difference in `ll_test`. Held-out chromosomes cannot reveal a source that is absent from every chromosome, and relatives among the sources share haplotypes on the held-out chromosomes, so such targets score optimistically.

**Bayesian MCMC** (only used when `method` includes `bayesian`):

| key | default | type | controls |
|-----|---------|------|----------|
| `mixture.mcmc_iter` | `200000` | int | Total MCMC iterations. |
| `mixture.burnin` | `20000` | int | Burn-in iterations. |
| `mixture.thin` | `5` | int | Thinning interval. |
| `mixture.proposal_scale` | `200` | number | Dirichlet proposal concentration. |
| `mixture.mcmc_chains` | `1` | int | Number of chains. |
| `mixture.adapt_burnin_frac` | `0.0` | fraction | Fraction of burn-in used for proposal adaptation. |
| `mixture.adapt_interval` | `200` | int | Adaptation update interval. |
| `mixture.adapt_target_accept` | `0.01` | fraction | Target acceptance rate for adaptation. |
| `mixture.local_move_prob` | `1.0` | prob | Probability of a local (vs global) move. |
| `mixture.mean_active_sources` | `6.0` | number | Prior mean number of active sources. |
| `mixture.active_eps` | `1e-4` | number | Threshold below which a source counts as inactive. |
| `mixture.hybrid_active_search` | `1` | 0/1 | Toggle the hybrid active-source search. |
| `mixture.max_active_sources` | `0` | int | Cap on simultaneously active sources (`0` = unlimited). |
| `mixture.active_search_slots` | `200` | int | Candidate slots in the active search. |
| `mixture.active_search_iter` | `1500` | int | Active-search iterations. |
| `mixture.active_search_burnin` | `500` | int | Active-search burn-in. |
| `mixture.active_search_thin` | `10` | int | Active-search thinning. |
| `mixture.active_search_jump_prob` | `0.1` | prob | Active-search jump probability. |

### Auto source selection

Used when a mixture set is `auto` (`make_default_mixture_from_tvd.R`): pick source
populations from the TVD / neighbour-joining tree.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `mixture.auto_k_min` | `2` | int ≥ 2 | Min number of source clusters. |
| `mixture.auto_k_max` | `10` | int ≥ `auto_k_min` | Max number of source clusters (`2 ≤ min ≤ max` enforced). |
| `mixture.auto_source_pick_method` | `tree_spread` | `tree_spread` \| `farthest` \| `cluster_medoids` \| `differentiated` \| `differentiated_spread` \| `differentiated_unadmixed` | Source-picking strategy. `differentiated_unadmixed` screens out populations that look like a mixture of other candidates (triangle-inequality slack + greedy convex-mixture residual), then ranks survivors by drift. |
| `mixture.auto_source_broad_k` | `0` | int | Number of broad clades to partition into (`0` = auto). |
| `mixture.auto_source_max_per_broad_clade` | `1` | int ≥ 1 | Max sources per broad clade. |
| `mixture.auto_source_min_tree_dist_quantile` | `0.0` | [0, 1] | Min pairwise tree-distance quantile between chosen sources. |
| `mixture.auto_source_label_prefix_parts` | `0` | int ≥ 0 | Label prefix parts used for dedup. |
| `mixture.auto_source_max_per_label_prefix` | `0` | int ≥ 0 | Max sources sharing a label prefix. |
| `mixture.auto_source_min_cluster_size` | `1` | int ≥ 1 | (differentiated methods) drop tips below this size. |
| `mixture.auto_source_relative_pendant` | `false` | bool | Score pendant length relative to root-to-tip depth. |
| `mixture.auto_source_admix_slack_quantile` | `0.5` | [0, 1] | (`differentiated_unadmixed`) triangle-slack quantile; lower is stricter about calling a population admixed. |
| `mixture.auto_source_drift_weight` | `0.5` | number ≥ 0 | (`differentiated_unadmixed`) weight of drift vs differentiation when ranking survivors. |

### Residual diagnostic

Post-hoc, no re-fit. Writes `*.cluster_residuals.tsv`, `*.source_flags.tsv` and `*.source_sink_by_stratum.tsv` per mixture panel.

| key | default | type | meaning |
|---|---|---|---|
| `mixture.diag_distal_quantile` | `0.5` | [0, 1] | A source must reach this quantile of the distality distribution (mean TVD from the target mass) before the `absorber` / `poor_fit` flags apply. |
| `mixture.diag_sink_strata` | `12` | int ≥ 1 | Number of strata the target clusters are cut into (ward.D2 on the TVD between palette profiles — data-driven, no metadata region column). Too few and a large stratum becomes a catch-all whose foreign sub-blocks misattribute the flag; too many and each stratum's dominant source starts tracking its internal cline. |
| `mixture.diag_sink_min_r` | `0.4` | (0, 1] | Correlation with `res_norm`, within a stratum, at which a source counts as a sink — weight that buys down misfit rather than describing ancestry. |
| `mixture.diag_sink_min_n` | `15` | int ≥ 3 | Minimum targets in a stratum before its correlations are trusted. |
| `mixture.diag_sink_min_p` | `0.01` | [0, 1] | Minimum mean weight in the stratum; below this a source has no material influence there. |
| `mixture.diag_sink_min_gap` | `0.15` | [0, 2] | A sink must be the **only** source over `diag_sink_min_r` in its stratum and lead the runner-up by this margin. Ties mean the stratum hides a badly-fit sub-block that pulls every distant source at once, reported as `shared_gradient` instead. |

Endogamy coupling (`r_endog`, weight vs within-cluster IBD enrichment) is reported but never sufficient alone: a proximate source legitimately takes more weight in the more inbred, less admixed members of its own stratum.

## `ibd_window_peaks`

Stage 7 (`ibd_window_peaks.smk`): genome-window scan for population-specific IBD
peaks. **Disabled by default.**

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `ibd_window_peaks.enabled` | `false` | bool | Master switch. |
| `ibd_window_peaks.window_size` | `100000` | int (bp) | Scan window size; also names outputs (`<kb>kb_<norm_mode>`). |
| `ibd_window_peaks.norm_mode` | `pop_size` | `none` \| `pop_size` | Per-window normalization. |
| `ibd_window_peaks.methods` | `[kl, jsd, chi2, max_mad_z]` | list | Divergence scores for outlier detection. |
| `ibd_window_peaks.tail_prob` | `0.01` | prob | Tail probability for the outlier threshold. |
| `ibd_window_peaks.threshold_scope` | `genome` | string | Scope over which the threshold is computed. |
| `ibd_window_peaks.min_window_total` | `0` | number | Min per-window total to consider. |
| `ibd_window_peaks.min_cluster_n` | `6` | int | Min cluster size to include. |
| `ibd_window_peaks.contrib_top_n` | `5` | int | Top-N contributing populations reported. |
| `ibd_window_peaks.point_size` | `0.8` | number | Plot point size. |
| `ibd_window_peaks.max_concurrent_coverage_jobs` | `2` | int | Global-resource cap on coverage jobs. |
| `ibd_window_peaks.sample_file` | `config/individuals.tsv` | path | Sample sheet for the peak scan. |
| `ibd_window_peaks.panels` | (falls back to `aggregation.panels`) | list | Panels to scan. |
| `ibd_window_peaks.min_l_cm` / `min_lod` | (fall back to `aggregation.ibd_params`) | number | Segment filters for the scan. |

---

## Enable/disable flags and derived tags

- Setting `clustering.enabled: false` prunes the entire default panel family;
  only custom panels are built. When `clustering.enabled` is true,
  `base_height` and `gate_height` must both be set (equal = no gating) — there
  is exactly one default panel, never zero and never more than one.
- `mixture.enabled: false` prunes all mixture outputs; `ibd_window_peaks.enabled:
  false` prunes the peak scan.
- Clustering intermediates are cached under
  `results/cluster_cache/<tag>/`, where the tag is
  `d<dist_method>_n<transform>_m<clust_method>` (e.g.
  `dcosine_nzscore_mward_D2`). Changing `dist_method`, the feature transforms,
  or `clust_method` writes to a **new** cache directory instead of overwriting,
  and the generated panel is likewise named `cluster_h<base_height>_<tag>`
  (plain) or `cluster_h<base_height>g<gate_height>_<tag>` (gated).

### Validated enums / ranges (raise `ValueError` at load)

- `aggregation.color_map_embedding` ∈ {`tsne3`, `mds3`}
- `aggregation.color_map_mapping` ∈ {`radial`, `pca_axes`}
- `aggregation.color_tsne_hue_spread` ∈ {`raw`, `range`, `rank`}; `aggregation.color_tsne_lc_spread` ∈ {`raw`, `rank`}
- `aggregation.color_tsne_chroma_max` > `..._chroma_min`; `..._lum_max` > `..._lum_min`
- `aggregation.color_tsne_gamma_c` > 0 and `..._gamma_l` > 0
- `aggregation.color_map_shapes` — non-empty list of integers
- `mixture.method` ∈ {`nnls`, `bayesian`, `both`} (or a list of `nnls`/`bayesian`)
- `mixture.auto_source_pick_method` ∈ {`tree_spread`, `farthest`, `cluster_medoids`, `differentiated`, `differentiated_spread`, `differentiated_unadmixed`}
- `2 ≤ mixture.auto_k_min ≤ mixture.auto_k_max`
- `mixture.auto_source_min_tree_dist_quantile` ∈ [0, 1]
- `mixture.auto_source_max_per_broad_clade` ≥ 1; `..._min_cluster_size` ≥ 1
- `clustering.standardize_features` and `clustering.scale_features` are mutually exclusive
- `clustering.base_height`/`gate_height`: both required whenever `clustering.enabled` is true; both numeric; `gate_height` ≤ `base_height`; `min_sharing` ≥ 0 (checked only when they differ)
- `aggregation.full_cluster_pop_overrides.<panel>.include_pops`/`exclude_pops`: lists of pop_ids containing no commas or quotes
- `mixture.diag_sink_strata` ≥ 1; `diag_sink_min_r` ∈ (0, 1]; `diag_sink_min_n` ≥ 3; `diag_sink_min_p` ∈ [0, 1]; `diag_sink_min_gap` ∈ [0, 2]
- `ibd_window_peaks.norm_mode` ∈ {`none`, `pop_size`}
