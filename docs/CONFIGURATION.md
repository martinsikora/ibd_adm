# Configuration reference

Every knob lives in [`config/config.yml`](../config/config.yml) and is read by
[`workflow/Snakefile`](../workflow/Snakefile).
Values are validated at load time: an out-of-range enum or range raises a
`ValueError` before any job runs. The default column is the code fallback used
when a key is omitted; the shipped `config/config.yml` overrides several of them.

For how to interpret the output columns and diagnostics, see [`DIAGNOSTICS.md`](DIAGNOSTICS.md).

Contents: [Top level](#top-level) · [`ref`](#ref) · [`input_data`](#input_data) ·
[`masking`](#masking) · [`clustering`](#clustering) ·
[`aggregation`](#aggregation) · [colour-map knobs](#colour-map-knobs) ·
[`mixture`](#mixture) · [auto source selection](#auto-source-selection) ·
[`ibd_window_peaks`](#ibd_window_peaks) · [Enable/disable & derived tags](#enabledisable-flags-and-derived-tags)

---

## Top level

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `prefix` | `example_dataset` | string | Dataset/run name used verbatim as `PREFIX` in **every** output filename. Keep it filesystem-safe. |
| `tmpdir` | `/path/to/scratch/ibd_adm` | path or unset | If set, every shell command runs with `TMPDIR` exported here (dir created first). Omit to use the system temp dir. |

## `ref`

| key | default | type | controls |
|-----|---------|------|----------|
| `ref.fasta` | `/path/to/reference/genome.fa` | path | **Reserved.** Read into the config but not used by any rule. Leave it as a placeholder. |
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
| `masking.max_concurrent_ibd_jobs` | `8` | int | Cap on simultaneous masking coverage jobs (`bedtools genomecov`). |
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
| `clustering.standardize_features` | `false` | bool | Z-score each feature (centre and scale). Keep this off with cosine distance: centring turns a low row total into a negative offset in every coordinate, so all low-sharing samples point the same way and cluster together regardless of ancestry. Mutually exclusive with `scale_features`. |
| `clustering.scale_features` | `false` | bool | Divide each feature by its SD without centring. This up-weights low-variance donor columns, as the z-score does, but avoids the offset. Together with the two keys above it sets the transform tag (`raw`/`norm`/`scale`/`scalenorm`/`zscore`/`zscorenorm`). |
| `clustering.cl_size` | `2` | int | `dynamicTreeCut` minimum cluster size. At `2` the cut emits clusters of 1-2 members. These are not resolved sub-populations, and they cause most of the apparent shredding of endogamous groups at finer cuts. |
| `clustering.deep_split` | `3` | int (0–4) | `dynamicTreeCut` `deepSplit` sensitivity. |
| `clustering.knn` | `1` | int | k for the k-NN majority vote that assigns `cluster_min_dist` samples to a cluster (`1` = single nearest neighbour). |
| `clustering.threads` | `24` | int | Threads for the matrix/distance/clustering rules. |
| `clustering.default_panel` | `default` (fallback) | string | Name of the default panel whose clusters seed aggregation. |
| `clustering.base_height` | *(required)* | number | Coarse (fallback) cut height of the default clustering panel. Required when `clustering.enabled` is true; the workflow raises an error at load time if this or `gate_height` is unset. |
| `clustering.gate_height` | *(required)* | number | Fine cut height, at or below `base_height`. If equal to `base_height` there is no gating and the panel is a plain cut at that height. If lower, the panel combines both cuts (`cluster_h{base}g{gate}_{tag}`): a coarse cluster keeps its `base_height` label unless its median per-sample genome-wide IBD reaches `min_sharing`, in which case its finer subdivision is used. Either way there is one panel, not one per height. |
| `clustering.min_sharing` | `250000` | cM ≥ 0 | Median per-sample genome-wide IBD a coarse cluster needs before its finer subdivision is used. This checks that there is enough data, not that the split is good. Unused when `gate_height == base_height`. |

## `aggregation`

Stage 4 (`aggregate_ibd.smk`): per-population IBD sharing + TVD + colour map.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `aggregation.max_concurrent_ibd_jobs` | `8` | int | Global-resource cap on concurrent `aggregate_ibd` jobs. |
| `aggregation.ibd_params.min_l_cm` | `1` | number | Min segment length (cM) for the masked total-IBD pass (stage 2 `ibd_tot`). |
| `aggregation.ibd_params.max_l_cm` | `16` | number | Max segment length (cM). |
| `aggregation.ibd_params.min_lod` | `3` | number | Min LOD/score. |
| `aggregation.full_cluster_pop_overrides` | `{}` | mapping: panel name → options | Per-panel overrides, keyed by a custom panel's directory name under `config/panels/` or by `default` for the raw clustering panel. A panel that is not listed gets none. Do not set this globally, because a pop_id can mean different things in different panels. Options are listed below the table. |
| `aggregation.panels` | `[example_panel]` | list of names | Custom panels; each must be a `config/panels/<name>/` directory. Drives custom aggregation, mixture, PCA, and peaks. Omit/empty to use only the default clustering panels. |

Options of `aggregation.full_cluster_pop_overrides` (per panel):

- `include_recipient_only_pops` (bool, default `false`): treat every pop_id with no `donor_recipient` sample as a full cluster. It stays in the TVD matrix (`tvd_matrix.py`), and its `_r` suffix is dropped in the mixmodel sample_map (`make_mix_sample_map.py`), so it appears under its own name in the TVD tree, PCA and mixture plots.
- `include_pops` (list of pop_ids, default `[]`): the same treatment, by name.
- `exclude_pops` (list of pop_ids, default `[]`): applied after the includes. In the TVD matrix these pops are dropped. In the sample_map they are only removed from the include set and their samples are kept (with `_r`), so a catch-all bin is never silently dropped from the mixture model.

The same pop_id can mean different things in different panels. For example, `unassigned` can be a deliberate catch-all of recipients in one curated panel and an ordinary, mostly `donor_recipient` terminal cluster of `cut_tree` in another.

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
| `aggregation.color_tsne_lc_spread` | `raw` | `raw` \| `rank` | Chroma/luminance spreading mode. With `raw` (min-max), a skewed embedding axis bunches most clusters at the dark, desaturated end and the palette loses its reds and yellows to maroon and olive. |

## `mixture`

Stage 5 (`mixmodel_ibd.smk`): admixture / mixture modelling. Disable with
`enabled: false`. Requires a marker file (`mixture.marker_file` or `ref.marker_file`).

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `mixture.enabled` | `true` | bool | Master switch. |
| `mixture.method` | `nnls` | `nnls` \| `bayesian` \| `both` \| list | Which estimator(s) to run. `both` = NNLS + Bayesian. |
| `mixture.threads` | `12` | int | `future` workers per `run_models` job (parallel across targets). |
| `mixture.seed` | `-1` | int | RNG seed for both estimators (the hybrid slot search and the NNLS jackknife resampling are stochastic too). Negative = unseeded. |
| `mixture.r_flag_warn` | `2.5` | number > 1 | Source-level R scale QC. A source whose emitted IBD R differs from the panel median by at least this fold-change is flagged `WARN`. Validated cohorts start to show a real offset at about 2.5 (Morocco at 5.3x is a confirmed offset; the Steppe/WHG/EEF axis below 2.4x is clean). |
| `mixture.r_flag_severe` | `10.0` | number > `r_flag_warn` | Fold-change at which the flag becomes `SEVERE` (separates those cases from African sources at 20-500x). |
| `mixture.r_target_pmin` | `0.002` | [0, 1) | Per-target R risk: raw weight below which a flagged source cannot raise a target's tier. Stops a weight close to zero divided by a small R from producing a spurious share. |
| `mixture.r_target_high` | `0.02` | (0, 1] | R-corrected share on `SEVERE` sources at which a target is `HIGH` risk. |
| `mixture.r_target_moderate` | `0.05` | (0, 1] | R-corrected share on `WARN` sources at which a target is `MODERATE` risk. |
| `mixture.r_target_low` | `0.01` | (0, 1] | Total R-corrected share on flagged sources at which a target is `LOW` risk. |
| `mixture.cv` | `none` | `none` \| `evenodd` \| `loco` \| `k<K>` \| `test:<chroms>` | Chromosome hold-out CV for the NNLS fits: `evenodd` fits even and scores odd chromosomes and the reverse, `loco` holds out each chromosome in turn, `k<K>` uses K marker-balanced blocks, `test:1,3-5` holds out the listed chromosomes. Writes `<out>.cv.tsv`; the main table is unchanged. Not applied to Bayesian fits; run those by hand with `--cv evenodd --cv_only 1`. See "Chromosome hold-out CV" below. |
| `mixture.marker_file` | (falls back to `ref.marker_file`) | path | Per-chromosome marker counts for IBD length weighting. |

### Mixture output tables

`mixmodel/<set>/tables/*.tsv` has one row per target x source. Columns that are not obvious:

| column | meaning |
|---|---|
| `p`, `se` | Mixture weight and its standard error (NNLS: per-chromosome block jackknife; Bayesian: posterior). `p` is an **IBD share**, not a genome share; see the source R flag below. |
| `res_norm` | Fit residual. **Not comparable between estimators**: Bayesian is an RMSE, NNLS is a plain L2 norm (larger by sqrt of the number of donors). |
| `res_norm_rmse` | The same figure on the RMSE scale for both estimators. |
| `res_norm_ex_self` | RMSE after dropping the target's own-cluster donor row and renormalising both vectors. Use it to compare fits across targets: a large endogamous cohort has a dominant own-cluster row that inflates `res_norm`. `NA` if the target's cluster is not a donor. |
| `self_share` | Fraction of the target's palette in its own cluster row. |
| `self_is_source` | Whether that cluster is itself a source. If TRUE the model can fit the own-cluster column, so `res_norm_ex_self` measures the fit away from home instead of the part of the palette no source can reach. |
| `rhat_median`, `rhat_max` | Chain convergence. Judge it on `rhat_median` and `ess_*`. `rhat_max` becomes large for sources at the simplex corner (near-zero weight) without indicating a problem, and is `Inf` when chains are frozen at different values. `NA` for NNLS rows. |
| `accept_rate*`, `ess_*`, `active_sources_median` | Bayesian sampler diagnostics (`NA` for NNLS). |

**Source R flags** (`diagnostics/<prefix>.source_R_flags.tsv`) mark sources whose total emitted IBD R is far from the panel median. A low R deflates a source's proportions and a high R inflates them. The flag is relative to the panel, is computed before any fit, and is never a reason to drop a source. It gives the direction of a possible scale offset but not its size; no per-source correction was found that generalises. See `workflow/scripts/awk/mixmodel_source_r_flags.awk`.

**Target R risk** (`diagnostics/<prefix>.target_R_flags.tsv`) shows, for each target, how much of its estimate rests on sources flagged in `source_R_flags.tsv`. It is computed from the Bayesian table and the source flags, with no extra IBD pass. The tier uses the R-corrected share `q = (p / R)`, renormalised over the sources. It uses this share because a source with a low R has a deflated raw `p`, so tiering on `p` would miss the worst cases. Dividing by R over-corrects as an estimator, so `q_flagged` is an upper bound on what the flagged sources could contribute. It is not a corrected proportion.

| column | meaning |
|---|---|
| `p_severe`, `p_warn` | Raw weight on `SEVERE` and `WARN` sources. |
| `q_severe`, `q_warn`, `q_flagged` | The same after the R correction (`q_flagged` is their sum). |
| `top_flagged_source` | The flagged source with the largest corrected share (`-` if none). |
| `risk` | `HIGH` if `q_severe >= r_target_high`; else `MODERATE` if `q_warn >= r_target_moderate`; else `LOW` if `q_flagged >= r_target_low`; else `none`. |

### Chromosome hold-out CV (`mixture.cv`)

With `cv` set, the NNLS run also writes `<out>.cv.tsv` next to the main table (`--cv_out` overrides the path), with one row per target and fold. Weights are fitted on the training chromosomes and the target's palette is predicted on the held-out chromosomes, using source palettes rebuilt from those chromosomes.

| column | meaning |
|---|---|
| `fold`, `n_chr_train`, `n_chr_test` | Fold name and chromosome counts. |
| `cm_test` | Held-out IBD (cM, excluding own cluster) the score is based on. |
| `ll_test` | **Ranking statistic.** Held-out multinomial log-likelihood per unit IBD (`sum(y * log(pred))`); higher is better. The prediction is mixed with the mean training palette at weight 1e-3 so a zero-mass donor cannot give -Inf. |
| `res_test_ex_self`, `res_train_ex_self` | RMSE on the held-out / training chromosomes. `res_test - res_train` measures fold size, not overfitting. |
| `n_active`, `p_lost_test`, `n_src_zero_test` | Active sources; weight on sources with no held-out IBD; number of such sources. |

Compare models on the same folds, as a paired per-target difference in `ll_test`. Held-out chromosomes cannot reveal a source that is missing from every chromosome. Relatives among the sources share haplotypes on the held-out chromosomes, so targets with such relatives score optimistically.

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
| `mixture.two_stage_se` | `false` | bool | Bayesian two-stage SE. `false`: `se` is the posterior SD of the single fit, whose likelihood has a fixed 20000 observations; it is too narrow (about 3 times in simulations). `true`: a second fit on the same sources, with `mixture.genome_length_cm` observations, gives the `se`; the weights stay those of the fixed fit. Doubles the Bayesian run time. |
| `mixture.genome_length_cm` | `3500` | number | Length of the genome covered by the IBD data in cM (about 3500 for human autosomes; set it for other species or partial genomes). Used only when `mixture.two_stage_se` is true; it is the number of trials of the likelihood behind the SE, as in SOURCEFIND. |
| `mixture.palette_scale` | `normalized` | `normalized`, `raw` | How palettes are scaled before fitting; see the palette scale section of DIAGNOSTICS.md. `raw` is for re-estimating proportions when the sources differ strongly in total IBD, and is not suitable with single, low-sharing sources. `raw` does not support `mixture.cv`. |
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

Post-hoc, no re-fit. Runs for every mixture panel whenever `bayesian` is in `mixture.method`, and writes:

- `*.cluster_residuals.tsv`: per target cluster, the leftover after the Bayesian fit and the unused population it most resembles. Bayesian fit only.
- `*.source_sink_by_stratum.tsv`: the per-stratum sink test. Bayesian fit only.
- `*.source_flags.tsv`: per-source flags (`absorber`, `poor_fit`, `sink`). Written only when `nnls` is also in `mixture.method`, because it reports the NNLS/Bayesian discordance `disc`, which feeds the `absorber` flag.

The source R flags (`*.source_R_flags.tsv`) are scheduled under the same condition (`bayesian` enabled) because they share the profile step, although they use no model output.

| key | default | type | meaning |
|---|---|---|---|
| `mixture.diag_distal_quantile` | `0.5` | [0, 1] | A source must reach this quantile of the distality distribution (mean TVD from the target mass) before the `absorber` / `poor_fit` flags apply. |
| `mixture.diag_sink_strata` | `12` | int ≥ 1 | Number of strata the target clusters are cut into (ward.D2 on the TVD between palette profiles; no metadata region column is used). With too few, a large stratum becomes a catch-all whose foreign sub-blocks misattribute the flag. With too many, each stratum's dominant source starts tracking its internal cline. |
| `mixture.diag_sink_min_r` | `0.4` | (0, 1] | Correlation with `res_norm`, within a stratum, at which a source counts as a sink: weight that buys down misfit instead of describing ancestry. |
| `mixture.diag_sink_min_n` | `15` | int ≥ 3 | Minimum targets in a stratum before its correlations are trusted. |
| `mixture.diag_sink_min_p` | `0.01` | [0, 1] | Minimum mean weight in the stratum; below this a source has no material influence there. |
| `mixture.diag_sink_min_gap` | `0.15` | [0, 2] | A sink must be the only source over `diag_sink_min_r` in its stratum and lead the runner-up by this margin. Ties mean the stratum contains a badly fitted sub-block that pulls every distant source at once; these are reported as `shared_gradient`. |

Endogamy coupling (`r_endog`, weight against within-cluster IBD enrichment) is reported but is never sufficient alone: a proximate source can legitimately take more weight in the more inbred, less admixed members of its own stratum.

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
  `base_height` and `gate_height` must both be set (equal = no gating): there
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
- `aggregation.color_map_shapes`: non-empty list of integers
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
