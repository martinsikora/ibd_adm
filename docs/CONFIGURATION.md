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
| `clustering.heights` | `[0.75, 1.0]` | list of numbers | Adaptive tree-cut heights; **one panel per height**. Empty list also disables the default pipeline. |
| `clustering.clust_method` | `ward.D2` | hclust method | Agglomeration method; part of the output cache tag. |
| `clustering.dist_method` | `cosine` | distance metric | Distance on the IBD feature vectors; part of the cache tag. |
| `clustering.normalize_ibd_vectors` | `false` | bool | L2-normalize each sample's row vector before the distance. |
| `clustering.standardize_features` | `false` | bool | Z-score each feature (centre **and** scale). **Keep off under cosine** — centring turns a low row total into a systematic negative offset in every coordinate, so all low-sharing samples point the same way and cluster together regardless of ancestry. Mutually exclusive with `scale_features`. |
| `clustering.scale_features` | `true` | bool | Divide each feature by its SD **without** centring. Keeps the useful half of the z-score (up-weighting low-variance donor columns) without the offset. Sets the transform tag together with the two above (`raw`/`norm`/`scale`/`scalenorm`/`zscore`/`zscorenorm`). |
| `clustering.cl_size` | `3` | int | `dynamicTreeCut` minimum cluster size. At `2` the cut emits 1–2 member clusters, which are not resolved sub-populations and account for most of the apparent shredding of endogamous groups at finer cuts. |
| `clustering.deep_split` | `3` | int (0–4) | `dynamicTreeCut` `deepSplit` sensitivity. |
| `clustering.knn` | `7` | int | k for the k-NN majority vote that assigns `cluster_min_dist` samples to a cluster (`1` = single nearest neighbour). |
| `clustering.threads` | `32` | int | Threads for the matrix/distance/clustering rules. |
| `clustering.default_panel` | `default` (fallback) | string | Name of the default panel whose clusters seed aggregation. |
| `clustering.gated_cut.enabled` | `false` | bool | Build one labelling from a coarse and a fine cut of the same tree, taking fine labels only where a coarse cluster carries enough IBD to support them. Emits panel `cluster_h{base}g{fine}_{tag}`. |
| `clustering.gated_cut.base_height` | `0.5` | number | Coarse cut height (the fallback label). |
| `clustering.gated_cut.fine_height` | `0.2` | number | Fine cut height; must be **below** `base_height`. Neither height need appear in `heights`. |
| `clustering.gated_cut.min_sharing` | `250000` | cM ≥ 0 | Median per-sample genome-wide IBD a coarse cluster needs before its finer subdivision is used. A data-sufficiency gate, not a split-quality test. |

## `aggregation`

Stage 4 (`aggregate_ibd.smk`): per-population IBD sharing + TVD + colour map.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `aggregation.max_concurrent_ibd_jobs` | `8` | int | Global-resource cap on concurrent `aggregate_ibd` jobs. |
| `aggregation.ibd_params.min_l_cm` | `1` | number | Min segment length (cM) for the masked total-IBD pass (stage 2 `ibd_tot`). |
| `aggregation.ibd_params.max_l_cm` | `16` | number | Max segment length (cM). |
| `aggregation.ibd_params.min_lod` | `3` | number | Min LOD/score. |
| `aggregation.panels` | `[example_panel]` | list of names | Custom panels; each must be a `config/panels/<name>/` directory. Drives custom aggregation, mixture, PCA, and peaks. Omit/empty to use only the default clustering panels. |

### Colour-map knobs

Consumed by `make_color_map_mds.R` (`default_color_map` rule) to derive a
perceptually spread colour + shape map from a 3-D embedding of the TVD matrix.

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `aggregation.color_map_embedding` | `tsne3` | `tsne3` \| `mds3` | 3-D embedding of the TVD matrix. |
| `aggregation.color_map_mapping` | `radial` | `radial` \| `pca_axes` | How the embedding maps to hue. |
| `aggregation.color_map_shapes` | `[0..19]` | non-empty int list | Plotting shape ids cycled across populations. |
| `aggregation.color_tsne_chroma_min` | `40` | number | Min HCL chroma. |
| `aggregation.color_tsne_chroma_max` | `115` | number | Max HCL chroma (**must exceed** the min). |
| `aggregation.color_tsne_lum_min` | `15` | number | Min HCL luminance. |
| `aggregation.color_tsne_lum_max` | `95` | number | Max HCL luminance (**must exceed** the min). |
| `aggregation.color_tsne_gamma_c` | `0.5` | number > 0 | Chroma gamma. |
| `aggregation.color_tsne_gamma_l` | `0.6` | number > 0 | Luminance gamma (tuned for `lc_spread: rank`; use `0.8` with `raw`). |
| `aggregation.color_tsne_hue_scale` | `0.9` | number | Hue scaling factor. |
| `aggregation.color_tsne_hue_rotate` | `25` | degrees | Hue rotation. |
| `aggregation.color_tsne_hue_spread` | `rank` | `raw` \| `range` \| `rank` | Hue spreading mode. |
| `aggregation.color_tsne_lc_spread` | `rank` | `raw` \| `rank` | Chroma/luminance spreading mode. With `raw` (min–max), a skewed embedding axis bunches most clusters at the dark, desaturated end and the palette loses its reds and yellows to maroon and olive. |

## `mixture`

Stage 5 (`mixmodel_ibd.smk`): admixture / mixture modelling. Disable with
`enabled: false`. Requires a marker file (`mixture.marker_file` or `ref.marker_file`).

| key | default | type / allowed | controls |
|-----|---------|----------------|----------|
| `mixture.enabled` | `true` | bool | Master switch. |
| `mixture.method` | `both` | `nnls` \| `bayesian` \| `both` \| list | Which estimator(s) to run. `both` = NNLS + Bayesian. |
| `mixture.threads` | `48` | int | `future` workers per `run_models` job (parallel across targets). |
| `mixture.marker_file` | (falls back to `ref.marker_file`) | path | Per-chromosome marker counts for IBD length weighting. |

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
| `mixture.auto_k_min` | `15` | int ≥ 2 | Min number of source clusters. |
| `mixture.auto_k_max` | `20` | int ≥ `auto_k_min` | Max number of source clusters (`2 ≤ min ≤ max` enforced). |
| `mixture.auto_source_pick_method` | `differentiated_unadmixed` | `tree_spread` \| `farthest` \| `cluster_medoids` \| `differentiated` \| `differentiated_spread` \| `differentiated_unadmixed` | Source-picking strategy. `differentiated_unadmixed` screens out populations that look like a mixture of other candidates (triangle-inequality slack + greedy convex-mixture residual), then ranks survivors by drift. |
| `mixture.auto_source_broad_k` | `0` | int | Number of broad clades to partition into (`0` = auto). |
| `mixture.auto_source_max_per_broad_clade` | `1` | int ≥ 1 | Max sources per broad clade. |
| `mixture.auto_source_min_tree_dist_quantile` | `0.75` | [0, 1] | Min pairwise tree-distance quantile between chosen sources. |
| `mixture.auto_source_label_prefix_parts` | `2` | int ≥ 0 | Label prefix parts used for dedup. |
| `mixture.auto_source_max_per_label_prefix` | `1` | int ≥ 0 | Max sources sharing a label prefix. |
| `mixture.auto_source_min_cluster_size` | `3` | int ≥ 1 | (differentiated methods) drop tips below this size. |
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

- Setting `clustering.enabled: false` (or an empty `heights` list) prunes the
  entire default panel family; only custom panels are built.
- `mixture.enabled: false` prunes all mixture outputs; `ibd_window_peaks.enabled:
  false` prunes the peak scan.
- Clustering intermediates are cached under
  `results/cluster_cache/<tag>/`, where the tag is
  `d<dist_method>_n<transform>_m<clust_method>` (e.g.
  `dcosine_nzscore_mward_D2`). Changing `dist_method`, the feature transforms,
  or `clust_method` writes to a **new** cache directory instead of overwriting,
  and generated panels are likewise named `cluster_h<height>_<tag>`.

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
- `clustering.gated_cut`: `fine_height` < `base_height`, both numeric, `min_sharing` ≥ 0; requires the default clustering pipeline
- `mixture.diag_sink_strata` ≥ 1; `diag_sink_min_r` ∈ (0, 1]; `diag_sink_min_n` ≥ 3; `diag_sink_min_p` ∈ [0, 1]; `diag_sink_min_gap` ∈ [0, 2]
- `ibd_window_peaks.norm_mode` ∈ {`none`, `pop_size`}
