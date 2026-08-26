def mix_group_file_default(mix_panel, height):
    if mix_panel == "auto":
        return f"{cluster_mix_dir(height)}/panels/mixture_auto.tsv"
    return f"{PANELS_CFG_DIR}/default/mixture_{mix_panel}.tsv"


def cluster_mix_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/mixmodel"

def cluster_agg_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/aggregation"


def panel_mix_dir(panel):
    return f"{PANELS_DIR}/{panel}/mixmodel"


def custom_mix_file(panel, mix_panel):
    return f"{PANELS_CFG_DIR}/{panel}/mixture_{mix_panel}.tsv"


# custom_color_file is defined in pca_ibd.smk (included before this file)


def mixmodel_tag(method):
    return f"mixmodel_{method}"


# residual diagnostic needs both estimators (nnls loadings + bayesian loadings
# for the discordance flag); only wired when both are produced
MIX_DIAG_ENABLED = (
    MIX_ENABLED and bool(MIX_MARKER_FILE)
    and {"nnls", "bayesian"}.issubset(set(MIX_METHODS))
)


def mix_diag_prefix(mix_dir, mix_panel):
    return f"{mix_dir}/{mix_panel}/diagnostics/{PREFIX}.residual_diagnostic"


def mixmodel_extra_args(method):
    args = f"--method {method}"
    if MIX_SEED >= 0:
        args += f" --seed {MIX_SEED}"
    if method == "bayesian":
        args += (
            f" --mcmc_iter {MIX_MCMC_ITER}"
            f" --burnin {MIX_BURNIN}"
            f" --thin {MIX_THIN}"
            f" --proposal_scale {MIX_PROPOSAL_SCALE}"
            f" --mcmc_chains {MIX_MCMC_CHAINS}"
            f" --adapt_burnin_frac {MIX_ADAPT_BURNIN_FRAC}"
            f" --adapt_interval {MIX_ADAPT_INTERVAL}"
            f" --adapt_target_accept {MIX_ADAPT_TARGET_ACCEPT}"
            f" --local_move_prob {MIX_LOCAL_MOVE_PROB}"
            f" --mean_active_sources {MIX_MEAN_ACTIVE_SOURCES}"
            f" --active_eps {MIX_ACTIVE_EPS}"
            f" --hybrid_active_search {MIX_HYBRID_ACTIVE_SEARCH}"
            f" --max_active_sources {MIX_MAX_ACTIVE_SOURCES}"
            f" --active_search_slots {MIX_ACTIVE_SEARCH_SLOTS}"
            f" --active_search_iter {MIX_ACTIVE_SEARCH_ITER}"
            f" --active_search_burnin {MIX_ACTIVE_SEARCH_BURNIN}"
            f" --active_search_thin {MIX_ACTIVE_SEARCH_THIN}"
            f" --active_search_jump_prob {MIX_ACTIVE_SEARCH_JUMP_PROB}"
        )
    return args


MIXMODEL_OUTPUTS = []
if MIX_ENABLED and MIX_MARKER_FILE:
    MIX_SAMPLE_MAPS = []
    if ENABLE_DEFAULT_PIPELINE:
        MIX_SAMPLE_MAPS += expand(
            f"{cluster_mix_dir('{height}')}/sample_map.tsv",
            height=CLUSTER_HEIGHTS,
        )
    for _panel in CUSTOM_PANELS:
        MIX_SAMPLE_MAPS.append(
            f"{panel_mix_dir(_panel)}/sample_map.tsv"
        )
    MIXMODEL_TABLES = []
    if ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
        MIXMODEL_TABLES += expand(
            f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv",
            height=CLUSTER_HEIGHTS,
            mix_panel=DEFAULT_MIX_PANELS,
            mix_method=MIX_METHODS,
        )
    for _panel, _sets in CUSTOM_MIX_PANELS.items():
        for _mix in _sets:
            for _method in MIX_METHODS:
                MIXMODEL_TABLES.append(
                    f"{panel_mix_dir(_panel)}/{_mix}/tables/{PREFIX}.{mixmodel_tag(_method)}.tsv"
                )

    MIXMODEL_PLOTS = []
    if ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
        MIXMODEL_PLOTS += (
            expand(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.pdf",
                height=CLUSTER_HEIGHTS,
                mix_panel=DEFAULT_MIX_PANELS,
                mix_method=MIX_METHODS,
            )
            + expand(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}_grid.pdf",
                height=CLUSTER_HEIGHTS,
                mix_panel=DEFAULT_MIX_PANELS,
                mix_method=MIX_METHODS,
            )
            + expand(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.source_legend.pdf",
                height=CLUSTER_HEIGHTS,
                mix_panel=DEFAULT_MIX_PANELS,
                mix_method=MIX_METHODS,
            )
        )
    for _panel, _sets in CUSTOM_MIX_PANELS.items():
        if _panel not in CUSTOM_PANELS_WITH_COLOR:
            continue
        for _mix in _sets:
            for _method in MIX_METHODS:
                MIXMODEL_PLOTS.append(
                    f"{panel_mix_dir(_panel)}/{_mix}/plots/{PREFIX}.{mixmodel_tag(_method)}.pdf"
                )
                MIXMODEL_PLOTS.append(
                    f"{panel_mix_dir(_panel)}/{_mix}/plots/{PREFIX}.{mixmodel_tag(_method)}_grid.pdf"
                )
                MIXMODEL_PLOTS.append(
                    f"{panel_mix_dir(_panel)}/{_mix}/plots/{PREFIX}.{mixmodel_tag(_method)}.source_legend.pdf"
                )

    MIXMODEL_DIAG = []
    if MIX_DIAG_ENABLED:
        _diag_suffixes = ("source_flags", "cluster_residuals", "source_sink_by_stratum")
        if ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
            for _s in _diag_suffixes:
                MIXMODEL_DIAG += expand(
                    f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.{_s}.tsv",
                    height=CLUSTER_HEIGHTS,
                    mix_panel=DEFAULT_MIX_PANELS,
                )
        for _panel, _sets in CUSTOM_MIX_PANELS.items():
            for _mix in _sets:
                for _s in _diag_suffixes:
                    MIXMODEL_DIAG.append(
                        f"{panel_mix_dir(_panel)}/{_mix}/diagnostics/{PREFIX}.residual_diagnostic.{_s}.tsv"
                    )
        # source-level R scale flag: derived from the same profiles, no IBD pass
        if ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
            MIXMODEL_DIAG += expand(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.source_R_flags.tsv",
                height=CLUSTER_HEIGHTS,
                mix_panel=DEFAULT_MIX_PANELS,
            )
        for _panel, _sets in CUSTOM_MIX_PANELS.items():
            for _mix in _sets:
                MIXMODEL_DIAG.append(
                    f"{panel_mix_dir(_panel)}/{_mix}/diagnostics/{PREFIX}.source_R_flags.tsv"
                )

    MIXMODEL_OUTPUTS = MIX_SAMPLE_MAPS + MIXMODEL_TABLES + MIXMODEL_PLOTS + MIXMODEL_DIAG


if MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule make_mix_sample_map_default:
        input:
            panel=f"{cluster_agg_dir('{height}')}/panels/default.tsv"
        output:
            f"{cluster_mix_dir('{height}')}/sample_map.tsv"
        params:
            # "default": see the note above full_cluster_pop_flags() in the Snakefile --
            # the raw clustering panel is a different pop_id namespace from any
            # custom panel and must not share its overrides.
            pop_flags=full_cluster_pop_flags("default"),
        priority:
            70
        shell:
            """
            mkdir -p $(dirname {output})
            python3 workflow/scripts/python/make_mix_sample_map.py \
            -i {input.panel} {params.pop_flags} -o {output}
            """

if MIX_MARKER_FILE and CUSTOM_PANELS:
    rule make_mix_sample_map_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS)
        input:
            panel=lambda wc: f"{PANELS_CFG_DIR}/{wc.agg_panel}/aggregate.tsv",
        output:
            f"{panel_mix_dir('{agg_panel}')}/sample_map.tsv"
        params:
            pop_flags=lambda wc: full_cluster_pop_flags(wc.agg_panel),
        priority:
            70
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.agg_panel}/mixmodel
            python3 workflow/scripts/python/make_mix_sample_map.py \
            -i {input.panel} {params.pop_flags} -o {output}
            """

if MIX_ENABLED and MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule run_models_default:
        wildcard_constraints:
            mix_method="|".join(MIX_METHODS),
        input:
            ibd_files=lambda wc: expand(
                f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                height=wc.height,
                chrom=CHROMS,
            ),
            sample_file=f"{cluster_mix_dir('{height}')}/sample_map.tsv",
            individuals=INDIVIDUALS,
            marker_file=MIX_MARKER_FILE,
            group_file=lambda wc: mix_group_file_default(wc.mix_panel, wc.height),
        output:
            f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
        threads:
            MIX_THREADS
        priority:
            90
        params:
            extra_args=lambda wc: mixmodel_extra_args(wc.mix_method)
        shell:
            """
            mkdir -p $(dirname {output})
            Rscript workflow/scripts/r/mixmodel_ibd.R -g {input.group_file} -s {input.sample_file} -i {input.individuals} -l {input.marker_file} -o {output} -t {threads} {params.extra_args} {input.ibd_files}
            """

if MIX_ENABLED and MIX_MARKER_FILE and CUSTOM_PANELS:
    rule run_models_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
            mix_method="|".join(MIX_METHODS),
        input:
            ibd_files=lambda wc: expand(
                f"{PANELS_DIR}/{{agg_panel}}/aggregation/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                agg_panel=wc.agg_panel,
                chrom=CHROMS,
            ),
            sample_file=f"{panel_mix_dir('{agg_panel}')}/sample_map.tsv",
            individuals=INDIVIDUALS,
            marker_file=MIX_MARKER_FILE,
            group_file=lambda wc: custom_mix_file(wc.agg_panel, wc.mix_panel),
        output:
            f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
        threads:
            MIX_THREADS
        ## Memory has to be set here rather than in profiles/slurm because it
        ## depends on mix_method, and --set-resources expressions are evaluated
        ## with only {input, attempt, threads} in scope (snakemake/resources.py
        ## :530) -- a `wildcards` reference there raises NameError, which that
        ## code SWALLOWS, silently assigning the expression string as the value.
        ##
        ## nnls needs more than bayesian, which is counter-intuitive but
        ## measured: on 2026-08-20 all three nnls fits were killed at 96 GB with
        ## MaxRSS 96.8-97.9 GB, while the bayesian fit that completed peaked at
        ## 74.3 GB. The jackknife resamples over 22 leave-one-chromosome-out
        ## matrices, and those are held alongside the genome-wide ones.
        ## Floors are set where each method actually lands so the retry ladder
        ## is insurance, not the mechanism -- climbing from 48 GB cost three
        ## dead attempts per fit last time.
        resources:
            ## 2026-08-22: bayesian 96000 -> 80000. The only true peak we have
            ## is slurm MaxRSS 74.3 GB, from the fit that completed at 48
            ## threads; at 16 threads the mid-MCMC footprint is 34-48 GB. It is
            ## NOT cut further because peak RSS is the IBD read at startup, not
            ## the MCMC -- that read is single-threaded and loads the same 22
            ## tables whatever `threads` is, so the peak does not scale with
            ## thread count and 74.3 GB may well recur at 16 threads. 80000
            ## clears the known peak by 8%.
            ##
            ## nnls stays at 144000 until measured. Its only datapoint is the
            ## OUT_OF_MEMORY kill at 96 GB (MaxRSS 96.8-97.9), which is a lower
            ## bound on the peak, not the peak. The three nnls fits run locally
            ## on 2026-08-21/22 succeeded but local runs have no MaxRSS;
            ## workflow/scripts/tune/peak_rss.py is sampling world_base_2 nnls
            ## for a real figure. Do not lower this on the strength of "it
            ## worked locally" -- the local node has 755 GB and never capped it.
            mem_mb=lambda wc, attempt: (144000 if wc.mix_method == "nnls" else 80000) * attempt,
            ## 2026-08-21: 2880 -> 1800. At mixture.threads 16 and mcmc_iter
            ## 100000 a bayesian fit is ~17.4h (19501 targets x 4 chains x
            ## 100k iter x 125us measured, / 16 threads), so 48h was a ~2.8x
            ## over-request. That matters here because fairshare is 0.059 and
            ## priority ~5.9k against a pending median of 19.3k -- backfill is
            ## the realistic way in, and it favours short walltime. 30h keeps
            ## ~70% margin; retries escalate to 60h/90h if that is ever wrong.
            runtime=lambda wc, attempt: 1800 * attempt,
        priority:
            90
        params:
            extra_args=lambda wc: mixmodel_extra_args(wc.mix_method)
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.agg_panel}/mixmodel/{wildcards.mix_panel}/tables
            Rscript workflow/scripts/r/mixmodel_ibd.R -g {input.group_file} -s {input.sample_file} -i {input.individuals} -l {input.marker_file} -o {output} -t {threads} {params.extra_args} {input.ibd_files}
            """

if MIX_ENABLED and MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule make_default_mixture_auto:
        input:
            tvd=f"{cluster_agg_dir('{height}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
            sample_file=f"{cluster_mix_dir('{height}')}/sample_map.tsv",
        output:
            mix=f"{cluster_mix_dir('{height}')}/panels/mixture_auto.tsv",
        params:
            k_min=MIX_AUTO_K_MIN,
            k_max=MIX_AUTO_K_MAX,
            source_pick_method=MIX_AUTO_SOURCE_PICK_METHOD,
            source_broad_k=MIX_AUTO_SOURCE_BROAD_K,
            source_max_per_broad_clade=MIX_AUTO_SOURCE_MAX_PER_BROAD_CLADE,
            source_min_tree_dist_quantile=MIX_AUTO_SOURCE_MIN_TREE_DIST_QUANTILE,
            source_label_prefix_parts=MIX_AUTO_SOURCE_LABEL_PREFIX_PARTS,
            source_max_per_label_prefix=MIX_AUTO_SOURCE_MAX_PER_LABEL_PREFIX,
            source_min_cluster_size=MIX_AUTO_SOURCE_MIN_CLUSTER_SIZE,
            source_relative_pendant=MIX_AUTO_SOURCE_RELATIVE_PENDANT,
            source_admix_slack_quantile=MIX_AUTO_SOURCE_ADMIX_SLACK_QUANTILE,
            source_drift_weight=MIX_AUTO_SOURCE_DRIFT_WEIGHT,
        priority:
            85
        shell:
            """
            mkdir -p $(dirname {output.mix})
            Rscript workflow/scripts/r/make_default_mixture_from_tvd.R -i {input.tvd} -s {input.sample_file} \
              --k_min {params.k_min} --k_max {params.k_max} --source_pick_method {params.source_pick_method} \
              --source_broad_k {params.source_broad_k} --source_max_per_broad_clade {params.source_max_per_broad_clade} \
              --source_min_tree_dist_quantile {params.source_min_tree_dist_quantile} \
              --source_label_prefix_parts {params.source_label_prefix_parts} --source_max_per_label_prefix {params.source_max_per_label_prefix} \
              --source_min_cluster_size {params.source_min_cluster_size} {params.source_relative_pendant} \
              --source_admix_slack_quantile {params.source_admix_slack_quantile} \
              --source_drift_weight {params.source_drift_weight} \
              -o {output.mix}
            """

if MIX_ENABLED and MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule plot_models_default:
        wildcard_constraints:
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            sample_file=f"{cluster_mix_dir('{height}')}/sample_map.tsv",
            color_file=(
                f"{cluster_agg_dir('{height}')}/panels/default_color_map.tsv"
            ),
        output:
            f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.pdf"
        priority:
            60
        shell:
            """
            mkdir -p $(dirname {output})
            Rscript workflow/scripts/r/plot_mixmodel.R -i {input.model} -s {input.sample_file} -c {input.color_file} -o {output}
            """

if MIX_ENABLED and MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule plot_models_grid_default:
        wildcard_constraints:
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            sample_file=f"{cluster_mix_dir('{height}')}/sample_map.tsv",
            color_file=(
                f"{cluster_agg_dir('{height}')}/panels/default_color_map.tsv"
            ),
        output:
            f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}_grid.pdf"
        priority:
            60
        shell:
            """
            mkdir -p $(dirname {output})
            Rscript workflow/scripts/r/plot_mixmodel.R -i {input.model} -s {input.sample_file} -c {input.color_file} -o {output} --source_grid
            """

if MIX_ENABLED and MIX_MARKER_FILE and ENABLE_DEFAULT_PIPELINE:
    rule plot_source_legend_default:
        wildcard_constraints:
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            color_file=(
                f"{cluster_agg_dir('{height}')}/panels/default_color_map.tsv"
            ),
        output:
            f"{cluster_mix_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.source_legend.pdf"
        priority:
            60
        shell:
            """
            mkdir -p $(dirname {output})
            Rscript workflow/scripts/r/plot_source_legend.R -i {input.model} -c {input.color_file} -o {output}
            """

if MIX_ENABLED and MIX_MARKER_FILE and CUSTOM_PANELS_WITH_COLOR:
    rule plot_models_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS_WITH_COLOR),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            sample_file=f"{panel_mix_dir('{agg_panel}')}/sample_map.tsv",
            color_file=lambda wc: custom_color_file(wc.agg_panel),
        output:
            f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.pdf"
        priority:
            60
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.agg_panel}/mixmodel/{wildcards.mix_panel}/plots
            Rscript workflow/scripts/r/plot_mixmodel.R -i {input.model} -s {input.sample_file} -c {input.color_file} -o {output}
            """

    rule plot_models_grid_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS_WITH_COLOR),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            sample_file=f"{panel_mix_dir('{agg_panel}')}/sample_map.tsv",
            color_file=lambda wc: custom_color_file(wc.agg_panel),
        output:
            f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}_grid.pdf"
        priority:
            60
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.agg_panel}/mixmodel/{wildcards.mix_panel}/plots
            Rscript workflow/scripts/r/plot_mixmodel.R -i {input.model} -s {input.sample_file} -c {input.color_file} -o {output} --source_grid
            """

    rule plot_source_legend_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS_WITH_COLOR),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
            mix_method="|".join(MIX_METHODS),
        input:
            model=(
                f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_{{mix_method}}.tsv"
            ),
            color_file=lambda wc: custom_color_file(wc.agg_panel),
        output:
            f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/plots/{PREFIX}.mixmodel_{{mix_method}}.source_legend.pdf"
        priority:
            60
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.agg_panel}/mixmodel/{wildcards.mix_panel}/plots
            Rscript workflow/scripts/r/plot_source_legend.R -i {input.model} -c {input.color_file} -o {output}
            """


## ==========================================================================
## residual diagnostic: flag deep sources behaving as bad proxies (absorber /
## poor_fit) and name the population each is standing in for. Post-hoc on the
## existing nnls + bayesian model outputs; no model re-fit. A streaming awk pass
## over the panel's ibd_pop tables builds the genome-total palette / source /
## validation profiles, then the R script reconstructs residuals and projects
## the unexplained sharing onto every unused population.

if MIX_DIAG_ENABLED and CUSTOM_PANELS:
    rule mixmodel_residual_profiles_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
        input:
            ibd_files=lambda wc: expand(
                f"{PANELS_DIR}/{{agg_panel}}/aggregation/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                agg_panel=wc.agg_panel,
                chrom=CHROMS,
            ),
            group_file=lambda wc: custom_mix_file(wc.agg_panel, wc.mix_panel),
        output:
            palette=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/pop_prof.tsv",
            src=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            val=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/val_prof.tsv",
            valids=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/val.ids",
            srcn=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/src_n.tsv",
        priority:
            55
        shell:
            """
            d=$(dirname {output.palette})
            mkdir -p "$d"
            awk -F'\\t' 'NR>1 && $2=="source"{{print $1}}' {input.group_file} > "$d/src.ids"
            awk -F'\\t' 'NR>1 && $2=="target"{{print $1}}' {input.group_file} | head -3 > {output.valids}
            zcat {input.ibd_files} | gawk -v SRCF="$d/src.ids" -v VALF={output.valids} \
              -v POPF={output.palette} -v SRCPF={output.src} -v VALPF={output.val} \
              -v SRCNF={output.srcn} \
              -f workflow/scripts/awk/ibd_residual_profiles.awk
            """

    rule mixmodel_source_r_flags_custom:
        # Source-level R QC flag. Reads only the profiles written by
        # mixmodel_residual_profiles_custom, so it costs no IBD pass and cannot
        # affect any fit. See workflow/scripts/awk/mixmodel_source_r_flags.awk for
        # what a tier does and does not mean.
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
        input:
            src=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            srcn=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/src_n.tsv",
        output:
            flags=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/{PREFIX}.source_R_flags.tsv",
        params:
            panel=lambda wc: f"{wc.agg_panel}/{wc.mix_panel}",
            warn=MIX_R_FLAG_WARN,
            severe=MIX_R_FLAG_SEVERE,
        priority:
            55
        shell:
            """
            gawk -v NF_FILE={input.srcn} -v PANEL={params.panel} \
              -v WARN={params.warn} -v SEVERE={params.severe} \
              -f workflow/scripts/awk/mixmodel_source_r_flags.awk {input.src} > {output.flags}
            """

    rule mixmodel_residual_diagnostic_custom:
        wildcard_constraints:
            agg_panel="|".join(CUSTOM_PANELS),
            mix_panel="|".join(CUSTOM_MIX_ALL) if CUSTOM_MIX_ALL else ".*",
        input:
            palette=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/pop_prof.tsv",
            src=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            val=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/val_prof.tsv",
            valids=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/val.ids",
            sample_file=f"{panel_mix_dir('{agg_panel}')}/sample_map.tsv",
            nnls=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_nnls.tsv",
            bayesian=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_bayesian.tsv",
        output:
            flags=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.source_flags.tsv",
            resid=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.cluster_residuals.tsv",
            sink=f"{panel_mix_dir('{agg_panel}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.source_sink_by_stratum.tsv",
        params:
            prefix=lambda wc: mix_diag_prefix(panel_mix_dir(wc.agg_panel), wc.mix_panel),
            distal_q=MIX_DIAG_DISTAL_QUANTILE,
            sink_strata=MIX_DIAG_SINK_STRATA,
            sink_min_r=MIX_DIAG_SINK_MIN_R,
            sink_min_n=MIX_DIAG_SINK_MIN_N,
            sink_min_p=MIX_DIAG_SINK_MIN_P,
            sink_min_gap=MIX_DIAG_SINK_MIN_GAP,
        priority:
            55
        shell:
            """
            Rscript workflow/scripts/r/mixmodel_residual_diagnostic.R \
              {input.palette} {input.src} {input.val} {input.sample_file} \
              {input.nnls} {input.bayesian} {input.valids} {params.prefix} {params.distal_q} \
              {params.sink_strata} {params.sink_min_r} {params.sink_min_n} \
              {params.sink_min_p} {params.sink_min_gap}
            """

if MIX_DIAG_ENABLED and ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
    rule mixmodel_residual_profiles_default:
        wildcard_constraints:
            mix_panel="|".join(DEFAULT_MIX_PANELS),
        input:
            ibd_files=lambda wc: expand(
                f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                height=wc.height,
                chrom=CHROMS,
            ),
            group_file=lambda wc: mix_group_file_default(wc.mix_panel, wc.height),
        output:
            palette=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/pop_prof.tsv",
            src=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            val=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/val_prof.tsv",
            valids=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/val.ids",
            srcn=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/src_n.tsv",
        priority:
            55
        shell:
            """
            d=$(dirname {output.palette})
            mkdir -p "$d"
            awk -F'\\t' 'NR>1 && $2=="source"{{print $1}}' {input.group_file} > "$d/src.ids"
            awk -F'\\t' 'NR>1 && $2=="target"{{print $1}}' {input.group_file} | head -3 > {output.valids}
            zcat {input.ibd_files} | gawk -v SRCF="$d/src.ids" -v VALF={output.valids} \
              -v POPF={output.palette} -v SRCPF={output.src} -v VALPF={output.val} \
              -v SRCNF={output.srcn} \
              -f workflow/scripts/awk/ibd_residual_profiles.awk
            """

    rule mixmodel_source_r_flags_default:
        # Default-panel twin of mixmodel_source_r_flags_custom.
        wildcard_constraints:
            mix_panel="|".join(DEFAULT_MIX_PANELS),
        input:
            src=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            srcn=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/src_n.tsv",
        output:
            flags=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.source_R_flags.tsv",
        params:
            panel=lambda wc: f"h{wc.height}/{wc.mix_panel}",
            warn=MIX_R_FLAG_WARN,
            severe=MIX_R_FLAG_SEVERE,
        priority:
            55
        shell:
            """
            gawk -v NF_FILE={input.srcn} -v PANEL={params.panel} \
              -v WARN={params.warn} -v SEVERE={params.severe} \
              -f workflow/scripts/awk/mixmodel_source_r_flags.awk {input.src} > {output.flags}
            """

    rule mixmodel_residual_diagnostic_default:
        wildcard_constraints:
            mix_panel="|".join(DEFAULT_MIX_PANELS),
        input:
            palette=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/pop_prof.tsv",
            src=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/src_prof.tsv",
            val=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/val_prof.tsv",
            valids=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/val.ids",
            sample_file=f"{cluster_mix_dir('{height}')}/sample_map.tsv",
            nnls=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_nnls.tsv",
            bayesian=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.mixmodel_bayesian.tsv",
        output:
            flags=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.source_flags.tsv",
            resid=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.cluster_residuals.tsv",
            sink=f"{cluster_mix_dir('{height}')}/{{mix_panel}}/diagnostics/{PREFIX}.residual_diagnostic.source_sink_by_stratum.tsv",
        params:
            prefix=lambda wc: mix_diag_prefix(cluster_mix_dir(wc.height), wc.mix_panel),
            distal_q=MIX_DIAG_DISTAL_QUANTILE,
            sink_strata=MIX_DIAG_SINK_STRATA,
            sink_min_r=MIX_DIAG_SINK_MIN_R,
            sink_min_n=MIX_DIAG_SINK_MIN_N,
            sink_min_p=MIX_DIAG_SINK_MIN_P,
            sink_min_gap=MIX_DIAG_SINK_MIN_GAP,
        priority:
            55
        shell:
            """
            Rscript workflow/scripts/r/mixmodel_residual_diagnostic.R \
              {input.palette} {input.src} {input.val} {input.sample_file} \
              {input.nnls} {input.bayesian} {input.valids} {params.prefix} {params.distal_q} \
              {params.sink_strata} {params.sink_min_r} {params.sink_min_n} \
              {params.sink_min_p} {params.sink_min_gap}
            """
