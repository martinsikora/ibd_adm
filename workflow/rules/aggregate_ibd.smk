def cluster_agg_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/aggregation"

def cluster_clust_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/clustering"


def panel_agg_dir(panel):
    return f"{PANELS_DIR}/{panel}/aggregation"


def agg_sample_file(panel, height):
    # Default branches use generated panel files per clustering height.
    if panel == "default":
        return f"{cluster_agg_dir(height)}/panels/default.tsv"
    # Custom branches read user-provided panel definitions.
    if panel in CUSTOM_PANELS:
        return f"{PANELS_CFG_DIR}/{panel}/aggregate.tsv"
    raise ValueError(f"Unknown aggregation panel: {panel}")


DEFAULT_AGG_PANELS = []
if ENABLE_DEFAULT_PIPELINE:
    DEFAULT_AGG_PANELS = [
        f"{cluster_agg_dir(_height)}/panels/default.tsv"
        for _height in CLUSTER_HEIGHTS
    ]

AGG_IBD_POP = []
if ENABLE_DEFAULT_PIPELINE:
    AGG_IBD_POP += expand(
        f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
        height=CLUSTER_HEIGHTS,
        chrom=CHROMS,
    )
AGG_IBD_POP += expand(
    f"{panel_agg_dir('{panel}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
    panel=CUSTOM_PANELS,
    chrom=CHROMS,
)

AGG_TVD_OUTPUTS = []
DEFAULT_COLOR_MAPS = []
DEFAULT_TVD_PLOTS = []
if ENABLE_DEFAULT_PIPELINE:
    AGG_TVD_OUTPUTS += [
        f"{cluster_agg_dir(_height)}/tables/{PREFIX}.ibd_pop_tvd.tsv"
        for _height in CLUSTER_HEIGHTS
    ]
    DEFAULT_COLOR_MAPS = [
        f"{cluster_agg_dir(_height)}/panels/default_color_map.tsv"
        for _height in CLUSTER_HEIGHTS
    ]
    DEFAULT_TVD_PLOTS = [
        f"{cluster_agg_dir(_height)}/plots/{PREFIX}.ibd_pop_tvd.pdf"
        for _height in CLUSTER_HEIGHTS
    ] + [
        f"{cluster_agg_dir(_height)}/plots/{PREFIX}.ibd_pop_tvd_heatmap.pdf"
        for _height in CLUSTER_HEIGHTS
    ]
if CUSTOM_PANELS_WITH_COLOR:
    AGG_TVD_OUTPUTS += (
        expand(
            f"{panel_agg_dir('{panel}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
            panel=CUSTOM_PANELS_WITH_COLOR,
        )
        + expand(
            f"{panel_agg_dir('{panel}')}/plots/{PREFIX}.ibd_pop_tvd.pdf",
            panel=CUSTOM_PANELS_WITH_COLOR,
        )
        + expand(
            f"{panel_agg_dir('{panel}')}/plots/{PREFIX}.ibd_pop_tvd_heatmap.pdf",
            panel=CUSTOM_PANELS_WITH_COLOR,
        )
    )

AGG_OUTPUTS = (
    DEFAULT_AGG_PANELS
    + AGG_IBD_POP
    + AGG_TVD_OUTPUTS
    + DEFAULT_COLOR_MAPS
    + DEFAULT_TVD_PLOTS
)


if ENABLE_DEFAULT_PIPELINE:
    rule make_default_agg_panel:
        input:
            clusters=(
                f"{cluster_clust_dir('{height}')}/{CLUSTER_DEFAULT_PANEL}.{PREFIX}.clusters.tsv"
            )
        output:
            panel=(
                f"{cluster_agg_dir('{height}')}/panels/default.tsv"
            )
        shell:
            """
            mkdir -p $(dirname {output.panel})
            python3 workflow/scripts/python/make_agg_panel_from_clusters.py \
                --clusters {input.clusters} \
                --height {wildcards.height} \
                --out {output.panel}
            """


if ENABLE_DEFAULT_PIPELINE:
    rule aggregate_ibd_default:
        input:
            ibd_file=RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_tot.tsv.gz",
            sample_file=lambda wc: agg_sample_file("default", wc.height),
        output:
            f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz"
        resources:
            aggregation_ibd_jobs=1
        priority:
            80
        shell:
            """
            mkdir -p $(dirname {output})
            python3 workflow/scripts/python/aggregate_ibd.py -i {input.ibd_file} -s {input.sample_file} -o {output}
            """


if CUSTOM_PANELS:
    rule aggregate_ibd_custom:
        wildcard_constraints:
            panel="|".join(CUSTOM_PANELS)
        input:
            ibd_file=RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_tot.tsv.gz",
            sample_file=lambda wc: agg_sample_file(wc.panel, "custom"),
        output:
            f"{panel_agg_dir('{panel}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz"
        resources:
            aggregation_ibd_jobs=1
        priority:
            80
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/aggregation/tables
            python3 workflow/scripts/python/aggregate_ibd.py -i {input.ibd_file} -s {input.sample_file} -o {output}
            """


if ENABLE_DEFAULT_PIPELINE:
    rule tvd_default:
        input:
            ibd_files=lambda wc: expand(
                f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                height=wc.height,
                chrom=CHROMS,
            ),
            sample_file=lambda wc: agg_sample_file("default", wc.height),
        output:
            tsv=f"{cluster_agg_dir('{height}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
        params:
            # scoping key "default": the raw clustering panel, keyed by bare
            # tree cluster_label, is a different pop_id namespace from any
            # custom panel -- see the note above full_cluster_pop_flags() in the
            # Snakefile before adding an override here.
            pop_flags=full_cluster_pop_flags("default"),
        priority:
            75
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            python3 workflow/scripts/python/tvd_matrix.py -s {input.sample_file} {params.pop_flags} -o {output.tsv} {input.ibd_files}
            """

    rule default_color_map:
        input:
            tsv=f"{cluster_agg_dir('{height}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
        output:
            cmap=f"{cluster_agg_dir('{height}')}/panels/default_color_map.tsv",
        params:
            chroma_min=COLOR_TSNE_CHROMA_MIN,
            chroma_max=COLOR_TSNE_CHROMA_MAX,
            lum_min=COLOR_TSNE_LUM_MIN,
            lum_max=COLOR_TSNE_LUM_MAX,
            gamma_c=COLOR_TSNE_GAMMA_C,
            gamma_l=COLOR_TSNE_GAMMA_L,
            hue_scale=COLOR_TSNE_HUE_SCALE,
            hue_rotate=COLOR_TSNE_HUE_ROTATE,
            hue_spread=COLOR_TSNE_HUE_SPREAD,
            lc_spread=COLOR_TSNE_LC_SPREAD,
            shapes=COLOR_MAP_SHAPES_ARG,
            embedding=COLOR_MAP_EMBEDDING,
            mapping=COLOR_MAP_MAPPING,
        priority:
            70
        shell:
            """
            mkdir -p $(dirname {output.cmap})
            Rscript workflow/scripts/r/make_color_map_mds.R -i {input.tsv} -o {output.cmap} -k 8 --shapes {params.shapes} --embedding {params.embedding} --mapping {params.mapping} --chroma_min {params.chroma_min} --chroma_max {params.chroma_max} --lum_min {params.lum_min} --lum_max {params.lum_max} --gamma_c {params.gamma_c} --gamma_l {params.gamma_l} --hue_scale {params.hue_scale} --hue_rotate {params.hue_rotate} --hue_spread_mode {params.hue_spread} --lc_spread_mode {params.lc_spread}
            """

    rule tvd_plot_default:
        input:
            tsv=f"{cluster_agg_dir('{height}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
            color_file=f"{cluster_agg_dir('{height}')}/panels/default_color_map.tsv",
        output:
            pdf=f"{cluster_agg_dir('{height}')}/plots/{PREFIX}.ibd_pop_tvd.pdf",
            pdf_hm=f"{cluster_agg_dir('{height}')}/plots/{PREFIX}.ibd_pop_tvd_heatmap.pdf",
        priority:
            65
        shell:
            """
            mkdir -p $(dirname {output.pdf})
            Rscript workflow/scripts/r/tvd_plot.R -i {input.tsv} -c {input.color_file} -p {output.pdf} -m {output.pdf_hm}
            """

if CUSTOM_PANELS_WITH_COLOR:
    rule tvd_custom:
        wildcard_constraints:
            panel="|".join(CUSTOM_PANELS_WITH_COLOR)
        input:
            ibd_files=lambda wc: expand(
                f"{panel_agg_dir('{panel}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                panel=wc.panel,
                chrom=CHROMS,
            ),
            sample_file=lambda wc: agg_sample_file(wc.panel, "custom"),
        output:
            tsv=f"{panel_agg_dir('{panel}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
        params:
            pop_flags=lambda wc: full_cluster_pop_flags(wc.panel),
        priority:
            75
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/aggregation/tables
            python3 workflow/scripts/python/tvd_matrix.py -s {input.sample_file} {params.pop_flags} -o {output.tsv} {input.ibd_files}
            """

    rule tvd_plot_custom:
        wildcard_constraints:
            panel="|".join(CUSTOM_PANELS_WITH_COLOR)
        input:
            tsv=f"{panel_agg_dir('{panel}')}/tables/{PREFIX}.ibd_pop_tvd.tsv",
            color_file=lambda wc: f"{PANELS_CFG_DIR}/{wc.panel}/color_map.tsv",
        output:
            pdf=f"{panel_agg_dir('{panel}')}/plots/{PREFIX}.ibd_pop_tvd.pdf",
            pdf_hm=f"{panel_agg_dir('{panel}')}/plots/{PREFIX}.ibd_pop_tvd_heatmap.pdf",
        priority:
            65
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/aggregation/plots
            Rscript workflow/scripts/r/tvd_plot.R -i {input.tsv} -c {input.color_file} -p {output.pdf} -m {output.pdf_hm}
            """
