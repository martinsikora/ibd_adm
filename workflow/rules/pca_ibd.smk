def cluster_pca_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/pca"

def cluster_agg_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/aggregation"

def cluster_mix_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/mixmodel"


def panel_pca_dir(panel):
    return f"{PANELS_DIR}/{panel}/pca"

def default_color_file(height):
    return f"{cluster_agg_dir(height)}/panels/default_color_map.tsv"


def default_mix_group_file(height, mix_panel):
    if mix_panel == "auto":
        return f"{cluster_mix_dir(height)}/panels/mixture_auto.tsv"
    return f"{PANELS_CFG_DIR}/default/mixture_{mix_panel}.tsv"

def default_sample_map(height):
    return f"{cluster_mix_dir(height)}/sample_map.tsv"


def custom_color_file(panel):
    return f"{PANELS_CFG_DIR}/{panel}/color_map.tsv"


def custom_mix_group_file(panel, mix_panel):
    return f"{PANELS_CFG_DIR}/{panel}/mixture_{mix_panel}.tsv"

def custom_sample_map(panel):
    return f"{PANELS_DIR}/{panel}/mixmodel/sample_map.tsv"


CUSTOM_PCA_PANELS = [
    p for p in CUSTOM_PANELS_WITH_COLOR
    if CUSTOM_MIX_PANELS.get(p)
]

PCA_OUTPUTS = []
if ENABLE_DEFAULT_PIPELINE:
    PCA_OUTPUTS += expand(
        f"{cluster_pca_dir('{height}')}/default/tables/{PREFIX}.pca.tsv",
        height=CLUSTER_HEIGHTS,
    )
    PCA_OUTPUTS += expand(
        f"{cluster_pca_dir('{height}')}/default/plots/{PREFIX}.pca.pdf",
        height=CLUSTER_HEIGHTS,
    )
    if DEFAULT_MIX_PANELS:
        PCA_OUTPUTS += expand(
            f"{cluster_pca_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.pca_proj.tsv",
            height=CLUSTER_HEIGHTS,
            mix_panel=DEFAULT_MIX_PANELS,
        )
        PCA_OUTPUTS += expand(
            f"{cluster_pca_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.pca_proj.pdf",
            height=CLUSTER_HEIGHTS,
            mix_panel=DEFAULT_MIX_PANELS,
        )

if CUSTOM_PCA_PANELS:
    PCA_OUTPUTS += expand(
        f"{panel_pca_dir('{panel}')}/default/tables/{PREFIX}.pca.tsv",
        panel=CUSTOM_PCA_PANELS,
    )
    PCA_OUTPUTS += expand(
        f"{panel_pca_dir('{panel}')}/default/plots/{PREFIX}.pca.pdf",
        panel=CUSTOM_PCA_PANELS,
    )
    for _panel, _sets in CUSTOM_MIX_PANELS.items():
        if _panel not in CUSTOM_PANELS_WITH_COLOR:
            continue
        for _mix in _sets:
            PCA_OUTPUTS.append(
                f"{panel_pca_dir(_panel)}/{_mix}/tables/{PREFIX}.pca_proj.tsv"
            )
            PCA_OUTPUTS.append(
                f"{panel_pca_dir(_panel)}/{_mix}/plots/{PREFIX}.pca_proj.pdf"
            )


if ENABLE_DEFAULT_PIPELINE:
    rule pca_default:
        input:
            ibd_files=lambda wc: expand(
                f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                height=wc.height,
                chrom=CHROMS,
            ),
            sample_file=lambda wc: default_sample_map(wc.height),
            color_file=lambda wc: default_color_file(wc.height),
        output:
            tsv=f"{cluster_pca_dir('{height}')}/default/tables/{PREFIX}.pca.tsv",
            pdf=f"{cluster_pca_dir('{height}')}/default/plots/{PREFIX}.pca.pdf",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {output.pdf})
            Rscript workflow/scripts/r/pca_ibd.R -s {input.sample_file} -c {input.color_file} -o {output.tsv} -p {output.pdf} {input.ibd_files}
            """


if ENABLE_DEFAULT_PIPELINE and DEFAULT_MIX_PANELS:
    rule pca_default_projected:
        input:
            ibd_files=lambda wc: expand(
                f"{cluster_agg_dir('{height}')}/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                height=wc.height,
                chrom=CHROMS,
            ),
            sample_file=lambda wc: default_sample_map(wc.height),
            group_file=lambda wc: default_mix_group_file(wc.height, wc.mix_panel),
            color_file=lambda wc: default_color_file(wc.height),
        output:
            tsv=f"{cluster_pca_dir('{height}')}/{{mix_panel}}/tables/{PREFIX}.pca_proj.tsv",
            pdf=f"{cluster_pca_dir('{height}')}/{{mix_panel}}/plots/{PREFIX}.pca_proj.pdf",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            mkdir -p $(dirname {output.pdf})
            Rscript workflow/scripts/r/pca_ibd.R -s {input.sample_file} -g {input.group_file} -c {input.color_file} -o {output.tsv} -p {output.pdf} --project {input.ibd_files}
            """


if CUSTOM_PCA_PANELS:
    rule pca_custom:
        wildcard_constraints:
            panel="|".join(CUSTOM_PCA_PANELS)
        input:
            ibd_files=lambda wc: expand(
                f"{PANELS_DIR}/{{panel}}/aggregation/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                panel=wc.panel,
                chrom=CHROMS,
            ),
            sample_file=lambda wc: custom_sample_map(wc.panel),
            color_file=lambda wc: custom_color_file(wc.panel),
        output:
            tsv=f"{panel_pca_dir('{panel}')}/default/tables/{PREFIX}.pca.tsv",
            pdf=f"{panel_pca_dir('{panel}')}/default/plots/{PREFIX}.pca.pdf",
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/pca/default/tables
            mkdir -p {PANELS_DIR}/{wildcards.panel}/pca/default/plots
            Rscript workflow/scripts/r/pca_ibd.R -s {input.sample_file} -c {input.color_file} -o {output.tsv} -p {output.pdf} {input.ibd_files}
            """

    if CUSTOM_MIX_ALL:
        rule pca_custom_projected:
            wildcard_constraints:
                panel="|".join(CUSTOM_PCA_PANELS),
                mix_panel="|".join(CUSTOM_MIX_ALL),
            input:
                ibd_files=lambda wc: expand(
                    f"{PANELS_DIR}/{{panel}}/aggregation/tables/{{chrom}}.{PREFIX}.ibd_pop.tsv.gz",
                    panel=wc.panel,
                    chrom=CHROMS,
                ),
                sample_file=lambda wc: custom_sample_map(wc.panel),
                group_file=lambda wc: custom_mix_group_file(wc.panel, wc.mix_panel),
                color_file=lambda wc: custom_color_file(wc.panel),
            output:
                tsv=f"{panel_pca_dir('{panel}')}/{{mix_panel}}/tables/{PREFIX}.pca_proj.tsv",
                pdf=f"{panel_pca_dir('{panel}')}/{{mix_panel}}/plots/{PREFIX}.pca_proj.pdf",
            shell:
                """
                mkdir -p {PANELS_DIR}/{wildcards.panel}/pca/{wildcards.mix_panel}/tables
                mkdir -p {PANELS_DIR}/{wildcards.panel}/pca/{wildcards.mix_panel}/plots
                Rscript workflow/scripts/r/pca_ibd.R -s {input.sample_file} -g {input.group_file} -c {input.color_file} -o {output.tsv} -p {output.pdf} --project {input.ibd_files}
                """
