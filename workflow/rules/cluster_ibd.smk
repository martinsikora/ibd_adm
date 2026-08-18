def cluster_heights(panel):
    heights = CLUSTER_PANELS[panel].get("heights")
    if heights is None:
        heights = CLUSTER_HEIGHTS
    return [str(h) for h in heights]

def cluster_clust_dir(height):
    return f"{PANELS_DIR}/{cluster_panel_name(height)}/clustering"


def gated_clust_dir():
    return f"{PANELS_DIR}/{gated_panel_name()}/clustering"


CLUSTER_OUTPUTS = []
if ENABLE_DEFAULT_PIPELINE:
    for _panel in CLUSTER_PANELS:
        for _height in cluster_heights(_panel):
            CLUSTER_OUTPUTS += [
                f"{cluster_clust_dir(_height)}/{_panel}."
                f"{PREFIX}.clusters.tsv",
                f"{cluster_clust_dir(_height)}/{_panel}."
                f"{PREFIX}.clusters_hierarchy.pdf",
                f"{cluster_clust_dir(_height)}/{_panel}."
                f"{PREFIX}.clusters_heatmap.pdf",
            ]

# When GATED_ENABLED, CLUSTER_HEIGHTS is already [CLUSTER_PANEL_HEIGHT] = the
# composite "{base}g{gate}" string (set in the Snakefile), so the loop above
# already added this panel's clusters.tsv + plots under the gated dir -- no
# separate append needed. gated_clust_dir()/gated_panel_name() stay defined
# below for cut_tree_gated's own input/output paths.


if ENABLE_DEFAULT_PIPELINE:
    # Stage 1: dense IBD feature matrix (shared across every config + height)
    rule build_feature_matrix:
        input:
            ibd=IBD_TOT_OUTPUTS,
            sample_file=CLUSTER_PANELS["default"]["sample_file"],
        output:
            rds=CLUSTER_MRAW_RDS,
        threads:
            CLUSTER_THREADS
        resources:
            mem_mb=120000,
            runtime=240,
        shell:
            """
            Rscript workflow/scripts/r/cluster_build_matrix.R \
            -s {input.sample_file} --out {output.rds} --threads {threads} {input.ibd}
            """

    # Stage 2: feature transforms + distance matrix (shared across heights)
    rule distance_matrix:
        input:
            rds=CLUSTER_MRAW_RDS,
        output:
            rds=CLUSTER_MD_RDS,
        params:
            dist_method=CLUSTER_DIST_METHOD,
            normalize_flag=(
                "--normalize_ibd_vectors" if CLUSTER_NORMALIZE_IBD_VECTORS else ""
            ),
            standardize_flag=(
                "--standardize_features" if CLUSTER_STANDARDIZE_FEATURES else ""
            ),
            scale_flag=("--scale_features" if CLUSTER_SCALE_FEATURES else ""),
        threads:
            CLUSTER_THREADS
        resources:
            mem_mb=120000,
            runtime=240,
        shell:
            """
            Rscript workflow/scripts/r/cluster_distance.R \
            --in {input.rds} --out {output.rds} --dist_method {params.dist_method} \
            {params.normalize_flag} {params.standardize_flag} {params.scale_flag} --threads {threads}
            """

    # Stage 3: hierarchical clustering + tree hierarchy (shared across heights)
    rule hclust:
        input:
            rds=CLUSTER_MD_RDS,
            sample_file=CLUSTER_PANELS["default"]["sample_file"],
        output:
            hc=CLUSTER_HC_RDS,
            hier=CLUSTER_HIER_RDS,
        params:
            clust_method=CLUSTER_CLUST_METHOD,
        resources:
            mem_mb=120000,
            runtime=240,
        shell:
            """
            Rscript workflow/scripts/r/cluster_hclust.R \
            --in {input.rds} -s {input.sample_file} --out_hc {output.hc} \
            --out_hier {output.hier} --clust_method {params.clust_method}
            """

    # Stage 4: adaptive tree cut at one height -> clusters.tsv (per height)
    rule cut_tree:
        input:
            hc=CLUSTER_HC_RDS,
            hier=CLUSTER_HIER_RDS,
            dist=CLUSTER_MD_RDS,
            sample_file=lambda wc: CLUSTER_PANELS[wc.panel]["sample_file"],
        output:
            tsv=f"{cluster_clust_dir('{height}')}/{{panel}}.{PREFIX}.clusters.tsv",
        params:
            height=lambda wc: wc.height,
            cl_size=CLUSTER_CL_SIZE,
            deep_split=CLUSTER_DEEP_SPLIT,
            knn=CLUSTER_KNN,
        resources:
            mem_mb=80000,
            runtime=120,
        shell:
            """
            Rscript workflow/scripts/r/cluster_cut.R \
            --hc {input.hc} --hier {input.hier} --dist {input.dist} \
            -s {input.sample_file} --out_tsv {output.tsv} --height {params.height} \
            --cl_size {params.cl_size} --deep_split {params.deep_split} --knn {params.knn}
            """

    # Stage 4b: sharing-gated cut -- one labelling combining the coarse and fine
    # cuts, taking fine labels only where a coarse cluster carries enough IBD to
    # support them. Cheap (reads two clusters.tsv + m_raw row sums), so it is a
    # separate rule rather than folded into cut_tree.
    #
    # The gated panel dir (cluster_h0.5g0.2_TAG) also matches cut_tree's
    # {height} wildcard, so both rules can claim its clusters.tsv. Make the
    # choice explicit rather than relying on resolution order.
    if GATED_ENABLED:
        ruleorder: cut_tree_gated > cut_tree

    if GATED_ENABLED:
        rule cut_tree_gated:
            input:
                base=f"{cluster_clust_dir(GATED_BASE_HEIGHT)}/default.{PREFIX}.clusters.tsv",
                fine=f"{cluster_clust_dir(GATED_FINE_HEIGHT)}/default.{PREFIX}.clusters.tsv",
                matrix=CLUSTER_MRAW_RDS,
            output:
                tsv=f"{gated_clust_dir()}/default.{PREFIX}.clusters.tsv",
            params:
                base_height=GATED_BASE_HEIGHT,
                fine_height=GATED_FINE_HEIGHT,
                min_sharing=GATED_MIN_SHARING,
            resources:
                mem_mb=40000,
                runtime=60,
            shell:
                """
                mkdir -p $(dirname {output.tsv})
                Rscript workflow/scripts/r/cluster_cut_gated.R \
                --base {input.base} --fine {input.fine} --matrix {input.matrix} \
                --base_height {params.base_height} --fine_height {params.fine_height} \
                --min_sharing {params.min_sharing} --out {output.tsv}
                """

    # Stage 5: dendrogram + heatmap plots (per height; isolated so re-cutting
    # never regenerates the large heatmap PDF)
    rule plot_clusters:
        input:
            hc=CLUSTER_HC_RDS,
            matrix=CLUSTER_MRAW_RDS,
            tsv=f"{cluster_clust_dir('{height}')}/{{panel}}.{PREFIX}.clusters.tsv",
            sample_file=lambda wc: CLUSTER_PANELS[wc.panel]["sample_file"],
        output:
            cl=f"{cluster_clust_dir('{height}')}/{{panel}}.{PREFIX}.clusters_hierarchy.pdf",
            hm=f"{cluster_clust_dir('{height}')}/{{panel}}.{PREFIX}.clusters_heatmap.pdf",
        params:
            height=lambda wc: wc.height,
        resources:
            mem_mb=80000,
            runtime=180,
        shell:
            """
            Rscript workflow/scripts/r/cluster_plot.R \
            --hc {input.hc} --matrix {input.matrix} --clusters {input.tsv} \
            -s {input.sample_file} --out_cl {output.cl} --out_hm {output.hm} \
            --height {params.height}
            """
