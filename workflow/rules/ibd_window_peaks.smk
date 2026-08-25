def peak_suffix():
    return f"{PEAK_WINDOW_SIZE // 1000}kb_{PEAK_NORM_MODE}"


def is_default_peak_panel(panel):
    return panel.startswith("cluster_h")


def peak_panel_base(panel):
    return f"{PANELS_DIR}/{panel}/ibd_window_peaks"


def peak_sample_file(panel):
    if is_default_peak_panel(panel):
        return f"{PANELS_DIR}/{panel}/aggregation/panels/default.tsv"
    return f"{PANELS_CFG_DIR}/{panel}/aggregate.tsv"


def peak_color_file(panel):
    if is_default_peak_panel(panel):
        return f"{PANELS_DIR}/{panel}/aggregation/panels/default_color_map.tsv"
    cf = f"{PANELS_CFG_DIR}/{panel}/color_map.tsv"
    return cf if os.path.exists(cf) else None


PEAK_PANELS = []
if PEAK_ENABLED:
    PEAK_PANELS += PEAK_CUSTOM_PANELS
    if ENABLE_DEFAULT_PIPELINE and PEAK_DEFAULT_PANELS:
        PEAK_PANELS += CLUSTER_DEFAULT_PANELS

PEAK_PANELS = sorted(set(PEAK_PANELS))
PEAK_PANELS_WITH_COLOR = [p for p in PEAK_PANELS if peak_color_file(p) is not None]

PEAK_COVERAGE = []
PEAK_WINDOWS = []
PEAK_SCORES = []
PEAK_OUTLIERS = []
PEAK_CONTRIB = []
PEAK_PLOTS = []
PEAK_OUTPUTS = []

if PEAK_ENABLED and PEAK_PANELS:
    PEAK_COVERAGE = expand(
        f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_by_pop.tsv.gz",
        panel=PEAK_PANELS,
        chrom=CHROMS,
    )

    PEAK_WINDOWS = expand(
        f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_window_{peak_suffix()}.tsv.gz",
        panel=PEAK_PANELS,
        chrom=CHROMS,
    )

    PEAK_SCORES = expand(
        f"{peak_panel_base('{panel}')}/tables/genome.scores_{{method}}_{peak_suffix()}.tsv",
        panel=PEAK_PANELS,
        method=PEAK_METHODS,
    )
    PEAK_OUTLIERS = expand(
        f"{peak_panel_base('{panel}')}/tables/genome.outliers_{{method}}_{peak_suffix()}.tsv",
        panel=PEAK_PANELS,
        method=PEAK_METHODS,
    )
    PEAK_CONTRIB = expand(
        f"{peak_panel_base('{panel}')}/tables/genome.contrib_{{method}}_{peak_suffix()}.tsv",
        panel=PEAK_PANELS,
        method=PEAK_METHODS,
    )

    PEAK_PLOTS += expand(
        f"{peak_panel_base('{panel}')}/plots/genome.outlier_scores_{peak_suffix()}.png",
        panel=PEAK_PANELS,
    )
    PEAK_PLOTS += expand(
        f"{peak_panel_base('{panel}')}/plots/genome.coverage_heatmap_{peak_suffix()}.{{chrom}}.png",
        panel=PEAK_PANELS_WITH_COLOR,
        chrom=CHROMS,
    )
    PEAK_PLOTS += expand(
        f"{peak_panel_base('{panel}')}/plots/genome.outlier_contributors_{{method}}_{peak_suffix()}.{{chrom}}.png",
        panel=PEAK_PANELS_WITH_COLOR,
        method=PEAK_METHODS,
        chrom=CHROMS,
    )

    PEAK_OUTPUTS = PEAK_COVERAGE + PEAK_WINDOWS + PEAK_SCORES + PEAK_OUTLIERS + PEAK_CONTRIB + PEAK_PLOTS


if PEAK_ENABLED and PEAK_PANELS:
    rule peak_coverage_by_pop:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS)
        input:
            ibd_file=lambda wc: IBD_INPUT.format(chrom=wc.chrom),
            sample_file=lambda wc: peak_sample_file(wc.panel),
        output:
            f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_by_pop.tsv.gz"
        resources:
            peak_coverage_jobs=1
        priority:
            20
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/tables
            Rscript workflow/scripts/r/ibd_coverage_by_pop.R \
              -i {input.ibd_file} \
              -s {input.sample_file} \
              -o {output} \
              --min_l_cm {PEAK_IBD_MIN_L} \
              --min_lod {PEAK_IBD_MIN_LOD}
            """


if PEAK_ENABLED and PEAK_PANELS:
    rule peak_window_average:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS)
        input:
            cov=f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_by_pop.tsv.gz",
            sample_file=lambda wc: peak_sample_file(wc.panel),
        output:
            f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_window_{peak_suffix()}.tsv.gz"
        params:
            sample_arg=lambda wc, input: f"--sample_file {input.sample_file}" if PEAK_NORM_MODE == "pop_size" else "",
        priority:
            15
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/tables
            Rscript workflow/scripts/r/ibd_window_average.R \
              -i {input.cov} \
              -o {output} \
              -w {PEAK_WINDOW_SIZE} \
              --norm_mode {PEAK_NORM_MODE} \
              {params.sample_arg}
            """


if PEAK_ENABLED and PEAK_PANELS:
    rule peak_outliers:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS),
            method="|".join(PEAK_METHODS),
        input:
            win=lambda wc: expand(
                f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_window_{peak_suffix()}.tsv.gz",
                panel=wc.panel,
                chrom=CHROMS,
            ),
        output:
            out=f"{peak_panel_base('{panel}')}/tables/genome.outliers_{{method}}_{peak_suffix()}.tsv",
            score=f"{peak_panel_base('{panel}')}/tables/genome.scores_{{method}}_{peak_suffix()}.tsv",
            contrib=f"{peak_panel_base('{panel}')}/tables/genome.contrib_{{method}}_{peak_suffix()}.tsv",
        priority:
            10
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/tables
            Rscript workflow/scripts/r/ibd_window_outliers.R \
              --method {wildcards.method} \
              --tail_prob {PEAK_TAIL_PROB} \
              --threshold_scope {PEAK_THRESHOLD_SCOPE} \
              --min_window_total {PEAK_MIN_WINDOW_TOTAL} \
              --min_cluster_n {PEAK_MIN_CLUSTER_N} \
              --contrib_top_n {PEAK_CONTRIB_TOP_N} \
              --contrib_file {output.contrib} \
              --score_file {output.score} \
              -o {output.out} \
              {input.win}
            """


if PEAK_ENABLED and PEAK_PANELS:
    rule peak_plot_scores:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS),
        input:
            scores=lambda wc: expand(
                f"{peak_panel_base('{panel}')}/tables/genome.scores_{{method}}_{peak_suffix()}.tsv",
                panel=wc.panel,
                method=PEAK_METHODS,
            ),
        priority:
            5
        output:
            f"{peak_panel_base('{panel}')}/plots/genome.outlier_scores_{peak_suffix()}.png"
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/plots
            Rscript workflow/scripts/r/plot_ibd_outlier_scores.R \
              --point_size {PEAK_POINT_SIZE} \
              --min_cluster_n {PEAK_MIN_CLUSTER_N} \
              -o {output} \
              {input.scores}
            """


if PEAK_ENABLED and PEAK_PANELS_WITH_COLOR:
    rule peak_merge_windows:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS_WITH_COLOR)
        input:
            lambda wc: expand(
                f"{peak_panel_base('{panel}')}/tables/{{chrom}}.coverage_window_{peak_suffix()}.tsv.gz",
                panel=wc.panel,
                chrom=CHROMS,
            )
        output:
            f"{peak_panel_base('{panel}')}/tables/genome.coverage_window_{peak_suffix()}.tsv.gz"
        priority:
            10
        shell:
            """
            python - <<'PY' {output} {input}
import gzip
import sys

out_file = sys.argv[1]
in_files = sys.argv[2:]

header_written = False
with gzip.open(out_file, "wt") as out:
    for f in in_files:
        with gzip.open(f, "rt") as inp:
            header = inp.readline()
            if not header:
                continue
            if not header_written:
                out.write(header)
                header_written = True
            for line in inp:
                out.write(line)
PY
            """


if PEAK_ENABLED and PEAK_PANELS_WITH_COLOR:
    rule peak_plot_heatmap:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS_WITH_COLOR)
        input:
            in_file=f"{peak_panel_base('{panel}')}/tables/genome.coverage_window_{peak_suffix()}.tsv.gz",
            color_file=lambda wc: peak_color_file(wc.panel),
        output:
            expand(
                f"{peak_panel_base('{{panel}}')}/plots/genome.coverage_heatmap_{peak_suffix()}.{{chrom}}.png",
                chrom=CHROMS,
            )
        params:
            out_base=f"{peak_panel_base('{panel}')}/plots/genome.coverage_heatmap_{peak_suffix()}.png",
            chrom_list=",".join(CHROMS)
        priority:
            5
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/plots
            Rscript workflow/scripts/r/plot_ibd_window_heatmap.R \
              -i {input.in_file} \
              -c {input.color_file} \
              --min_cluster_n {PEAK_MIN_CLUSTER_N} \
              --chrom_list {params.chrom_list} \
              -o {params.out_base}
            """


if PEAK_ENABLED and PEAK_PANELS_WITH_COLOR:
    rule peak_plot_outlier_contributors:
        wildcard_constraints:
            panel="|".join(PEAK_PANELS_WITH_COLOR),
            method="|".join(PEAK_METHODS),
        input:
            in_file=f"{peak_panel_base('{panel}')}/tables/genome.contrib_{{method}}_{peak_suffix()}.tsv",
            x_range_file=f"{peak_panel_base('{panel}')}/tables/genome.coverage_window_{peak_suffix()}.tsv.gz",
            color_file=lambda wc: peak_color_file(wc.panel),
        output:
            expand(
                f"{peak_panel_base('{{panel}}')}/plots/genome.outlier_contributors_{{{{method}}}}_{peak_suffix()}.{{chrom}}.png",
                chrom=CHROMS,
            )
        params:
            out_base=f"{peak_panel_base('{panel}')}/plots/genome.outlier_contributors_{{method}}_{peak_suffix()}.png",
            chrom_list=",".join(CHROMS)
        priority:
            5
        shell:
            """
            mkdir -p {PANELS_DIR}/{wildcards.panel}/ibd_window_peaks/plots
            Rscript workflow/scripts/r/plot_ibd_outlier_contributors.R \
              -i {input.in_file} \
              -x {input.x_range_file} \
              -c {input.color_file} \
              --chrom_list {params.chrom_list} \
              -o {params.out_base}
            """
