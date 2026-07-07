IBD_COVERAGE = expand(
    RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_coverage.tsv.gz",
    chrom=CHROMS,
)

IBD_MASK_BEDS = expand(
    RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_mask.bed",
    chrom=CHROMS,
)

IBD_MASK_PLOTS = expand(
    RESULTS_DIR + "/masking/plots/{chrom}." + PREFIX + ".ibd_mask.pdf",
    chrom=CHROMS,
)

MASK_OUTPUTS = IBD_COVERAGE + IBD_MASK_BEDS + IBD_MASK_PLOTS


rule get_genomecov:
    input:
        ibd=IBD_INPUT,
        inds=INDIVIDUALS,
        genome=GENOME,
    output:
        RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_coverage.tsv.gz"
    params:
        min_l=MASK_MIN_L,
        max_l=MASK_MAX_L,
        min_lod=MASK_MIN_LOD,
    shell:
        """
        set -euo pipefail
        mkdir -p {RESULTS_DIR}/masking/tables
        gzip -cd {input.ibd} | tail -n+2 | \
        awk 'FNR==NR {{ a[$1]; next }} ($1 in a && $2 in a && $6 >= {params.min_lod} && $9 >= {params.min_l} && $9 <= {params.max_l}) {{print $3"\\t"$4-1"\\t"$5}}' {input.inds} - | \
        bedtools genomecov -i stdin -g {input.genome} -bga | awk '$1 == "{wildcards.chrom}"' | gzip > {output}
        """


rule mask_ibd:
    input:
        RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_coverage.tsv.gz"
    output:
        bed=RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_mask.bed",
        pdf=RESULTS_DIR + "/masking/plots/{chrom}." + PREFIX + ".ibd_mask.pdf",
    params:
        ibd_trim=MASK_TRIM,
        ibd_sd=MASK_SD,
    shell:
        """
        mkdir -p {RESULTS_DIR}/masking/tables {RESULTS_DIR}/masking/plots
        Rscript workflow/scripts/r/mask_ibd_regions.R -i {input} -b {output.bed} -p {output.pdf} --trim {params.ibd_trim} --sd {params.ibd_sd}
        """
