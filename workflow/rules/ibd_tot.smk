IBD_EXCL_MASKS = expand(
    RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_excl_mask.bed",
    chrom=CHROMS,
)

IBD_TOT_OUTPUTS = expand(
    RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_tot.tsv.gz",
    chrom=CHROMS,
)


rule get_excl_mask:
    input:
        mask=RESULTS_DIR + "/masking/tables/{chrom}." + PREFIX + ".ibd_mask.bed"
    output:
        temp(RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_excl_mask.bed")
    shell:
        """
        mkdir -p {RESULTS_DIR}/ibd_tot/tables
        awk '$4 != "0"' {input.mask} > {output}
        """


rule get_ibd:
    input:
        ibd=IBD_INPUT,
        inds=INDIVIDUALS,
        mask=RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_excl_mask.bed",
    output:
        RESULTS_DIR + "/ibd_tot/tables/{chrom}." + PREFIX + ".ibd_tot.tsv.gz"
    params:
        min_l=AGG_MIN_L,
        max_l=AGG_MAX_L,
        min_lod=AGG_MIN_LOD,
    threads:
        2
    shell:
        """
        set -euo pipefail
        mkdir -p {RESULTS_DIR}/ibd_tot/tables
        gzip -cd {input.ibd} | tail -n+2 | \
        awk 'FNR==NR {{ a[$1]; next }} ($1 in a && $2 in a && $6 >= {params.min_lod} && $9 >= {params.min_l} && $9 <= {params.max_l}) {{print $3"\\t"$4"\\t"$5"\\t"$1"\\t"$2"\\t"$9}}' {input.inds} - | \
        bedtools intersect -a stdin -v -b {input.mask} | sort -S10% -k4,4 -k5,5 --parallel={threads} | datamash -g4,5 first 1 sum 6 | gzip > {output}
        """
