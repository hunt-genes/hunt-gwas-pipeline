from utils.helpers import title



rule manhattan_plot:
    input:
        final_results = "{prefix}/data/imputed/{pheno}.final_results.txt",
        candidate_regions = "{prefix}/data/imputed/{pheno}.imputed.candidate_regions"
    output:
        graph = "{prefix}/data/manhattan/{pheno}.imputed_Manhattan.png"
    params:
        prefix = "{prefix}/data/manhattan/{pheno}.imputed",
        title = lambda w: title(pk, w),
        cutoff = config["logp_cutoff"]
    benchmark:
        "{prefix}/benchmarks/manhattan_plot/{pheno}.benchmark"
    conda:
        "../../envs/manhattan.yaml"
    shell:
        "Rscript --vanilla scripts/ManhattanPlot.R "
        "--input {input.final_results} "
        "--prefix {params.prefix} "
        "--hitregion {input.candidate_regions} "
        "--chr CHROM "
        "--pos POS "
        "--pvalue LOG10P_GC "
        "--cutoff {params.cutoff} "
        "--log10p T "
        "--maintitle '{params.title}'"
