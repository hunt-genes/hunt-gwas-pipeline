
from utils.helpers import title

rule qqplot:
    input:
        "{prefix}/data/imputed/{pheno}.final_results.txt"
    output:
        "{prefix}/data/qqplot/{pheno}.QQ_Plot.png"
    params:
        prefix = "{prefix}/data/qqplot/{pheno}",
        title = lambda w: title(pk, w)
    benchmark:
        "{prefix}/benchmarks/qqplot/{pheno}.benchmark"
    conda:
        "../../envs/qq.yaml"
    shell:
	    "Rscript --vanilla ./scripts/QQplot.r "
	    "--input {input[0]} "
	    "--prefix {params.prefix} "
	    "--maf MAF "
	    "--pvalue LOG10P_GC "
	    "--log10p T "
	    "--maintitle '{params.title}' "
