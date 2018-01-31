
rule collect_chunks:
    input:
        expand("{prefix}/data/saige_step2/{{pheno}}_{chunk}.saige", chunk=dosages.keys(), prefix=prefix)
    output:
        f"{prefix}/data/collect_imputed_files/{{pheno}}.imputed.results"
    benchmark:
        f"{prefix}/benchmarks/prepare_for_plotting/{{pheno}}.benchmark"
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} " # only keep header from first file
        "| tr ' ' '\t' > {output[0]}" # replace spaces with tabs (if SAGE)
