rule saige_associate_snps_and_phenotype:
    input:
        pheno = lambda w: phenofile(w, config["phenotype_file"],
                                    pheno_contains_continuous, remove_related, pheno_trait_type[w.pheno]),
        rda = "{prefix}/data/saige_step1/{pheno}.rda",
        dosage = lambda w: dosages[w.chunk],
        variance_ratio_file = "{prefix}/data/saige_step1/{pheno}.varianceRatio.txt",
        spagmmat = f"{{prefix}}/data/saige_step1/{{pheno}}_{num_markers}markers.SAIGE.results.txt"
    output:
        "{prefix}/data/saige_step2/{pheno}_{chunk}.saige"
    threads:
        24
    benchmark:
        "{prefix}/benchmarks/saige_step2/{pheno}_{chunk}.benchmark"
    singularity:
        "shub://endrebak/singularity_recipes:saige_0_26"
    script:
        "../../scripts//03.1_SAIGE_Step2.R"
