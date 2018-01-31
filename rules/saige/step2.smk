rule saige_test_associate_snps_and_phenotype:
    input:
        rda = "{prefix}/data/saige_step1/{pheno}.rda",
        dosage = lambda w: dosages[w.chunk],
        pheno = config["phenotype_file"] if not remove_related else "{prefix}/data/extract_unrelated/{pheno}_unrelated.txt",
        variance_ratio_file = "{prefix}/data/saige_step1/{pheno}.varianceRatio.txt",
        spagmmat = f"{{prefix}}/data/saige_step1/{{pheno}}_{num_markers}markers.SAIGE.results.txt"
    output:
        "{prefix}/data/saige_step2/{pheno}_{chunk}.saige"
    params:
        trait_type = lambda w: pk.loc[pk.name == w.pheno].type.iloc[0],
        covars = lambda w: list(pk.loc[pk.name == w.pheno, "covar_name"].dropna())
    threads:
        1
    benchmark:
        "{prefix}/benchmarks/saige_step2/{pheno}_{chunk}.benchmark"
    singularity:
        "shub://endrebak/singularity_recipes:saige_0_24"
    script:
        "../../scripts//03.1_SAIGE_Step2.R"
