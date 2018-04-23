
def phenofile(w, pheno_file, pheno_contains_continuous, extract_unrelated, covar_type):

    if pheno_contains_continuous:
        f = "{prefix}/data/rntransform/{pheno}_rntransformed.txt"
    elif extract_unrelated:
        f = "{prefix}/data/extract_unrelated/{pheno}_unrelated.txt"
    else:
        return pheno_file

    return expand(f, **w)


rule saige_fit_null_logistic_mixed_model:
    input:
        pheno = lambda w: phenofile(w, config["phenotype_file"],
                                    pheno_contains_continuous, remove_related, pheno_trait_type[w.pheno]),
        plink_binary_pruned = config["plink_binary_prefix_pruned"] + ".bed"
    output:
        glmm_model_information = "{prefix}/data/saige_step1/{pheno}.rda",
        variance_ratio = "{prefix}/data/saige_step1/{pheno}.varianceRatio.txt",
        spagmmat = f"{{prefix}}/data/saige_step1/{{pheno}}_{num_markers}markers.SAIGE.results.txt"
    params:
        plink_prefix = config["plink_binary_prefix_pruned"],
        trait_type = lambda w: pheno_trait_type[w.pheno],
        covariates = lambda w: list(pk.loc[pk.name == w.pheno].covar_name),
        output_prefix = "{prefix}/data/saige_step1/{pheno}"
    threads:
        48
    benchmark:
        "{prefix}/benchmarks/saige_step1/{pheno}.benchmark"
    singularity:
        "shub://endrebak/singularity_recipes:saige_0_26"
    script:
        "../../scripts/03.1_SAIGE_Step1.R"
