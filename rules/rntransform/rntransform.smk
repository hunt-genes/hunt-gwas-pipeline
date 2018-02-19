

if pheno_contains_continuous:

    rule rntransform:
        input:
            config["phenotype_file"] if not remove_related else "{prefix}/data/extract_unrelated/{pheno}_unrelated.txt"
        output:
            "{prefix}/data/rntransform/{pheno}_rntransformed.txt"
        params:
            continuous_phenotypes = [k for k, v in pheno_trait_type.items() if v == "quantitative"]
        conda:
            "../../envs/rntransform.yaml"
        script:
            "../../scripts/rntransform.R"
