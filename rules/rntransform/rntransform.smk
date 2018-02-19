

# def rntransform_formula(w):

#     right_covar_type = (pk.covar_type.isin(["addQcovar", "addBcovar"]))
#     covariates = list(pk.loc[(pk.name == w.pheno) & right_covar_type].covar_name)

#     return w.pheno + " ~ " + " + ".join(covariates)


if pheno_continuous:
    rule rntransform:
        input:
            config["phenotype_file"]
        output:
            "{prefix}/data/rntransform/{pheno}_rntransformed.txt"
        conda:
            "../../envs/rntransform.yaml"
        script:
            "../../scripts/rntransform.R"
