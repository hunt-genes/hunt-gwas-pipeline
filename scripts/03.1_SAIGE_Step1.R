options(stringsAsFactors=F)

library("SAIGE")

phenoFile = snakemake@input[["pheno"]]

plinkFile = snakemake@params[["plink_prefix"]]
traitType = snakemake@params[["trait_type"]]
outputPrefix = snakemake@params[["output_prefix"]]

phenoCol = snakemake@wildcards[["pheno"]]
nThreads = snakemake@threads

sampleIDColinphenoFile = snakemake@config[["id_column"]]
numMarkers = snakemake@config[["number_of_markers"]]
skipModelFitting = snakemake@config[["skip_model_fitting"]]

covars = snakemake@params[["covariates"]]

centerVariables = snakemake@params[["center_variables"]]

#run step 1 to fit the null glmm model and estimate the constant of variance ratio
fitNULLGLMM(plinkFile = plinkFile,
            phenoFile = phenoFile,
            phenoCol = phenoCol,
            traitType = traitType,
            invNormalize = FALSE,
            covarColList = covars,
            qCovarCol = NULL, # not implemented yet
            sampleIDColinphenoFile = sampleIDColinphenoFile,
            tol=0.02,
            maxiter=20,
            tolPCG=1e-5,
            maxiterPCG=500,
            nThreads = nThreads,
            Cutoff = 2,
            numMarkers = numMarkers,
            skipModelFitting = skipModelFitting,
            outputPrefix = outputPrefix)
