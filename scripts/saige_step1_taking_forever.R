
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('test/data/extract_unrelated/VT_BO_unrelated.txt', "pheno" = 'test/data/extract_unrelated/VT_BO_unrelated.txt'),
    output = list('test/data/saige_step1/VT_BO.rds'),
    params = list('/mnt/scratch/hunt24chunks/genotyped', 'binary', c('BirthYear', 'PC1', 'PC2', 'PC3', 'PC4'), NULL, 'test/data/saige_step1/VT_BO', "plink_prefix" = '/mnt/scratch/hunt24chunks/genotyped', "trait_type" = 'binary', "covariates" = c('BirthYear', 'PC1', 'PC2', 'PC3', 'PC4'), "center_variables" = NULL, "output_prefix" = 'test/data/saige_step1/VT_BO'),
    wildcards = list('test', 'VT_BO', "prefix" = 'test', "pheno" = 'VT_BO'),
    threads = 48,
    log = list(),
    resources = list(),
    config = list("prefix" = 'test', "phenotype_file" = '/mnt/cargo/benb_in/VT_BO_for_Lars.txt', "key_file" = '/mnt/cargo/benb_in/VT_BO_key_file.txt', "id_column" = 'IID', "use_unrelated" = TRUE, "plink_binary_prefix" = '/mnt/scratch/hunt24chunks/genotyped', "king_kinship" = 'example_data/king.kin0', "prioritize" = TRUE, "flip_priority" = FALSE, "unrelated_threshold" = 0.08838, "number_of_runs" = 2, "seed" = 1260646394, "number_of_markers" = 30, "skip_model_fitting" = FALSE),
    rule = 'saige_step1'
)
######## Original script #########
options(stringsAsFactors=F)

library("SAIGE")

## option_list <- list(
##   make_option("--plinkFile", type="character",default="",
##     help="path to plink file to be used for the kinship matrix"),
##   make_option("--phenoFile", type="character", default="",
##     help="path to the phenotype file, a column 'IID' is required"),
##   make_option("--phenoCol", type="character", default="",
##     help="path to the phenotype file, a column 'IID' is required"),
##   make_option("--covarColList", type="character", default="",
##     help="list of covariates (comma separated)"),
##   make_option("--sampleIDColinphenoFile", type="character", default="",
##     help="Column name of the IDs in the phenotype file"),
##   make_option("--centerVariables", type="character", default="NULL",
##     help="Covariates that should be centered (comma separated) [default='NULL']"),
##   make_option("--skipModelFitting", type="logical", default=FALSE,
##     help="skip model fitting, [default='FALSE']"),
##   make_option("--traitType", type="character", default="binary",
##     help="binary/quantitative [default=binary]"),
##   make_option("--outputPrefix", type="character", default="~/",
##     help="path to the output files [default='~/']"),
##   make_option("--numMarkers", type="integer", default=30,
##     help="An integer greater than 0 Number of markers to be used for estimating the variance ratio [default=30]"),
##   make_option("--nThreads", type="integer", default=12,
##     help="Number of threads")
## )


## parser <- OptionParser(usage="%prog [options]", option_list=option_list)

## args <- parse_args(parser, positional_arguments = 0)
## opt <- args$options
## print(opt)

phenoFile = snakemake@input[["pheno"]]

plinkFile = snakemake@params[["plink_prefix"]]
traitType = snakemake@params[["trait_type"]]
outputPrefix = snakemake@params[["output_prefix"]]

phenoCol = snakemake@wildcards[["pheno"]]
nThreads = snakemake@threads
print(nThreads)

sampleIDColinphenoFile = snakemake@config[["id_column"]]
numMarkers = snakemake@config[["id_column"]]
skipModelFitting = snakemake@config[["skip_model_fitting"]]
#run step 1 to fit the null glmm model and estimate the constant of variance ratio

covars = snakemake@params[["covariates"]]

## if(opt$centerVariables == "NULL"){

centerVariables = snakemake@params[["center_variables"]]

fitNULLGLMM(plinkFile=plinkFile,
            phenoFile = phenoFile,
            phenoCol = phenoCol,
            traitType = traitType,
            invNormalize = FALSE,
            covarColList = covars,
            qCovarCol = NULL, # not implemented yet
            sampleIDColinphenoFile = sampleIDColinphenoFile,
            centerVariables = centerVariables,
            tol=0.02,
            maxiter=20,
            tolPCG=1e-5,
            maxiterPCG=500,
            nThreads = nThreads,
            Cutoff = 2,
            numMarkers = numMarkers,
            skipModelFitting = skipModelFitting,
            outputPrefix = outputPrefix)

## } else {
## 	fitNULLGLMM(plinkFile=opt$plinkFile,
## 		phenoFile = opt$phenoFile,
## 		phenoCol = opt$phenoCol,
## 		traitType = opt$traitType,
## 		invNormalize = FALSE,
## 		covarColList = covars,
## 		qCovarCol = NULL,
## 		sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
## 		centerVariables = strsplit(opt$centerVariables,",")[[1]],
## 		tol=0.02,
## 		maxiter=20,
## 		tolPCG=1e-5,
## 		maxiterPCG=500,
## 		nThreads = opt$nThreads,
## 		Cutoff = 2,
## 		numMarkers = opt$numMarkers,
## 		skipModelFitting = opt$skipModelFitting,
## 		outputPrefix = opt$outputPrefix)
## }
