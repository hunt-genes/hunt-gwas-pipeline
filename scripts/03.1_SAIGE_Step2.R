library(SAIGE)

# check if openblas installed R
## require(inline)

## openblas.set.num.threads <- cfunction(signature(ipt="integer"),
##      body = 'openblas_set_num_threads(*ipt);',
##      otherdefs = c ('extern void openblas_set_num_threads(int);'),
##      libargs = c ('-L/opt/openblas/lib -lopenblas'),
##      language = "C",
##      convention = ".C"
## )

## openblas.set.num.threads(snakemake@threads)

dosageFile = snakemake@input[["dosage"]]

SPAGMMATtest(dosageFile = dosageFile,
             sampleFile = snakemake@config[["samples_in_dosage"]],
             dosageFileNcolSkip=snakemake@config[["number_dosage_cols_to_skip"]],
             dosageFileNrowSkip=snakemake@config[["number_dosage_rows_to_skip"]],
             minMAC = snakemake@config[["min_minor_allele_count"]],
             minMAF = snakemake@config[["min_minor_allele_frequency"]],
             ## phenoFile = snakemake@input[["pheno"]],
             ## phenoCol = snakemake@wildcards[["pheno"]],
             ## covarColList = snakemake@params[["covars"]],
             ## sampleIDColinphenoFile = snakemake@config[["id_column"]],
             ## centerVariables = snakemake@config[["center_variables"]],
             GMMATmodelFile = snakemake@input[["rda"]],
             varianceRatioFile = snakemake@input[["variance_ratio_file"]],
             SAIGEOutputFile = snakemake@output[[1]],
             dosageFilecolnamesSkip = snakemake@config[["dosage_file_colnames_skip"]])


## SPAGMMATtest(
##   vcfFile = opt$vcfFile,
##   vcfFileIndex = paste0(opt$vcfFile,".tbi"),
##   vcfField = opt$vcfField,
##   dosageFile=opt$dosageFile,
##   dosageFileNrowSkip=opt$dosageFileNrowSkip,
##   dosageFileNcolSkip=opt$dosageFileNcolSkip,
##   sampleFile=opt$sampleFile,
##   idstoIncludeFile = opt$sampleFile,
##   minMAC = opt$minMAC,
##   minMAF = opt$minMAF,
##   minInfo = opt$minInfo,
##   phenoFile=opt$phenoFile,
##   phenoCol=opt$phenoCol,
##   traitType=opt$traitType,
##   covarColList=covars,
##   sampleIDColinphenoFile=opt$sampleIDColinphenoFile,
##   GMMATmodelFile=opt$GMMATmodelFile,
##   varianceRatioFile=opt$varianceRatioFile,
##   SAIGEOutputFile=opt$SAIGEOutputFile,
##   IsOutputAFinCaseCtrl = IsOutputAFinCaseCtrl)
