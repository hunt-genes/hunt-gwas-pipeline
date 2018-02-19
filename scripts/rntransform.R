library(GenABEL)

phenotype_file = snakemake@input[[1]]
phenotype = snakemake@wildcards[["pheno"]]


pheno_df = read.table(phenotype_file, sep="\t", header=TRUE, check.names=FALSE, comment.char="")

for (phenotype in snakemake@params[["continuous_phenotypes"]]) {

  transformed = rntransform(pheno_df[[phenotype]], pheno_df, family=gaussian)

  pheno_df[[phenotype]] = transformed

}


write.table(pheno_df, snakemake@output[[1]], row.names=FALSE)
