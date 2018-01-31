# Rscript written by Lars Fritsche to select relatively unrelated samples
# relatively unrelated = not related to 2nd degree or closer
# April 2017
# Copyright Lars Fritsche

options(stringsAsFactors=F)

## source("config.txt")

pickUnrel <- function(unrelated.threshold, number.of.runs, seed, file.kin0,iids,ex.iids=character(0)){
	options(stringsAsFactors=F)
	iids <- as.character(iids)
	ex.iids <- as.character(ex.iids)
	tmp0 <- read.table(file.kin0, header = TRUE)
	# fastindep doesn't like negative values
	tmp0$Kinship[which(tmp0$Kinship < 0)] <- 0
	tmp0 <- tmp0[which(tmp0$Kinship > 0.08839),]

	# find all samples that are related to samples of the exclusion list
	ex.iids2 <- unique(c(tmp0$ID2[which(tmp0$ID1 %in% ex.iids)],tmp0$ID1[which(tmp0$ID2 %in% ex.iids)],ex.iids))

	# reduce kinship matrix to samples that are related to "iids" but not related to excluded samples
	tmp0 <- tmp0[which(!tmp0$ID1 %in% ex.iids2 & !tmp0$ID2 %in% ex.iids2
						& tmp0$ID1 %in% iids & tmp0$ID2 %in% iids),]

	# determine remaining singletons in data (i.e. samples that have no related individual in the data set)
	singletons <- iids[which(!iids %in% unique(c(tmp0$ID1,tmp0$ID2)) & !iids %in% ex.iids2)]

	if(dim(tmp0)[1] == 0){
		return(singletons)
		next
	}

	uid01 <- as.character(unique(tmp0$ID1))
	uid02 <- as.character(unique(tmp0$ID2))
	ids <- as.character(unique(c(uid01, uid02)))
	nsamp <- length(ids)

    if(nsamp > 30000){
    	# create 4 matrices
 		# larger matrices cause some troubles
		file.temp <- tempfile(fileext=".dat")
    	parts <- split(ids, cut(1:length(ids), 4))
    	for(p in 1:4){
	    	kingMat <- matrix(0,nrow = length(parts[[p]]), ncol = nsamp)
	    	tmp0p <- tmp0[which(tmp0$ID1 %in% parts[[p]] | tmp0$ID2 %in% parts[[p]]),]
			colnames(kingMat) <- ids
			rownames(kingMat) <- parts[[p]]
			nr <- dim(tmp0p)[1]
			id1 <- as.character(tmp0p$ID1)
			id2 <- as.character(tmp0p$ID2)
			for (r in 1:nr) {
				if(id1[r] %in% parts[[p]]) kingMat[id1[r],id2[r]] <- tmp0p$Kinship[r]
				if(id2[r] %in% parts[[p]]) kingMat[id2[r],id1[r]] <- tmp0p$Kinship[r]
			}
			if(p == 1){
				colnames(kingMat)[1] <- paste0("\t",ids[1])
				write.table(kingMat,file.temp,sep="\t",row.names=T,col.names=T,quote=F)
			} else {
				write.table(kingMat,file.temp,sep="\t",row.names=T,col.names=F,quote=F,append=T)
			}
			rm("kingMat")
			gc()
		}
	} else {
		# standard matrix for smaller datasets
		kingMat <- matrix(0,nrow = nsamp, ncol = nsamp)
		nr <- dim(tmp0)[1]
		id1 <- as.character(tmp0$ID1)
		id2 <- as.character(tmp0$ID2)
		for (r in 1:nr) {
			x <- which(ids == id1[r])
			y <- which(ids == id2[r])
			kingMat[x,y] <- tmp0$Kinship[r]
			kingMat[y,x] <- tmp0$Kinship[r]
		}
		# for some reason fastindep requires a tab in the first column
		colnames(kingMat) <- c(paste0("\t",ids[1]),ids[2:length(ids)])
		rownames(kingMat) <- ids
		file.temp <- tempfile(fileext=".dat")
		write.table(kingMat,file.temp,sep="\t",row.names=T,col.names=T,quote=F)
	}

	# copy fastindep and parameters to tempdir; run fastindep
	## file.copy(fastindep,paste0(tempdir(),"/fastIndep"))
	## file.copy(parameterfile,paste0(tempdir(),"/parameter.txt"))

	## system(paste0("cd ",tempdir(),"; ",tempdir(),"/fastIndep ",gsub(".+/([^/]+)$","\\1",file.temp)," parameter.txt > /dev/null"))

  temp.outfile = paste0(tempdir(),"/outfile.txt")
  command = paste("fastindep", file.temp,unrelated.threshold, number.of.runs, seed, ">", temp.outfile)
  write(command, stderr())
  system(command)
  write("Writing massive amounts of data. Please be patient.", stderr())
	# extract maximal set of unrelated individuals; from greedy algorithm
	fastOut <- readLines(temp.outfile)
	set1 <- which(grepl("Set Size",fastOut))[1]
	fastUnrel <- strsplit(fastOut[set1]," +")[[1]]
	fastUnrel <- fastUnrel[6:length(fastUnrel)]
	# return list of unrelated samples (maximal set + singletons)
    return(c(singletons,fastUnrel))
}


## option_list <- list(
##   make_option("--file.kin0", type="character",default="./data/king.kin0",
##     help="Full path to kinship file from KING [default=''./data/king.kin0]"),
##   make_option("--file.pheno", type="character", default="",
##     help="Full path to pheno file (1st column 'IID', 2nd column 'Priority')"),
##   make_option("--col.priority", type="character", default="Priority",
##     help="Column name with priority information [default='Priority'])"),
##   make_option("--file.bed", type="character", default="/mnt/scratch/hunt24chunks/genotyped",
##     help="Prefix to binary PLINK files [default='/mnt/scratch/hunt24chunks/genotyped']"),
##   make_option("--priorization", type="logical", default=FALSE,
##     help="Prioritize cases; (larger numbers have priority) [default=FALSE]"),
##   make_option("--flip.priority", type="logical", default=FALSE,
##     help="Prioritize smaller numbers [default=FALSE]"),
##   make_option("--file.out", type="character", default="",
##     help="Full path to pheno file with unrelated samples")
## )



file.pheno = snakemake@input[["pheno"]]
file.kin0 = snakemake@input[["king_kinship"]]
## paramsfile = snakemake@input[["fastindep_params"]]

col.priority = snakemake@wildcards[["pheno"]]

plink_prefix = snakemake@params[["plink_binary_prefix"]]
prioritization = snakemake@params[["prioritization"]]
flip.priority = snakemake@params[["flip_priority"]]
unrelated.threshold = snakemake@params[["unrelated_threshold"]]
number.of.runs = snakemake@params[["number_of_runs"]]
seed = snakemake@params[["seed"]]

fam <- read.table(paste0(plink_prefix, ".fam"), sep=" ", header=F)
phenoIn <- read.table(file.pheno,sep="\t",header=T,quote="\"",comment.char="")
phenoIn <- phenoIn[which(phenoIn$IID %in% fam[[2]]),]


file.out = snakemake@output[[1]]


if(prioritization){
	selected.samples <- character(0)
	for(p in sort(unique(phenoIn[[col.priority]]),decreasing=!flip.priority)){
		selected.samples <- c(selected.samples,
								pickUnrel(unrelated.threshold, number.of.runs, seed, file.kin0,phenoIn$IID[which(phenoIn[[col.priority]] == p)],selected.samples))
	}
} else {
	selected.samples <- pickUnrel(unrelated.threshold, number.of.runs, seed, file.kin0, phenoIn$IID)
}


write.table(phenoIn[which(phenoIn$IID %in% selected.samples),],file.out,sep="\t",col.names=T,row.names=F,quote=F)
