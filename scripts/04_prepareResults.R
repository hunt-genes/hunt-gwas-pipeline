## library("plotrix")
library("data.table")
## library("RColorBrewer")


## summary.infile = snakemake@input[["summary"]]
imputed.infile = snakemake@input[["imputed"]]

final.results.outfile = snakemake@output[["final_results"]]
candidate.regions.outfile = snakemake@output[["candidate_regions"]]


# Settings for plot
plotLine <- T
gwthreshold <- 5E-8
yLine <- -log10(gwthreshold)
colLine <- "red"

chrcol = snakemake@config[["dosage_chr_col"]]
idcol = snakemake@config[["dosage_snp_id_col"]]
poscol = snakemake@config[["dosage_pos_col"]]
refcol = snakemake@config[["dosage_ref_allele_col"]]
altcol = snakemake@config[["dosage_alt_allele_col"]]

# Get some basic information about phenotypes to determine title and MAF threshold
# read phenotype summary
## phenoSummary <- read.table(summary.infile,sep="\t",comment.char="",
## 	header=T,check.names=F)
## traitname <- phenoSummary[1,1]
traitname <- snakemake@wildcards[["pheno"]]

# Read input
file_imp <- imputed.infile

analysis = snakemake@params[["analysis"]]

if(analysis == "BOLTLMM"){
	pcol <- "P_BOLT_LMM_INF"
	ycol <- "log10P_GC"
	chicol <- "CHISQ_BOLT_LMM_INF"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "BP"
	idcol <- "ID"
	refcol <- ""
	altcol <- ""
	accol <- ""
	maccol <- ""
	afcol <- "A1FREQ"
	af0col <- ""
	af1col <- ""
} else if (analysis == "SAIGE"){
	pcol <- "p.value"
	ycol <- "log10P_GC"
	chicol <- "CHISQ"
	betacol <- "BETA"
	secol <- "SE"
	## chrcol <- "CHR"
	## poscol <- "POS"
	## idcol <- "ID"
	## refcol <- ""
	## altcol <- ""
	accol <- ""
	maccol <- ""
	afcol <- "AF"
	af0col <- ""
	af1col <- ""
} else {
	stop("Unknown Analysis")
}


maf_threshold <- snakemake@params[["maf_threshold"]]


# convert -log10(P) values to as.character(P)
log10toP <- function(log10P){
		log10P <- abs(as.numeric(log10P))
		if(is.na(log10P)) return(NA)
		if(log10P > 300){
			part1 <- log10P%/%100*100
			part2 <- log10P-part1
			if(part2 != 0){
				P <- format(signif(10^-part2,3), scientific = T)
				P <- paste(as.numeric(gsub("e-.+","",P)),"e-",as.numeric(gsub(".+-","",P),sep="")+part1,sep="")
			} else {
				P <- paste("1e-",part1,sep="")
			}
		} else {
			P <- signif(10^-log10P,3)
		}
	return(as.character(P))
}

# Read input file; add -log10 P values
resIn <- fread(file_imp) # read.table(file_imp, header=T)

# Remove invalid results (SE ==0)
if(length(which(resIn[[secol]] == 0)) > 0){
    resIn <- resIn[which(resIn[[secol]] > 0),]
}

resIn$log10P <- -log10(resIn[[pcol]])
# If p values are too small for R (=0). Get log10 p-values from CHISQ or SCORE statistics
noP <- which(resIn[[pcol]] == 0)
# log10(P) from CHISQ
if(length(noP)>0 & chicol %in% colnames(resIn)){
    resIn$log10P[noP] <- signif(-pchisq(resIn[[chicol]][noP],df=1,lower.tail=F,log.p = T)/log(10),5)
    resIn[[pcol]][noP] <- sapply(resIn$log10P[noP],log10toP)
# log10(P) from SCORE statistics
} else if (length(noP)>0 & !chicol %in% colnames(resIn)){
	resIn$log10P[noP] <- rep(-log10(as.numeric(c(unlist(format(.Machine)))[3])),length(noP))
	#sapply(resIn$Tstat[noP],Score2log10P)
    resIn[[pcol]][noP] <- sapply(resIn$log10P[noP],log10toP)
}

# Add overall MAF

resIn$MAF <- ifelse(resIn[[afcol]] <= 0.5,resIn[[afcol]],1-resIn[[afcol]])
resIn <- resIn[which(resIn$MAF >= maf_threshold),]

# ADD CHISQ, calculate from p value;  could be faster with a parallel solution
if(!chicol %in% colnames(resIn)){
	resIn$CHISQ <- qchisq(-resIn$log10P*log(10), df=1, lower.tail = F,log.p=T)
}

# Minimal MAF
minMAF <- min(resIn$MAF,na.rm=T)

# Lambda estimation using SNPs with MAF > 0.01 and GC correction if necessary
denom<-qchisq(0.5, df=1)
Lambda <- qchisq(10^-median(resIn$log10P[which(resIn$MAF >= 0,01)],na.rm=T), df=1, lower.tail = F) / denom
print(Lambda)
if(Lambda > 1){
	resIn$log10P_GC <- -pchisq(resIn[[chicol]], df=1, lower.tail = F, log.p = T)/log(10)
	resIn$SEBETA_GC <- signif(resIn[[secol]]*sqrt(Lambda),5)
} else {
	resIn$log10P_GC <- resIn$log10P
	resIn$SEBETA_GC <- resIn[[secol]]
}

print(head(resIn))
print(dim(resIn))
print(head(resIn))
## Clean GWAS output

if(analysis == "BOLTLMM"){
	resOut <- data.frame(
		'CHROM'=resIn[[chrcol]],
		'POS'=resIn[[poscol]],
		'ID'=resIn$SNP,
		'REF'=resIn$ALLELE0,
		'ALT'=resIn$ALLELE1,
#		'AC'=NA,
#		'MAC'=NA,
		'MAF'=resIn$MAF,
		'AF'=resIn$A1FREQ,
#		'AF.CTRL'=NA,
#		'AF.CASE'=NA,
		'CHISQ'=resIn[[chicol]],
		'BETA'=resIn[[betacol]], #Approxmation
		'SEBETA_GC'=resIn$SEBETA_GC, #Approxmation
		'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		'LOG10P_GC'=resIn[[ycol]])
} else if (analysis == "SAIGE"){
	resOut <- data.frame(
		'CHROM'=resIn[[chrcol]],
		'POS'=resIn[[poscol]],
		'ID'=resIn[[idcol]],
		'REF'=resIn[[refcol]],
		'ALT'=resIn[[altcol]],
		'AC'=resIn$AC,
		'MAC'=ifelse(resIn$AC > resIn$N,resIn$N*2-resIn$AC,resIn$AC),
		'MAF'=resIn$MAF,
		'AF'=resIn[[afcol]],
		'CHISQ'=resIn[[chicol]],
		'BETA'=resIn[[betacol]],
		'SEBETA_GC'=resIn$SEBETA_GC,
		'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		'LOG10P_GC'=resIn[[ycol]])
}
write.table(resOut, final.results.outfile,
	sep="\t",col.names=T,row.names=F,quote=F)


## "# Determine candidate regions if not provided"
candidateHits <- which(resIn$log10P_GC > -log10(gwthreshold))
## "# BOLTLMM"
if(analysis =="BOLTLMM" & length(candidateHits) > 0){
	tophits <- resIn[candidateHits,]

	# merge candidate regions if they are within 0.025 cM
	tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
	x <- as.numeric(tophits$GENPOS)
	y <- tophits$numCHR
	start = c(1, which(diff(y) != 0 | diff(x) <= -0.025 | diff(x) >= 0.025) + 1)
	end = c(start - 1, length(x))
	candidateRegions <- data.frame(
		'CHROM'=tophits$CHR[start],
		'START'=tophits$BP[start] - 50000,
		'END'=tophits$BP[end] + 50000,
		'COL'="blue",
		'LEGENDTEXT'="Candidate Regions [Top Hits +-(0.025cM + 50kb)]"
		)
# SAIGE
} else if (analysis =="SAIGE" & length(candidateHits) > 0){
	tophits <- resIn[candidateHits,]
	tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
	x <- as.numeric(tophits$POS)
	y <- tophits$numCHR

	start = c(1, which(diff(y) != 0 | diff(x) <= -500000 | diff(x) >= 500000) + 1)
	end = c(start - 1, length(x))
	candidateRegions <- data.frame(
		'CHROM'=tophits$CHR[start],
		'START'=tophits$POS[start] - 500000,
		'END'=tophits$POS[end] + 500000,
		'COL'="blue",
		'LEGENDTEXT'="Candidate Regions [+-500kb]"
		)
} else {
	candidateRegions <- data.frame(
		'CHROM'=character(0),
		'START'=numeric(0),
		'END'=numeric(0),
		'COL'=character(0),
		'LEGENDTEXT'=character(0)
)
}

write.table(candidateRegions[c("CHROM","START","END","COL","LEGENDTEXT")],
	candidate.regions.outfile,sep="\t",col.names=T,row.names=F,quote=F)
