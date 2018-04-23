## CHR     SNP     CM      POS     EFFECT_ALLELE   ALT_ALLELE      AC      AF      N       BETA    SE      Tstat   p.value varT    varTstar
## 1:55389_T/C     1       55389   C       T       0       76.3040006353986        0.000695088185353798    54888   -0.214383996486021      0.180851760172143       -0.75036613989207       0.235854262201548  0.400688807472199        0.400231742821897
## 1:61743_G/C     1       61743   C       G       0       215.464002979919        0.00196276055768036     54888   -0.173705474638331      0.106399830612812       -1.04530799248364       0.102558865355471  0.409961810832324        0.409494168491906
## 1:63210_G/T     1       63210   T       G       0       149.058997494518        0.00135784686538513     54888   0.298280764067639       0.119013834953428       1.72484966910774        0.012201252558086  0.473638399400642        0.473098121346065
## 1:77831_C/A     1       77831   A       C       0.001   251.407002032851        0.00229018184332505     54888   -0.15508157950062       0.10346958310895        -0.913579051303828      0.133922092863902  0.371532861313694        0.371109054773143
## 1:79629_A/C     1       79629   C       A       0.001   197.618001500028        0.00180019313420081     54888   -0.109660278605297      0.122891126347682       -0.516529397728492      0.37221241692585   0.335067609737606        0.334685399011934
## 1:79665_C/A     1       79665   A       C       0       7.76500005857088        7.07349517068474e-05    54888   0.671583580310185       0.513422483122017       0.914280742165422       0.190855491233176  0.488549722626571        0.487992435265498
## 1:79898_A/G     1       79898   G       A       0.001   197.626001413213        0.00180026600908407     54888   -0.109704526232827      0.122875338929221       -0.516860147168945      0.371958080726129  0.335140149351771        0.334757855880356
## 1:84295_G/A     1       84295   A       G       0.001   245.862003767281        0.00223966990751422     54888   -0.157256682848316      0.10284461236136        -0.948200785827074      0.126246919386325  0.384543515259425        0.384104867500739
## 1:85011_C/T     1       85011   T       C       0       107.157999057323        0.0009761514270635      54888   0.169775096894572       0.138794225957926       0.851371074448565       0.221248725196917  0.48443175747252 0.483879167463277

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
Lambda <- qchisq(10^-median(resIn$log10P[which(resIn$MAF >= 0.01)],na.rm=T), df=1, lower.tail = F) / denom
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
