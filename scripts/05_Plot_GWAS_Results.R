
## CHROM   POS     ID      REF     ALT     AC      MAC     MAF     AF      CHISQ   BETA    SEBETA_GC       PVALUE_GC       LOG10P_GC
## 1:55389_T/C     C       1       T       0       45.9770003277808        45.9770003277808        0.00069581996984958     0.00069581996984958     0.620256737964737       -0.938151262862578      1.1914  0.431       0.365571203695509
## 1:61743_G/C     C       1       G       0       131.939001867315        131.939001867315        0.00199677646751188     0.00199677646751188     0.0333648146234227      0.125904675782828       0.6894  0.855       0.0680011217025361
## 1:63210_G/T     T       1       G       0       90.7519983953098        90.7519983953098        0.00137344873169244     0.00137344873169244     0.0592093303810495      -0.188291258354402      0.77394 0.808       0.092723023454143
## 1:77831_C/A     A       1       C       0.001   152.892001248896        152.892001248896        0.00231388100443272     0.00231388100443272     0.000987121332283073    -0.0210775720227709     0.67098 0.975       0.011023982823315
## 1:79629_A/C     C       1       A       0.001   124.60100105172 124.60100105172 0.00188572251727889     0.00188572251727889     2.07039998593853        1.09871702936133        0.76372 0.15    0.823381743201756
## 1:79898_A/G     G       1       A       0.001   124.606001027161        124.606001027161        0.00188579818734732     0.00188579818734732     2.07163983456953        1.09890587255635        0.76362 0.150.823734774258531
## 1:84295_G/A     A       1       G       0.001   149.504002380767        149.504002380767        0.00226260673135128     0.00226260673135128     0.000933350492394319    0.020360157782033       0.66655 0.976       0.0107158370534758
## 1:85011_C/T     T       1       C       0       62.718999342178 62.718999342178 0.000949194856561808    0.000949194856561808    0.128800949975553       -0.335748302235367      0.93568 0.72    0.142861466191716
## 1:85022_G/A     A       1       G       0.002   265.563001764938        265.563001764938        0.00401905384352773     0.00401905384352773     0.00360984204816058     0.0342513557618343      0.57017 0.952       0.0213218334725562

## library("stringr")
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

source("config.txt")


option_list <- list(
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),
  make_option("--analysis", type="character", default="BOLTLMM",
    help="Analysis software BOLTLMM or SAIGE [default='BOLTLMM]"),
  make_option("--minMAC", type="numeric", default=25,
    help="minimal minor allele count threshold [default=25]"),
  make_option("--hitregion", type="character", default="",
    help="File with candidate regions, CHROM;START;END;COL;LEGENDTEXT [default='']"),
  make_option("--outputClean", type="logical", default=F,
    help="Output MAF filtered and GC corrected results T/F [default=F]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
#print(opt)

# output files
file_prefix <- paste0(opt$prefix,"_",opt$analysis,".imputed")

file_top500 <- paste0(opt$prefix,"_",opt$analysis,".imputed.top500.txt")

# Settings for plot
plotLine <- T
gwthreshold <- 5E-8
yLine <- -log10(gwthreshold)
colLine <- "red"

summary.infile = snakemake@input[["summary"]]
imputed.infile = snakemake@input[["imputed"]]

# Get some basic information about phenotypes to determine title and MAF threshold
# read phenotype summary
phenoSummary <- read.table(summary.infile,sep="\t",comment.char="",
	header=T,check.names=F)
traitname <- phenoSummary[1,1]

binary <- ifelse(phenoSummary[1,"Summary_Controls"] == "",F,T)

# Read input
file_imp <- imputed.infile

if(opt$analysis == "BOLTLMM"){
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
} else if (opt$analysis == "SAIGE"){
	pcol <- "p.value"
	ycol <- "log10P_GC"
	chicol <- "CHISQ"
	betacol <- "BETA"
	secol <- "SE"
	chrcol <- "CHR"
	poscol <- "POS"
	idcol <- "ID"
	refcol <- ""
	altcol <- ""
	accol <- ""
	maccol <- ""
	afcol <- "AF"
	af0col <- ""
	af1col <- ""
} else {
	stop("Unknown Analysis")
}


# quantitative trait
if(!binary){
	cases <- as.integer(gsub("N=([0-9]+); Mean=.+","\\1",phenoSummary[1,"Summary_Cases"]))
	maintitle <- paste0("Phenotype '",traitname,"': ",
			format(cases,big.mark=",",scientific=FALSE),
			" Individuals")
	maf_threshold <- signif(opt$minMAC / (2*cases),2)
# binary trait
} else {
	cases <- as.integer(gsub("N=","",phenoSummary[1,"Summary_Cases"]))
	controls <- as.integer(gsub("N=","",phenoSummary[1,"Summary_Controls"]))
	maintitle <- paste0("Phenotype '",traitname,"': ",
			format(cases,big.mark=",",scientific=FALSE),
			" cases versus ",
			format(controls,big.mark=",",scientific=FALSE)," controls")
	maf_threshold <- signif(opt$minMAC / (2*min(cases,controls)),2)
}

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
resIn <- fread(file_imp)

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

## Clean GWAS output

if(opt$analysis == "BOLTLMM"){
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
} else if (opt$analysis == "SAIGE"){
	resOut <- data.frame(
		'CHROM'=resIn[[chrcol]],
		'POS'=resIn[[poscol]],
		'ID'=resIn$SNPID,
		'REF'=resIn$ALLELE0,
		'ALT'=resIn$ALLELE1,
		'AC'=resIn$AC,
		'MAC'=ifelse(resIn$AC > resIn$N,resIn$N*2-resIn$AC,resIn$AC),
		'MAF'=resIn$MAF,
		'AF'=resIn[[afcol]],
		'AF.CTRL'=resIn$AF.CTRL,
		'AF.CASE'=resIn$AF.CASE,
		'CHISQ'=resIn[[chicol]],
		'BETA'=resIn[[betacol]],
		'SEBETA_GC'=resIn$SEBETA_GC,
		'PVALUE_GC'=sapply(resIn[[ycol]],log10toP,USE.NAMES=F),
		'LOG10P_GC'=resIn[[ycol]])
}
write.table(resOut,paste0(file_prefix,".finalResults.txt"),
	sep="\t",col.names=T,row.names=F,quote=F)


# QQ plot Rscript downloaded from here: https://github.com/ilarsf/gwasTools
system(paste(
	"Rscript --vanilla","./QQplot.r",
	"--input",paste0(file_prefix,".finalResults.txt"),
	"--prefix",file_prefix,
	"--maf","MAF",
	"--pvalue","LOG10P_GC",
	"--log10p","T",
	"--maintitle", paste0("\'",gsub("[\"\']","`",maintitle),"\'")
	))

# Determine candidate regions if not provided
candidateHits <- which(resIn$log10P_GC > -log10(gwthreshold))
# BOLTLMM
if(opt$hitregion == "" & opt$analysis =="BOLTLMM" & length(candidateHits) >0 ){
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
		'LEGENDTEXT'="Candidate Regions [Top Hits ±(0.025cM + 50kb)]"
		)
# SAIGE
} else if (opt$hitregion == "" & opt$analysis =="SAIGE" & length(candidateHits) >0 ){
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
		'LEGENDTEXT'="Candidate Regions [±500kb]"
		)
} else if (opt$hitregion != ""){
	candidateRegions <- read.table(opt$hitregion)
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
	paste0(file_prefix,".hitregions.txt"),sep="\t",col.names=T,row.names=F,quote=F)

# Manhattan plot Rscript downloaded from here: https://github.com/ilarsf/gwasTools
system(paste(
	"Rscript --vanilla","./ManhattanPlot.r",
	"--input",paste0(file_prefix,".finalResults.txt"),
	"--prefix",file_prefix,
	"--hitregion",ifelse(dim(candidateRegions)[1] > 0,paste0(file_prefix,".hitregions.txt"),""),
	"--chr","CHROM",
	"--pos","POS",
	"--pvalue","LOG10P_GC",
	"--log10p","T",
	"--maintitle", paste0("\'",gsub("[\"\']","`",maintitle),"\'")
	))
