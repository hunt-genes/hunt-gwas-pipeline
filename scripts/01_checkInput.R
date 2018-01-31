options(stringsAsFactors=F)
library("data.table")
library("GenABEL")


pheno = snakemake@input[["pheno"]]
key = snakemake@input[["key"]]

outdir = snakemake@params[["outdir"]]
normalize = snakemake@params[["normalize"]]

outfile = snakemake@output[["checked"]]


phenoIn <- read.table(pheno,sep="\t",header=T,check.names=F,comment.char="")
phenoKey <- read.table(key,sep="\t",header=T,check.names=F,comment.char="",quote="\"")

# check columns of key file
keycols <- c("name","version","type","baseQcovar","baseBcovar","addQcovar","addBcovar","description","N","N_cases","N_controls")
if(length(which(!keycols %in% colnames(phenoKey)))>0) stop("Wrong / missing columns in key file")

setnames(phenoKey,c("N","N_cases","N_controls"),c("original.N","original.N_cases","original.N_controls"))
phenoKey$N <- NA
phenoKey$N_cases <- NA
phenoKey$N_controls <- NA

phenoKey <- phenoKey[which(phenoKey$type %in% c("binary","categorical","quantitative")),]

for(i in 1:dim(phenoKey)[1]){
    traitName <- phenoKey$name[i]

    pedcols <- c("FID","IID","PATID","MATID")
    bcovars <- as.character(na.omit(unlist(strsplit(unlist(phenoKey[i,c("baseBcovar","addBcovar")]),"\\|"))))
    qcovars <- as.character(na.omit(unlist(strsplit(unlist(phenoKey[i,c("baseQcovar","addQcovar")]),"\\|"))))

    ncovars <- as.character(na.omit(unlist(strsplit(unlist(phenoKey[i,c("addQcovar","addBcovar")]),"\\|"))))
    covars <- c(bcovars,qcovars)

    if(length(which(!c(pedcols,traitName,covars) %in% colnames(phenoIn))) > 0) next
    phenoTemp <- phenoIn[,which(colnames(phenoIn) %in% c(pedcols,traitName,covars))]
    phenoTemp <- na.omit(phenoTemp)

    if(dim(phenoTemp)[1] == 0) next
    if(phenoKey$type[i] == "binary"){
        x <- table(phenoTemp[[traitName]])
        if(length(x) == 1) next
        if(length(which(names(x) %in% 1:2))==2) phenoTemp[[traitName]] <- phenoTemp[[traitName]] - 1
    } else if (phenoKey$type[i] == "quantitative" & normalize){
        phenoTemp[[traitName]] <- rntransform(as.formula(paste(phenoKey$name[i],"~",paste(ncovars,collapse=" + "))),phenoTemp,family=gaussian)
    }

    phenoKey$N[i] <- dim(phenoTemp)[1]

    checkInput <- data.frame('Column'=c(traitName,covars),
                             'Type'=c("Phenotype",rep("bCovariate",length(bcovars)),rep("qCovariate",length(qcovars))),
                             'N_withPheno'=apply(phenoIn[which(!is.na(phenoIn[[traitName]])),c(traitName,bcovars,qcovars)],2,
                                function(x) length(which(!is.na(x)))),
                             'Summary_Cases'=character(1),
                             'Summary_Controls'=character(1),
                            check.names=F)

    if(phenoKey$type[i] %in% c("categorical","quantitative")){
        keepID <- which(phenoIn$IID %in% phenoTemp$IID)
        checkInput$Summary_Cases[1] <- paste("N=",length(keepID),"; Mean=",signif(mean(phenoIn[keepID,checkInput$Column[1]],na.rm=T),3),
                                  "; SD=",signif(sd(phenoIn[keepID,checkInput$Column[1]],na.rm=T),3),
                                  "; [",paste(signif(range(phenoIn[keepID,checkInput$Column[1]],na.rm=T),3),
                                  collapse=";"),"]",sep="")
        for(j in 2:dim(checkInput)[1]){
            x <- table(phenoTemp[,checkInput$Column[j]])
            if(length(x) > 2){
               checkInput$Summary_Cases[j] <- paste("Mean=",signif(mean(phenoTemp[,checkInput$Column[j]],na.rm=T),3),
                                              "; SD=",signif(sd(phenoTemp[,checkInput$Column[j]],na.rm=T),3),
                                              "; [",paste(signif(range(phenoTemp[,checkInput$Column[j]],na.rm=T),3),
                                              collapse=";"),"]",sep="")
            } else {
               checkInput$Summary_Cases[j] <- paste(names(x),as.numeric(x),collapse="; ",sep="=")
            }
        }
    } else {
        cases <- which(phenoTemp[,traitName] == 1)
        controls <- which(phenoTemp[,traitName] == 0)
        checkInput$Summary_Cases[1] <- paste0("N=",length(cases))
        checkInput$Summary_Controls[1] <- paste0("N=",length(controls))

        phenoKey$N_cases[i] <- length(cases)
        phenoKey$N_controls[i] <- length(controls)

        for(j in 2:dim(checkInput)[1]){
            cases <- which(phenoTemp$IID %in% phenoTemp$IID & phenoTemp[,traitName])
            x <- table(phenoTemp[,checkInput$Column[j]])
            if(length(x) > 2){
               checkInput$Summary_Cases[j] <- paste("Mean=",signif(mean(phenoTemp[cases,checkInput$Column[j]],na.rm=T),3),
                                              "; SD=",signif(sd(phenoTemp[cases,checkInput$Column[j]],na.rm=T),3),
                                              "; [",paste(signif(range(phenoTemp[cases,checkInput$Column[j]],na.rm=T),3),
                                              collapse=";"),"]",sep="")
               checkInput$Summary_Controls[j] <- paste("Mean=",signif(mean(phenoTemp[controls,checkInput$Column[j]],na.rm=T),3),
                                              "; SD=",signif(sd(phenoTemp[controls,checkInput$Column[j]],na.rm=T),3),
                                              "; [",paste(signif(range(phenoTemp[controls,checkInput$Column[j]],na.rm=T),3),
                                              collapse=";"),"]",sep="")
            } else {
               x <- table(phenoTemp[,traitName],phenoTemp[,checkInput$Column[j]])
               checkInput$Summary_Cases[j] <- paste(colnames(x),x[2,],collapse="; ",sep="=")
               checkInput$Summary_Controls[j] <- paste(colnames(x),x[1,],collapse="; ",sep="=")
            }
        }
    }

    file_phenoOut <- paste0(outdir,"/",traitName,"_pheno.txt")
    write.table(checkInput,paste0(outdir,"/",traitName,"_summary.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    write.table(phenoTemp[,c(pedcols,traitName,covars)],file_phenoOut,sep="\t",col.names=T,row.names=F,quote=F)
}

write.table(phenoKey,outfile,sep="\t",col.names=T,row.names=F,quote=F)
