
## get subset of SNPs for PCA based on even spatial separation.
cat("This program takes the full starting set of SNPs and selects")
cat(" a subset according to a specifiable proportion--> 1:ratio")
cat(" ratio should be adjusted to give the desired proportion of")
cat(" SNPs taking into account genome gaps in the array, e.g, metabochip\n\n")

library(bigmemory)
library(biganalytics)
library(bigalgebra)
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

## settings ##
ratio <- 18.57 # achieved by trial and error to get desired SNP proportion (5%)
# e.g, 20 would be for even coverage, but slightly reduced to 
# 18.57 due to genome gaps on metabochip
writeResultsToFile <- F
genotypeFilt <- F  # run on a subset of the snps from the main file, e.g, filter out maf, etc
altSnps <- F  # whether to use alternative sublist if using genotype filt = T
initial5pc <-  T  # run the main selection component
big.descr <- "LRRFiltSortdescrFile"
out.fn <- "listOfSnpsForPCA.txt"

#### COMMAND ARGS ####
# process command line options (if present)
if (!exists("skipCmdArgs")) {
  arg.list <- parse.args(commandArgs(),coms=c("WRITE","GENFILT","ALT","INIT","RATIO","BFILE","OFILE"),
                         def=c(paste(as.integer(c(writeResultsToFile,genotypeFilt,altSnps,initial5pc))),
                               paste(ratio),big.descr,out.fn))
  if (!is.null(arg.list)) {
    if (arg.list["WRITE",]=="1") { writeResultsToFile <- T } else { writeResultsToFile <- F }
    if (arg.list["GENFILT",]=="1") { genotypeFilt <- T } else { genotypeFilt <- F }
    if (arg.list["ALT",]=="1") { altSnps <- T } else { altSnps <- F }
    if (arg.list["INIT",]=="1") { initial5pc <- T } else { initial5pc <- F }
    if (!is.na(as.numeric(arg.list["RATIO",1]))) { ratio <- as.numeric(paste(arg.list["RATIO",][1])) }
    if (nchar(arg.list["BFILE",])>1) { big.descr <- (paste(arg.list["BFILE",][1])) }
    if (nchar(arg.list["OFILE",])>1) { out.fn <- (paste(arg.list["OFILE",][1])) }
  }
}
######################
preload.bioC()  # load bioC packages without annoying warnings, etc

# if not in memory then read from file
sample.info <- read.table(file=dir.fn.cmb(dir$ano,"sample.info.tab"))

create.pheno.file <- F
if(create.pheno.file) {
  ## add t1d phenotype ## might add phenotype using a text file?
  # only necessary if using the association test..
  cases.grps <- c(2)
  phenotype <- integer(nrow(sample.info))
  phenotype[sample.info$grp %in% cases.grps] <- 1
  names(phenotype) <- rownames(sample.info)
  write.table(as.data.frame(phenotype),file=dir.fn.cmb(dir$ano,"pheno.lookup.txt"),quote=F)
}
phenotype <- read.table(file=dir.fn.cmb(dir$ano,"pheno.lookup.txt"))

sample.info <- add.to.sample.info(sample.info,phenotype,"phenotype")
cat("\nSample.info file preview:\n"); print(head(sample.info,4)) ; cat("\n")
write.table(sample.info,file=paste(dir$ano,"sample.info.tab",sep=""),quote=F)
cat("loaded object:",paste(load(paste(dir$ano,"samplensnp.RData",sep=""))),"\n")
snp.info <- sample.snp[[2]]
cat("\nSNP.info file preview:\n"); print(head(snp.info,4))
big.fn <- "combinedBigMat.RData"
snp.sub.fn <- "pca.snp.subset.txt"
#samp.fn <- "combined.samples.txt"
snps.to.keep <- extract.snp.subset(snp.info,sample.info,pc.to.keep=.05,assoc=F,
  writeResultsToFile=T,big.fn=big.fn,out.fn=snp.sub.fn,dir=dir)

bigMat <- getBigMat(big.fn,dir)
if(length(snps.to.keep)>100) {
  ##writeLines(colnames(bigMat),paste(dir$ano,samp.fn,sep=""))
  subset.descr <- big.exclude.sort(big.fn,dir,T,tranMode=1,pref="PCAMatrix",f.snp=snps.to.keep,verbose=F)
} else {
  stop("Error: list of snps to keep is too small - trying running again with a higher pc.to.keep\n")
}

save(subset.descr,file=paste(dir$big,"pcaSubMat.RData",sep=""))

#load.big.pointer.list(descr.list,dir)
  
  
#must.use.package("foreach")
#results[[dd]] <- foreach(i=1:ncolz, .combine='c') %do% ph.test(pheno,test.seg[,i])

# ph.test3 <- function(col,pheno) {
#  pv <- coefficients(summary(glm(pheno~col,family=binomial(logit))))[2,4]
#  return(pv)}
# ph.test2 <- function(col,pheno) { pv <- t.test(col~pheno)$p.value;   return(pv)}

system("R --slave < PCAPlusCorrect.R > F2.txt")