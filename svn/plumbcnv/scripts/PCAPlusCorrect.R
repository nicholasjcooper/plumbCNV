library(xtable)
library(bigmemory)
library(biganalytics)
#library(biglm) one day maybe!

source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

if (!exists("actuallyDoPCA")) { actuallyDoPCA <- T }
pcaCentre <- T  # identify core samples from 2 dimensions of PCA
doPCAcorrection <- F
doStatsNGraphs <- F  # run some stats n' graphs immediate after PCA/PC correction
labelPCOutliersByPlatePost <- F
norm.pl.col <- "grey"
SD.thresh <- 2  # for 'core' definition of PCs in PCVersusPLATEPlots.R
# sample list for the QC-ed PCA # most likely un-needed
ifn <- "SampleListLessStringentExclPlatesAndDuplicatesDLRS.txt"
# PCs on the original (un-qc-ed) dataset
opc.fn <- "PCsEVsFromOrigPCA.RData" # "PCs.orig.data.allSamples.csv"  #don't have at the moment?
# info on plates, regions, etc for each sample
sinfo.fn <- "SubjectPlateRegionTypeSexEthnic.txt"
# file for PC data from this analysis (or to load for graphing, etc)
pcs.fn <- paste("PCsEVsFromPCA",".RData",sep="") #file name to save principle comps
if (!exists("num.pcs")) { num.pcs <- 9 } else { print(num.pcs) }
pre.n.pc <- 9  # number of PCs in the original PCA
post.n.pc <- num.pcs   # number of PCs in the current PCA
pcs.to.keep <- 100
pref <- "PCA"   # prefix for the 5% snp subset data file for PCA
prefC <- "PCCOR" # prefix for the full snp file for correction
regenHisto <- F

#### COMMAND ARGS ####
# process command line options (if present)
if (!exists("skipCmdArgs")) {
	arg.list <- parse.args(commandArgs(),coms=c("HISTO","LAB","PCA","COR","FIG","CORE","RAW","PCS","STORE", "PRE","PREC","SFILE","OFILE","DMFILE","PCFILE"),def=c(as.character(as.integer(c(regenHisto,labelPCOutliersByPlatePost, actuallyDoPCA, doPCAcorrection, doStatsNGraphs, pcaCentre, pre.n.pc, post.n.pc, pcs.to.keep))), pref, prefC, ifn, opc.fn, sinfo.fn, pcs.fn))
	  if (!is.null(arg.list)) {
	  if (arg.list["HISTO",]=="1") { regenHisto <- T } else { regenHisto <- F }
	  if (arg.list["LAB",]=="1") { labelPCOutliersByPlatePost <- T } else { labelPCOutliersByPlatePost <- F }
	  if (arg.list["PCA",]=="1") { actuallyDoPCA <- T } else { actuallyDoPCA <- F }
	  if (arg.list["COR",]=="1") { doPCAcorrection <- T } else { doPCAcorrection <- F }
	  if (arg.list["FIG",]=="1") { doStatsNGraphs <- T } else { getCatData <- F }
	  if (arg.list["CORE",]=="1") { pcaCentre <- T } else { pcaCentre <- F }
	  if (arg.list["RAW",1] %in% paste(c(1:200))) { pre.n.pc <- as.integer(paste(arg.list["RAW",][1])) }
	  if (arg.list["PCS",1] %in% paste(c(1:200))) { post.n.pc <- as.integer(paste(arg.list["PCS",][1])) }
	  if (arg.list["STORE",1] %in% paste(c(1:1000))) { pcs.to.keep <- as.integer(paste(arg.list["STORE",][1])) }
	  if (nchar(arg.list["PRE",])>1) { pref <- (paste(arg.list["PRE",][1])) }
	  if (nchar(arg.list["PREC",])>1) { prefC <- (paste(arg.list["PREC",][1])) }
	  if (nchar(arg.list["SFILE",])>1) { ifn <- (paste(arg.list["SFILE",][1])) }
	  if (nchar(arg.list["OFILE",])>1) { opc.fn <- (paste(arg.list["OFILE",][1])) }
	  if (nchar(arg.list["DMFILE",])>1) { sinfo.fn <- (paste(arg.list["DMFILE",][1])) }
	  if (nchar(arg.list["PCFILE",])>1) { pcs.fn <- (paste(arg.list["PCFILE",][1])) }
	}
}
######################

preload.bioC()  # load bioC packages without annoying warnings, etc

pcs.fn <- paste(pcs.fn) #file name to save principle comps
sinfo.fn <- paste(dir$ano,sinfo.fn,sep="")
#opc.fn <- paste(dir$pcO,opc.fn,sep="")  #don't have at the moment?
opc.fn <- paste(opc.fn)  #don't have at the moment?
snp.info.fn <- paste(dir$ano,"rawdata.map",sep="")

descr.fn <- "combinedBigMat.RData"
subDescr <- get(paste(load(file=dir.fn.cmb(dir$big,"pcaSubMat.RData"))))

pca.result <- LRR.PCA(subDescr,pcs.to.keep=10,SVD=F,LAP=F,saver=T,pcs.fn="PCsEVsFromPCA.RData") 
varz <- do.scree.plots(pca.result$Evalues,dir,elbow=num.pcs)
cat("",round(cumsum(varz)[num.pcs]*100,1),"% LRR variance explained by first",num.pcs,"components.\n")
corrected.ref <- LRR.PCA.correct(pca.result,descr.fn,num.pcs=9,pref="corrected",write=T)
corrected.bigMat <- getBigMat(corrected.ref,dir) ; rm(corrected.ref) #(in case it was a bigmat)

out <- lrr.sample.dists(corrected.bigMat,snp.info.fn,mode="map3",gc=F,dlrs=F,dir=dir,pref="PostPC")

cleaned.stat.table <- out$stat.table # could here compare with the pre-cleaned one.

sample.info <- read.table(file=paste(dir$ano,"sample.info.tab",sep=""))

## optionally creates 3 x 3 plots of raw,sample.qc,pc.corrected data, returns the 3 stat tables
twc <- three.way.comparison(cleaned.stat.table,sample.info,dir=dir,batch.comps=c("plate"),
                     height=10,width=10) 
   #,"grp","phenotype"),
  

system("cat F0.txt F1.txt F2.txt > F_success.txt")
