## plotting / tabulating stats on the LRR means, medians, SDs 
# across samples
cat("\nThis module converts text LRR data to a bigmatrix file, applies exclusions")
cat("from call-rate QC, then performs sample-QC, evaluates and produces cleaned")
cat(" matrices for each input datafile as a list.\n\n")

library(limma)
library(xtable)
library(bigmemory)
library(biganalytics)
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

if (exists("mM")) { mainMode <- mM } else { mainMode <- 1 }

## DEFAULT SETTINGS ##
doPlots <- F  # set to F to skip generating many plots
writeExclFiles <- F  # write exclusions to files or just generate the lists.
wav.filt.mode <- 2  # whether to use: 1: R_GC ; 2:S_WF [wave from all sources]  3:S_GCWF [wave from GC sources]
scl <- 1000000 # graph scale factor (bases)
# dataset to use [nb: pc plate plots will want this as 3 different opts hence 'exists' condition]
if(!exists("big.descr")) {
  big.descr <- "LRRFiltSortdescrFile"   }  # bigData used for plotting extreme samps + chrinfo

pref <- "1"
med.chunk.fn <- "chunkMedianInd.RData"
snp.info.fn <- paste(dir$ano,"rawdata.map",sep="")

des.fn <- "LRR1SortdescrFile.RData"  
dlrs.pref <- "DLRSq"
#qc.fn <- paste(dir$qc.lrr,"QCstats",pref,".txt",sep="")
hib <- 2.5 #2.575 #3.29 # define stringent z score boundary 
lob <- 2 #1.96 #2.575 # define lenient z score boundary
input.dir <- dir$col #"/chiswick/data/ncooper/ImmunochipReplication/LRRDATA/ColumnData/"
input.suf <- "LRR.dat"
sample.fn <- "subIdsALL.txt"
snp.fn <- "snpNames.txt"
doPlateCols <- F

#### COMMAND ARGS ####
# process command line options (if present)
if (!exists("skipCmdArgs")) {
	arg.list <- parse.args(commandArgs(),coms=c("M","GCF","PLOTS","WRITE","COLS","BFILE"),
         def=c(paste(as.integer(c(mainMode,wav.filt.mode,doPlots,writeExclFiles,doPlateCols))),
                        big.descr))
	if (!is.null(arg.list)) {
	  if (arg.list["M",1] %in% paste(c(1,2,3))) { mainMode <- as.integer(paste(arg.list["M",][1])) }
	  if (arg.list["GCF",1] %in% paste(c(1,2))) { wav.filt.mode <- as.integer(paste(arg.list["GCF",][1])) }
	  if (arg.list["PLOTS",]=="1") { doPlots <- T } else { doPlots <- F }
	  if (arg.list["WRITE",]=="1") { writeExclFiles <- T } else { writeExclFiles <- F }
	  if (arg.list["COLS",]=="1") { doPlateCols <- T } else { doPlateCols <- F }
	  if (nchar(arg.list["BFILE",])>1) { big.descr <- (paste(arg.list["BFILE",][1])) }
	}
}
######################

###########################################
## PART ONE - PROCESS DATA / ANNOTATION  ##
###########################################

process.cohort.qc <- function(grp,of,dir,snp.fn,mode="map3",snp.info.fn,dlrs.pref,
                              med.chunk.fn,plate.lookup,gotMedian=F,gotBig=F,badPlateThresh=0.33,...)
{
  ### loop through once for each file
  ### passes arguments to lrr.sample.dists which passes to calc.chr.info
  dir <- validate.dir.for(dir,c("big","col","ano","qc.cab","qc.pl","qc.lrr","qc.gc"),warn=T)
  med.chunk.fn <- paste(dir$qc.gc,pref,grp,med.chunk.fn,sep="")
  if(of==1) {
    pref <- "LRR" #  sample.fn <- "subIdsALL.txt"
  } else {
    pref <- paste("LRR",grp,sep="")
  }

  des.fn <- paste(pref,"descrFile",sep="")
  ## LOAD RAW DATA FROM SOURCE##
  if(!gotBig | !file.exists(dir.fn.cmb(dir$big,des.fn))) {
    des.fn <- plain.vec.to.big.matrix(dir,grp,snp.fn=snp.fn,pref=pref)
  } else {
    cat("Loading big matrix from descriptor file",des.fn,":\n")
  }
  ## FILTER ON SNP-QC, call rate analysis prepared earlier, etc ##
  descr <- big.exclude.sort(des.fn, dir=dir, pref=pref)
  bigMat2 <- getBigMat(descr,dir)
  printBigSummary(bigMat2,"Snp_QC_Mat")
  
  datlist <- lrr.sample.dists(bigMat2,snp.info.fn,mode="map3",gc=T,dlrs=T,dir=dir,pref=pref,med.chunk.fn=med.chunk.fn,...)
  
  s.tab <- lrr.stats.tab(datlist$stat.table)
  cat("\nCohort distribution statistics for LRR-stats:\n"); print(s.tab); cat("\n")
  rez <- get.excl.list.from.stats(s.tab,datlist$stat.table,dir,writeFiles=T,append=T) # append!
  do.venn.QC(rez,dir,pref) # DRAW VENN DIAGRAM
    
  ## do chromosomal abberations
  cat("\nDetecting chromosomal abberations\n")
  chr.stat <- get.chr.stats(bigMat2,datlist$CHR.INFO,dir)
  chrWarns <- do.chr.ab.contrasts(chr.stat$chr.mean,lob,hib)
  excluded.samps <- sort(table(unlist(chrWarns)))
  chr.ab.report(chr.stat, chrWarns, excluded.samps, dir=dir, 
                makeGraphs=F, writeExclList=T,append=T)
  
  ### PLATE STATS ###
  new.fr <- update.plate.bad.count.table(dir,plate.lookup,filt.list=colnames(bigMat2))
  ## yields here total count set for this plate only (because of filt.list above):
  this.fr <- new.fr[!is.na(new.fr$SIZE),]
  cat("\nTable of Plate failure counts:\n"); print(this.fr)
  ofn <- paste(dir$qc.pl,"CountsPerPlate",pref,".tab",sep="")
   write.table(this.fr,sep="\t",col.names=T,row.names=F,quote=F,file=ofn)
  cat("written to",ofn,"\n")
  # the update.plate.... function....may also be useful when combining all file results
  
  ## PLATEWISE EXCLUSION - writes excluded samples to annotation directory
  # removes plates with more than one third samples excluded by LRR sample QC
  excl.info <- excl.bad.plates(this.fr,dir,badPlateThresh=0.33,writeExclList=T,append=T)
    
  ### PLOTS / REPORTS ###
  # now compile venn list and do failing example plots
  plot.extreme.samples(rez, datlist$stat.table, bigMat2, datlist$CHR.INFO, dir, pref)
  ## function to do those scatter plots
  lrr.boundary.scatter(datlist$stat.table,s.tab,dir,fn.pre=pref,fn.suf="",col="darkblue")
  fancy.dens.plot <- F
  ## do those density plots
  if (fancy.dens.plot) {
    plt <- plate.lookup
    grop <- plt[plt$id %in% colnames(bigMat2),"plate"]  # SUBINFO[match(rownames(stat.table),SUBINFO[,1]),5]
    batch.effects.density.plot(stat.table, grp.dat=grop, grp.lab="Plate", npg=8, 
                             ylz=c(50,50), xllz=c(-.2,.0), xlhz=c(.1,.5),dir=dir,extrapref=grp)
  }
  plate.lrr.stats <- round(get.plate.lrr.stats(plate.lookup,datlist$stat.table,colnames(bigMat2)),3)
  batch.box.plots(datlist$stat.table,s.tab,plate.lookup,batch="plate",pref=pref,dir=dir)
  ofn <- paste(dir$qc.pl,"StatsPerPlate",pref,".tab",sep="")
   write.table(plate.lrr.stats,sep="\t",col.names=T,row.names=T,quote=F,file=ofn)
  cat("table of LRR-statistics by Plate, written to:\n",ofn,"\n")
  descr.samp.filt <- descr #big.exclude.sort(descr, dir=dir, pref=paste("SampQC_LRR",grp,sep=""))
  return(descr.samp.filt)
}


preload.bioC()  # load bioC packages without annoying warnings, etc

# delete any existing sample exclusion files
sure <- T
if(sure) {initialise.excl.files(dir) }

## gather sample information (firstly loading from LoadInSNPs.R result)
sample.snp <- get(load(file=dir.fn.cmb(dir$ano,"samplensnp.RData",must.exist=T)))

## make sure sample info is nice, in particular that 'grp' matches 'GRP' in 'file.spec.txt'
sample.info <- validate.samp.info(sample.snp$sample.info,dir)

# load plate info
plate.lookup <- get.plate.info(dir,1,2,3,dup.action="trim")  # 1,4,5
####plate.lookup <- trim.plate.duplicates(plate.lookup,by="counts")
# add plate annotation to sample.info
###sample.info <- add.plates.to.sample.info(sample.info,plate.lookup)
sample.info <- add.to.sample.info(sample.info,plate.lookup[[1]],to.add=c("plate","well"))

## write unified get.SubIDs function that syncs with file.spec.txt in allowed ways
# should also allow to retrieve file names, separate lists, combined list, grp list
# read/write sample.info and snp.info , facilitate using these objects as input always.
# file.spec.txt should be a pinion too

cat("\nRunning sample QC for cohort(s)\n")

file.info <- get.file.specs(dir)

ngrp <- length(unique(file.info$GRP))  #length(grp.fnz)
descr.list <- vector("list", ngrp)

for (grp in 1:ngrp) {
  descr.list[[grp]] <- process.cohort.qc(grp,of=ngrp,dir,snp.fn,mode="map3",snp.info.fn,
      dlrs.pref,med.chunk.fn,plate.lookup[[1]],gotBig=T)
}

descr.fn <- paste(dir$big,"grpsBigMats.RData",sep="")
save(descr.list,file=descr.fn)
  
new.descr <- exclude.combine.big.mats(descr.list,dir) 
save(new.descr,file=paste(dir$big,"combinedBigMat.RData",sep=""))

  
#BML <- (load.big.pointer.list(descr.list,dir))$bigMatLst

#####make sure files/pics produced each time have unique names, make sure QC info is appended not overwritten

## run combined QC script as a function and generate table, final sample lists to pass to PC.

summary.list <- make.QC.summary.table(sample.list=rownames(sample.info),dir)
passQCsamples <- summary.list$pass.samples
  
sample.info$QCfail[rownames(sample.info) %in% get.all.samp.fails(dir)] <- 1

cat("\nFirst 20 samples excerpt:\n")
print(head(sample.info,20))
ofn <- paste(dir$ano,"sample.info.tab",sep="")
 write.table(sample.info,file=ofn,quote=F)
cat(paste("==> full sample-set written to:",ofn,"\n"))


system("R --slave < SubsetOfMarkers.R > F1.txt")
