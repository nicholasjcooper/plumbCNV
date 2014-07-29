
library(bigmemory)
library(biganalytics)
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

## SETTINGS ##

## GC for each marker
vc.R.fn <- "vcfDataWithGC.RData"  # R binary version (optional)
big.descr <- "PCCORSortdescrFile"
# BAF bigmatrix file names 
pref <- "BAF"
# BAF for each marker
doBAF <- F
writeBAFfile <- F
writeGCfile <- F
writeOpticallFileFrag <- F
baf.out.fn <- paste(dir.cnv,"BAF.pfb",sep="")
gc.out.fn <- paste(dir.cnv,"marker.gcm",sep="")
oc.out.fn <- paste(dir.cnv,"opticallFirstRows.tab",sep="")
#snp.list.fn <- paste(dir.ano,"snp_rownames.txt",sep="") # 188052 snps, 3 = 160000; 2=140000

#### COMMAND ARGS ####
# process command line options (if present)
if (!exists("skipCmdArgs")) {
	arg.list <- parse.args(commandArgs(),coms=c("BAF","GC","OPTI","REBAF","LRRFILE","BAFPRE","VCF"),
                       def=c(paste(as.integer(c(writeBAFfile,writeGCfile,writeOpticallFileFrag,doBAF))),
                        big.descr,pref,vc.R.fn))
	if (!is.null(arg.list)) {
	  if (arg.list["BAF",]=="1") { writeBAFfile <- T } else { writeBAFfile <- F }
	  if (arg.list["GC",]=="1") { writeGCfile <- T } else { writeGCfile <- F }
	  if (arg.list["OPTI",]=="1") { writeOpticallFileFrag <- T } else { writeOpticallFileFrag <- F }
	  if (arg.list["REBAF",]=="1") { doBAF <- T } else { doBAF <- F }
	  if (nchar(arg.list["LRRFILE",])>1) { big.descr <- (paste(arg.list["LRRFILE",][1])) }
	  if (nchar(arg.list["BAFPRE",])>1) { pref <- (paste(arg.list["BAFPRE",][1])) }
	  if (nchar(arg.list["VCF",])>1) { vc.R.fn <- (paste(arg.list["VCF",][1])) }
	}
}
######################


bck.fn <- paste(pref,"bckfile",sep="")
des.fn <- paste(pref,"descrFile",sep="")
pass.samp <- paste(dir.ano,"SampleListLessStringentExclPlatesAndDuplicatesDLRS.txt",sep="")
#vcf.fn <- paste(dir.ano,"MetaboChip.markers.ncbi36.20101203.vcf",sep="")

##snp.list <- read.table(snp.list.fn,header=T,quote="\"",stringsAsFactors=F)

bigMat2 <- getBigMat(big.descr,dir.big)

snp.list <- rownames(bigMat2)

if(snp.list[1]=="Columns") { snp.list <- snp.list[-1] }

if ((writeBAFfile|writeGCfile|writeOpticallFileFrag)) {
    # load pre-existing vcf object from binary RData file
    print("Reading vcf file...")
    vcf.file <- reader(vc.R.fn,dir.ano) # should give vcf.file
} else  {
    stop("check options, VCF file not loaded")
}


if (doBAF | (writeBAFfile & (!"markerBAFrowMeans.RData" %in% dir.ano))) {
  # attach datafile bigmemory object
  bigMat2 <- getBigMat(des.fn,dir.big)
  print("Regenerating BAF average data from BAF matrix...")
  to.keep <- colnames(bigMat2) %in% reader(pass.samp)
  #mean.sel <- function(X,to.keep) { mean(X[to.keep],na.rm=T) }
  #marker.means <- apply(bigMat2,1,mean.sel,to.keep)
  marker.means <- rowMeans(bigMat2[,to.keep],na.rm=T)
  ofn <- paste(dir.ano,"markerBAFrowMeans.RData",sep="")
  save(marker.means,file=ofn)
  print(paste("Saved BAF rowmeans as:",ofn))
}


if(writeBAFfile) {
  # header: "Name","Chr","Position","PFB"
  # markername<TAB>chromosome<TAB>position<TAB>meanBAF
  print(load(paste(dir.ano,"markerBAFrowMeans.RData",sep=""))) # should give marker.means
  row.select <- match(snp.list,vcf.file$marker)
  baf.file <- vcf.file[row.select,c(2,1,4,3)]
  baf.file[,4] <- marker.means[snp.list]
  colnames(baf.file) <- c("Name","Chr","Position","PFB")
  print(paste("writing file baf means for each SNP to file:",baf.out.fn))
  write.table(baf.file,file=baf.out.fn,row.names=F,col.names=T,quote=F,sep="\t")
}

if(writeOpticallFileFrag) {
  # just write headers - then split separate file for each Chromosome and add X,Y
  # Chr SNP  Coor	Alleles	sample1A	sample1B	sample2A	sample2B	sample3A	sample3B
  # chr<TAB>markername<TAB>position<TAB>Alleles<TAB>X<TAB>Y<TAB>.... NaN for missing
  row.select <- match(paste(snp.list),vcf.file$marker)
  oc.file <- vcf.file[row.select,c(1,2,4,5,6)]
  oc.file[,4] <- paste(oc.file[,4],oc.file[,5],sep="")
  oc.file <- oc.file[,-5]
  colnames(oc.file) <- c("Chr","SNP","Coor","Alleles")
  print(paste("writing header columns for opticall file:",oc.out.fn))
  write.table(oc.file,file=oc.out.fn,row.names=F,col.names=T,quote=F,sep="\t")
}


if(writeGCfile) {
  # Name<TAB>Chr<TAB>Position<TAB>GC
  # markername<TAB>chromosome<TAB>position<TAB>meanGC
  #print(load(paste(dir.ano,"vcfDataWithGC.RData",sep=""))) # should give vcf.file
  #vcf.file...
  row.select <- match(snp.list,vcf.file$marker)
  gc.file <- vcf.file[row.select,c(2,1,4,3)]
  colnames(gc.file) <- c("Name","Chr","Position","GC")
  print(paste("writing SNP GC averages to file:",gc.out.fn))
  write.table(gc.file,file=gc.out.fn,row.names=F,col.names=T,quote=F,sep="\t")
}

print("Done.")




