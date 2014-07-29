## read in big matrix, remove unwanted samples/snps and order
print("This program reads in a big matrix for LRR and BAF, combines them into sample files for penn cnv")

library(bigmemory)
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

# prepare files for penn cnv
# format is:
# Name<TAB>sampleid.Log R Ratio<TAB>sampleid.B Allele Freq
# filenames will be 1001 in each dir: sampleid.txt
#### BIGMATRIX FILES CURRENTLY IN USE HAVE DIMENSIONS AND INDEX FILES:
# all samples, e.g, "subIds.txt"
# BAF: 188052   5841 snp_rownames.txt
# QC passing, non dupl, not bad plt: "SampleListLessStringentExclPlatesAndDuplicatesDLRS.txt"
# toPCMCor9/BAF9:  162443   5005 snp_rownames3.txt
# toPC2Cor9:   137011   5005 snp_rownames2.txt

dir.sub <- "PENNRAWFILES"  #paste(dir.cnv,"FilteredMAF0/",sep="")
baf.txt <- ".B Allele Freq"  
lrr.txt <- ".Log R Ratio"
BAF.fn <- c("BAFdescrFile")  # should have all snps and samples in LRR.fn (but also may have others)
LRR.fn <- c("PCCORSortdescrFile")  # this file drives the samples and snps included (inc order)

#### COMMAND ARGS ####
# process command line options (if present)
if (!exists("skipCmdArgs")) {
	arg.list <- parse.args(commandArgs(),coms=c("DIR","BAF","LRR","N","TEST"),
                       def=c(paste(as.integer(c(custom.n))),
                        dir.sub,BAF.fn,LRR.fn,ff))
	if (!is.null(arg.list)) {
	  if (arg.list["N",1] %in% paste(c(0:10^6))) { custom.n <- as.integer(paste(arg.list["N",][1])) }
	  if (arg.list["TEST",]=="1") { ff <- T } else { ff <- F }
	  if (nchar(arg.list["DIR",])>1) { dir.sub <- (paste(arg.list["DIR",][1])) }
	  if (nchar(arg.list["BAF",])>1) { BAF.fn <- (paste(arg.list["BAF",][1])) }
	  if (nchar(arg.list["LRR",])>1) { LRR.fn <- (paste(arg.list["LRR",][1])) }
	}
}
######################


prepare.penncnv.samples <- function(LRR.fn,BAF.fn,dir,low.ram=T) 
{
  dir <- validate.dir.for(dir,c("cnv","big"))
  # default PENNCNV header suffixes
  baf.txt <- ".B Allele Freq"  ;  lrr.txt <- ".Log R Ratio"
  subdir <- paste(dir$cnv,"PENNRAWFILES",sep="") 
  bigBAF <- getBigMat(BAF.fn, dir$big)
  bigLRR <- getBigMat(LRR.fn, dir$big)
  printBigSummary(bigBAF,"bigBAF") ; printBigSummary(bigLRR,"bigLRR")
  rB <- rownames(bigBAF); cB <- colnames(bigBAF)
  rL <- rownames(bigLRR); cL <- colnames(bigLRR)
  n <- ncol(bigLRR)
  ndirs <- round(n/1000)
  dir.sizes <- (n %/% ndirs) + c(rep(1,n %% ndirs),c(rep(0,ndirs-(n %% ndirs))))
  print(paste("Will create",ndirs,"penn folders with sizes",paste(dir.sizes,collapse=",")))
  
  ##BAF.nms <- paste(cB,baf.txt,sep="") #used where??
  ##LRR.nms <- paste(cL,lrr.txt,sep="")  #used where??
  baf.col.sel <- match(cL,cB) ; cat(length(baf.col.sel[!is.na(baf.col.sel)]),"cols from LRR in BAF\n")
  baf.row.sel <- match(rL,rB) ; cat(length(baf.row.sel[!is.na(baf.row.sel)]),"rows from LRR in BAF\n")
  lrr.col.sel <- which(!is.na(baf.col.sel))
  cat("BAF missing",length(which(is.na(baf.col.sel))),"samples in LRR file\n")
  lrr.row.sel <- which(!is.na(baf.row.sel))
  cat("BAF missing",length(which(is.na(baf.col.sel))),"snps in LRR file\n")
  grand.rn <- rL[lrr.row.sel]
  ff <- FALSE # initialise failure flag
  if(!all(grand.rn==rB[baf.row.sel])) 
    { print("warning: row name selections don't match") ; ff <- T }
  tot <- n.to.process <- sum(dir.sizes)
  if(n.to.process!=length(lrr.col.sel)) {   
    ff <- T ;  cat("Warning: number of samples per directory and number of directories")
    cat(" does not imply the number of common columns in the BAF,LRR sample matrices\n") }
  if(!ff){
    for (nn in 1:length(dir.sizes))
    {
      next.dir <- paste(dir,"/p",nn,sep="")
      dir.create(next.dir)
      set.offset <- sum(dir.sizes)-sum(dir.sizes[5:nn])
      for (mm in 1:dir.sizes[nn])
      {
        #start.time <- proc.time()[3]
        loc <- set.offset+mm
        baf.sel <- baf.col.sel[loc]
        lrr.sel <- lrr.col.sel[loc]
        sampleid <- cL[lrr.sel]; s2 <- cB[baf.sel] 
        if(sampleid!=s2) { ff <- T ; stop("Error: col indexes failed") ; break }
        cat("Writing #",loc," dir ",nn," sample ",mm," ID:",sampleid,"\n")
        header <- paste("Name",paste(sampleid,lrr.txt,sep=""),
                           paste(sampleid,baf.txt,sep=""),sep="\t")
        to.write <- paste(grand.rn,bigLRR[lrr.row.sel,lrr.sel],
                              bigBAF[baf.row.sel,baf.sel],sep="\t")
        conn <- file(paste(next.dir,"/",sampleid,".txt",sep=""),"w")
        writeLines(header,con=conn)
        writeLines(to.write,con=conn)
        close(conn)
        n.to.process <- n.to.process-1
        loop.tracker((tot-n.to.process),tot)
      }
      if(low.ram) {
        # flush and reload bigmatrices to clear RAM
        rm(bigBAF); rm(bigLRR); gc()
        bigBAF <- getBigMat(BAF.fn, dir$big)
        bigLRR <- getBigMat(LRR.fn, dir$big)
      }
    }
  }
  return(dir.sizes)
}


LRR.fn <- paste("describePCcorrect",num.pcs,".RData",sep="")

dir.sizes <- prepare.penncnv.samples(LRR.fn,BAF.fn,dir,low.ram=T) 
cat("number of samples written per PENN-CNV directory")
print(dir.sizes)


