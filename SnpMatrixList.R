###NAMESPACE ADDITIONS###
#' @importMethodsFrom "snpStats"  effect.sign
#' @importFrom "snpStats"  row.summary  col.summary  read.pedfile  snp.imputation  impute.snps  single.snp.tests
#' @importClassesFrom "snpStats"  SnpMatrix  XSnpMatrix  SingleSnpTests  SingleSnpTestsScore
#' @imports NCmisc reader bigpca humarray bigmemory biganalytics annotSnpStats snpStats GenomicRanges IRanges genoset

library(annotSnpStats)
library(humarray)
library(reader)
library(bigpca)

#### internals ####
# rbind.sml # untested, could export if it were good
# get.grp.level.snp.stats
# ld.prune.big
# colmeansel
# correct.by
# rmv.common.missing
# samp.from.annot
# snp.from.annot
# is.aSnpMatrix
# split.vec
# get.split.type
# apply.snp.thresholds
# sync.snpmat.with.info
# get.ichip.1000 # internal or specific
# string.in.common
# Get.chr
# get.chr
# make.split
# randomize.missing2
# logsum
# comma  # NCmisc?
# sortna
# sumna
# sdna
# medianna
# meanna
# maxna
# minna
# finitize
# rsnp
# rsnpid
# rsampid
# ldfun
# snpify.cont
# get.biggest
# get.top.n
# clean.fn
# read.ped.file
# validate.dir.for
# get.allele.counts
# Y_2
# X_2
# Likelihood_Lj
# LLikelihood_L
# get.SnpMatrix.in.file
# sync.snpmat.with.info
# fun.snpMatLst
# clean.snp.ids
# add.missing.snps
# choose.dups.to.drop

## eventually turn this into a class 'snpMatrixList'

## FUNCTIONS ##
#cbind3 - cbind for more than 2 [a]SnpMatrix objects at once
#rbind3 - rbind for more than 2 [a]SnpMatrix objects at once
#snp.mat.list.type - determine whether a snpMatLst is a list of filenames, SnpMatrix objects, or just a [a]SnpMatrix
#sampSel - select all SNPs for a subset of samples in a snpMatLst
#snpSel - select all samples for a subset of SNPs in a snpMatLst
#string.in.common - take a set of strings and extract any text common to all, or only the unique text. e.g, unique text can be used to find differing pre-fixes or suffixes from the same root file name
#read.impute.list - read impute file(s) into a SnpMatrix[List] format multiple files assumed to be same samples, different SNP sets (e.g, chromsomes)
#fun.snpMatLst - apply a function to each element of a snpMatLst
#get.SnpMatrix.in.file - specific to me? reads snp matrices from RData binary file. if it finds multiple, will attempt to join them together can also handle an XSnpMatrix
#check.snpmat.size - check whether a specific number of samples and snps will fit in a SnpMatrix object
#snpMatLst.collate - convert a snpMatLst to a single SnpMatrix
#get.snpMat.spec - get the dimensions of all of the SnpMatrix objects in a snpMatLst, as a matrix, where cols are datasets, first row is nsamp, second nsnp
#colnamesL - get a list of the 'colnames' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, combined set of colnames)
#rownamesL - get a list of the 'rownames' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, combined set of rownames)
#nrowL - get a list of the 'nrow' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, total nrow)
#ncolL - get a list of the 'ncol' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, total ncol)
#smlapply - 'lapply' function for SnpMatrixList, also works for normal SnpMatrix or aSnpMatrix or XSnpMatrix
#sync.snpmat.with.info - for a snpMatLst, sync with snp.info and sample.info to exclude snps/samples not in the support files
#sample.info.from.sml - internal, make a sample.info object using either just names, or using a snpMatLst too
#doSampQC - run sample quality control on a SnpMatrix or snpMatLst
#doSnpQC - run SNP quality control on a SnpMatrix or snpMatLst
#apply.snp.thresholds - apply the thresholds from doSnpQC() to a dataset to get exclusion lists
#get.split.type - detect the format of a snpMatLst (SnpMatrix list) returns snpSubGroups, sampleSubGroups, or singleEntry (else error)
#list.rowsummary - snpStats::row.summary for a snpMatLst (sample QC summary)
#list.colsummary - snpStats::col.summary for a snpMatLst (SNP-qc sumamry)
#convert.smp.to.chr22 - convert snpMatLst which is split by sample subsets to one split by chromosome
#convert.chr22.to.smp - convert snpMatLst which is split by chromsome to one split by sample subsets
#is.aSnpMatrix - detect whether object is aSnpMatrix type, even if within a snpMatLst
#sample.info.from.annot - create a full sample.info object from a list of aSnpMatrix objects
#snp.info.from.annot - create a full snp.info object from a list of aSnpMatrix objects
#snp.from.annot - create a (partial) snp.info object from a single aSnpMatrix object
#samp.from.annot - create a (partial) sample.info object from a single aSnpMatrix object 
# annot.sep.support - Create aSnpMatrix from SnpMatrix plus SNP/sample support files
# split.pq - split an aSnpMatrix into a long and short arm (arm.p and arm.q)
###############
# lambdas Calculate Lambda inflation factors for SNP dataset This function calculates SNP-wise or overall Lambda and Lambda_1000 statistics for inflation due to population structure. It works on a SnpMatrix object or dataframe coded 0,1,2,NA (autodetects which).
# SnpMatrix.to.data.frame Convert a snpStats SnpMatrix object to a dataframe Converts a snpStats::SnpMatrix object to a dataframe where coding becomes 0,1,2,NA, which represents genotypes as the number of copies of the reference allele.
# df.to.SnpMatrix Convert a data.frame to a snpStats SnpMatrix object Converts a dataframe to a snpStats::SnpMatrix object where the object contains genotypes coded as number of copies of the reference allele: 0,1,2, and missing=NA. This is an alternative to using new("SnpMatrix",data.frame()). Using 'new' the required format for the 'data.frame' argument is not as intuitive, as NA's are not allowed, and 2 copies of the reference allele must be coded as 3, and 1 copy as 2, 0 copies as 1. Note that this function will also accept data.frames/matrices coded in that way, and will detect the coding automatically.
# majmin Determine major or minor allele status for a set of SNPs For a snpStats object or data.frame containing values 0,1,2, NA representing genotypes AA, AB, BB and no-call. Determines whether the reference allele is the major or minor allele Where the homozygous genotype coded as highest value = reference, e.g, if AA=0, AB=1, BB=2, then B is considered the reference here, and by the snpStats package. Combines this with frequencies of the alleles to evaluate whether 'BB' is major or minor. Note that default behaviour for a SnpMatrix is to code alleles alphabetically, so usually the reference allele is the letter later in the alphabet, e.g, it is never an 'A' allele.
# caseway Find the direction of GWAS effects between cases and controls After conducting an case-control association analysis on a SnpMatrix, e.g, GWAS using for example snp.rhs.tests from the snpStats package, it is not always trivial to determine the direction of the effects with respect to the reference allele. This function calculates which way around the data is for case vs control (pheno) and will indicate with respect to cases whether they have more copies of the reference allele, or less, or also can highlight whether the heterozygous is the affected genotype. This function works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which). Note that using a SnpMatrix, with het.effects=FALSE can be much faster (5-100x) than using a data.frame and/or setting het.effects=TRUE.
# randomize.missing Multicore randomised replacement of missing genotypes snpStats imputation only works if there are correlated SNPs with non-missing values that can be used to interpolate missing SNPs. If any correlated SNPs are missing 'impute.missing' will leave these blank. This function mops up the remainder by randomly inserting values consistent with the minor allele frequency of each SNP. This can be run using multiple cores to speed up progress for large matrices.
# impute.missing Replace missing values in a SnpMatrix object with imputed values This function is a wrapper for the snpStats snp.imputation() and impute.snps() functions. It allows a full imputation with one simple command, and facilitates stratified imputation for subsets of samples, and the parameter 'by' allows subsetting by SNPs, which speeds up the imputation for large sets by running it in smaller batches (as large ones are really slow). The standard use of the snpStats imputation functions will still leave some NA's behind, whereas the option 'random' ensures each missing value is replaced, even if just by an allele chosen at random (using existing frequencies). (chris)
# rSnpMatrix Create a SNP matrix with simulated data A simple function to simulation random SnpMatrix objects for testing purposes. Does not produce data with an 'LD' structure so is not very realistic.
# logistic.summary Function to produce a clean table from logistic regression done via GLM With input as a glm model result object, returns a clean dataframe with coefficients, p values, confidence intervals and standard errors. Multiple options for which columns to return and digits to display
# meta.me Meta-analysis using odds ratio and standard error from 2 datasets This function calculates meta analysis odds ratios, standard errors and p-values using results from a table containing odds ratio and standard error data for analyses of 2 different datasets (typically logistic regression, but other analyses can be incorporated if an odds-ratio and SE can be derived, for instance one analysis might be a case control logistic regression GWAS and the other a family TDT analysis).
# lambda_nm Normalize Lambda inflation factors to specific case-control count Lambda inflation statistics are influenced by the size of the generating datasets. To facilitate comparison to other studies, this function converts a given lambda from nr cases and mr controls, to n cases and m controls, where m=n=1000 is the most common normalization. All values other than 'Lnm' are forced within this function to be positive integers.
# abf Calculate approximate Bayes factors from p values and MAF This is a function to calculate approximate Bayes factors from p values and MAF - for reference see Wakefield, J (2009) Bayes factors for genome-wide association studies: comparison with p-values.
# read.pedData Import a ped file to pedData format used by snpStats PLINK ped files (family information files) are not in the same format as ped files required for use with snpStats, for instance for the tdt.snp() function to conduct a transmission disequilibrium test. This function will import a PLINK style ped file and return a 'pedData' object in the correct form for snpStats and other rpackages. The plink file format is: column 1: family id, column 2: individual id: column 3: father's ID or 0, column 4: mother's ID or 0, column 5: sex of subject, column 7: phenotype of subject.


### INTERNALS ####


#internal
# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}

finitize <- function(X) {
  if(is.data.frame(X)) { X <- as.matrix(X) }
  return(X[is.finite(X)])
}

minna <- function(...) {
  if(length(list(...))==1) { 
    min(finitize(...),na.rm=TRUE)
  } else {
    min(...,na.rm=TRUE)
  }
}

maxna <- function(...) {
  if(length(list(...))==1) { 
    max(finitize(...),na.rm=TRUE)
  } else {
    max(...,na.rm=TRUE)
  }
}

meanna <- function(...) {
  if(length(list(...))==1) { 
    mean(finitize(...),na.rm=TRUE)
  } else {
    mean(...,na.rm=TRUE)
  }
}

medianna <- function(...) {
  if(length(list(...))==1) { 
    median(finitize(...),na.rm=TRUE)
  } else {
    median(...,na.rm=TRUE)
  }
}

sdna <- function(...) {
  if(length(list(...))==1) { 
    sd(finitize(...),na.rm=TRUE)
  } else {
    sd(...,na.rm=TRUE)
  }
}

sumna <- function(...) {
  if(length(list(...))==1) { 
    sum(finitize(...),na.rm=TRUE)
  } else {
    sum(...,na.rm=TRUE)
  }
}

sortna <- function(...) {
  sort(..., na.last=TRUE)
}


#internal
comma <- function(...) {
  paste(...,collapse=",")
}


##' log sum function @author Claudia Giambartolomei - internal
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}



## functions from Chris W - adapted

# internal
# create a factor that can split a group of 'size' entries into categories size 'by'
# merge last two groups if final group size is less than min.pc of 'by'
# if fac.out is false, return start/end ranges rather than a grouping factor
# (Nick)
make.split <- function(size,by,fac.out=T,min.pc=0.5) {
  stepz <- round(seq(from=1,to=(size+1),by=by))
  if((tail(stepz,1)) != (size+1)) { stepz <- c(stepz,(size+1)) }
  split.to <- length(stepz)-1
  out1 <- cbind(stepz[1:split.to],(stepz[2:(split.to+1)]-1))
  repz <- 1+apply(out1,1,diff)
  out <- rep(1:split.to,repz)
  lrg <- max(out)
  if(!length(which(out==lrg))>(by*min.pc)) {
    out[out==lrg] <- lrg-1
  }
  if(!fac.out) {
    out <- cbind(start(Rle(out)),end(Rle(out)))
  }
  return(out)
}



# internal
# snpStats imputation only works if there are correlated SNPs with non-missing values
# that can be used to interpolate missing SNPs. If any correlated SNPs are missing
# 'impute.missing' will leave these blank. This function mops up the remainder
# by randomly inserting values consistent with the minor allele frequency of each SNP
# (Nick)
randomize.missing2 <- function(X,verbose=FALSE) {
  miss1 <- function(x) { 
    TX <- table(c(round(x),0,1,2))-c(1,1,1) # to force zero counts to be in the table
    naz <- which(is.na(x))
    if(length(naz)>0 & length(TX)>0) {
      x[naz] <- sample(as.numeric(names(TX)),size=length(naz),
                       replace=T,prob=as.numeric(TX))
    }
    return(x)
  }
  # randomly generate replacements for missing values using current distribution for each column of X
  if(is.null(dim(X))) { warning("not a matrix/data.frame") ; return(X) }
  count.miss <- function(x) { length(which(is.na(x))) }
  nmiss <- apply(X,2,count.miss)
  FF <- nrow(X)
  select <- nmiss>0
  if(length(which(select))<1) { return(X) }
  if(verbose) { cat(sum(nmiss),"missing values replaced with random alleles\n") }
  if(length(which(select))==1) { X[,select] <- miss1(X[,select]); return(X) }
  X[,select] <- apply(X[,select],2,miss1)
  return(X)
}


# internal, function generate a random MAF, then random SNP
rsnp <- function(n,A.freq.fun=runif,cr=.95, A.freq=NA) { 
  if(is.na(A.freq)) {  m <- A.freq.fun(1) } else { m <- A.freq }
  x <- sample(0:3, n, replace=T, prob=c(1-cr,cr*(m^2),cr*(2*(1-m)*m),cr*((1-m)^2)))
  return(x) 
}


# internal
rsnpid <- function(n) { 
  id.len <- sample(c(3:8),n,replace=T,prob=c(0.01, 0.01, 0.01, 0.10, 0.50, 0.37))
  each.id <- function(l) { sapply(l,function(n) { paste(replicate(n,sample(1:9,1)),collapse="",sep="") }) }
  sufz <- each.id(id.len)
  ids <- paste0("rs",sufz)
  return(ids)
}

# internal
ldfun <- function(n) { x <- runif(n); r2 <- runif(n); x[x<.6 & r2>.1] <- x[x<.6 & r2>.1]+.4 ; return(x) }

# internal
rsampid <- function(n,pref="ID0") { paste0(pref,pad.left(1:n,"0")) }

# internal - convert continuous values (by rank) to SNP-like data, coded 0:3
snpify.cont <- function(x,call.rate=.95,A.freq.fun=runif) { 
  if(length(x)<2) { stop("x must be longer than 1 element") }
  if(length(x)<5) { warning("for good simulation x should be fairly large, e.g, at least 10, better >100") }
  n <- length(x)
  if(all(x %in% 0:3)) { return(x) } # these are already snp-coded
  rr <- rank(x)
  sim1 <- rsnp(n,A.freq.fun=A.freq.fun,cr=call.rate) # A.freq=NA
  x[rr] <- sort(sim1)
  return(x)
}

# internal - get largest correlate(s) in an LD R square matrix
get.biggest <- function(r2.mat) {
  mm <- max(r2.mat,na.rm=T)
  coord <- which(r2.mat==mm,arr.ind=T)
  if(!is.null(dim(coord))){ coord <- coord[1,] }
  return(coord)
}

# internal - experimental function, was trying to 
# create SnpMatrices randomly with LD, and this
# was an attempt to extract the most intercorrelated
# SNPs from a larger set
get.top.n <- function(mat,n=10) {
  cr <- cor(mat,use="pairwise.complete")^2
  diag(cr) <- 0
  nn <- NULL
  while(length(nn)<n) { 
    coord <- get.biggest(r2.mat=cr)
    nn <- unique(c(nn,coord))
    cr[coord[1],coord[2]] <- NA
    cr[coord[2],coord[1]] <- NA
    #cr[,coord[1]] <- NA; cr[coord[1],] <- NA
    #cr[,coord[2]] <- NA; cr[coord[2],] <- NA
  }
  nn <- nn[1:n]
  new.mat <- mat[,nn]
  return(new.mat)
}

# internal
# allows an sapply style function to only work on valid values
clean.fn <- function(x,fail=NA,fn=function(x) { x }) {
  if(!is.null(x)) { 
    x <- x[!is.na(x)]; x <- x[(x!="")]; return(fn(x)) 
  } else {  return(fail) } 
}



#Internal: Read in a plink formatted pedigree/family file
# This function will import a PLINK style
# ped file and return a data.frame object in the same form
read.ped.file <- function(fn,keepsix=TRUE) {
  rr1 <- reader(fn,header=TRUE)
  if(ncol(rr1)<6) { warning("invalid ped/fam file, should have at least 6 columns"); return(NULL) }
  if(any(colnames(rr1) %in% c("X0","X1","X2"))) { rr1 <- reader(fn,header=FALSE) }
  colnames(rr1) <- gsub("X.","",colnames(rr1))
  if(any(colnames(rr1)[1] %in% unique(rr1[,1]))) { rr1 <- reader(fn,header=FALSE) }
  colnames(rr1)[1:6] <- c("family","sample","father","mother","sex","phenotype")
  if(keepsix) { rr1 <- rr1[,1:6] }
  return(rr1)
}




# internal function
validate.dir.for <- function(dir,elements,warn=F) {
  # in case the 'dir' input list object is not the standardised list form, convert
  # allows flexible use of list or regular directory specifications in plumbCNV functions
  if(is.null(dir)) { cat("directory empty\n"); return(NULL) }
  if(!is.list(dir)) {
    if(warn) { cat(elements[cc],"'dir' object wasn't a list\n")}
    dir <- as.list(dir); names(dir)[1:length(dir)] <- elements[1:length(dir)] 
  }
  for (cc in 1:length(elements)) {
    if(!elements[cc] %in% names(dir)) { 
      dir[[paste(elements[cc])]] <- "" ;
      if(warn) { stop(paste("dir$",elements[cc]," was empty.. set to current\n",sep="")) } 
    }
  }
  return(dir)
}

## internal function for the 'lambdas' function below
get.allele.counts <- function(myData,cc1000=FALSE) {
  ii <-  col.summary(myData)
  ii[["majmin"]] <- c("minor","major")[as.numeric(round(ii$RAF,3)!=round(ii$MAF,3))+1]
  if(cc1000) { ii$Calls <- rep(1000,nrow(ii)) }
  aa <- aA <- AA <- rep(0,nrow(ii))
  aa[which(ii$majmin=="minor")] <- (ii$P.BB*ii$Calls)[which(ii$majmin=="minor")]
  aa[which(ii$majmin=="major")] <- (ii$P.AA*ii$Calls)[which(ii$majmin=="major")]
  AA[which(ii$majmin=="minor")] <- (ii$P.AA*ii$Calls)[which(ii$majmin=="minor")] 
  AA[which(ii$majmin=="major")] <- (ii$P.BB*ii$Calls)[which(ii$majmin=="major")]
  aA <- ii$P.AB*ii$Calls
  ii[["aa"]] <- aa
  ii[["aA"]] <- aA
  ii[["AA"]] <- AA
  colnames(ii)[1] <- "TOTAL"
  return(ii[c("aa","aA","AA","TOTAL")])
}


# internal functions for lambdas #
Y_2 <- function(r1,r2,n1,n2,N,R) { (N*((N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((R*(N-R))*((N*(n1+(4*n2)))-(n1+(2*n2))^2)) }
X_2 <- function(r1,r2,n1,n2,N,R) { (2*N*((2*N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((4*R*(N-R))*((2*N*(n1+(2*n2)))-(n1+(2*n2))^2)) }
Likelihood_Lj <- function(c,Lnm) { rchisq(c/Lnm,df=1)/Lnm } # likelihood for one marker
LLikelihood_L <- function(Cj,LNMj) {
  # total likelihood across all K markers  : http://www.nature.com/ng/journal/v36/n4/full/ng1333.html
  tot <- 0
  for (cc in 1:length(Cj)) { 
    tot <- tot + log(Likelihood_Lj(Cj[cc],LNMj[cc])) 
  }
  return(tot) 
}



#internal
#' Choose which SNPs with duplicate locations to keep
#' 
#' Sometimes a microarray will have two or more instances of a SNP that
#' turns out to be at the same location, chromosome and position. It is
#' often desirable to remove duplicates, but to keep the one with the best
#' quality control measurements. This function streamlines this procedure.
#' @param snp.info snp support object for the microarray
#' @param snp.excl character vector, a list of snp ids already excluded for QC
#' @param ind.fn a datafile containing a list of duplicate pairs on the chip
#' @param col.name the name of the column to maximise in the kept duplicates
# not exported
choose.dups.to.drop <- function(snp.info, snp.excl=NULL, ind.fn="../iChipDupPairs.RData", col.name="call.rate") {
  dup.mat <- reader(ind.fn)
  if(!is.null(snp.excl)) {
    dup.mat <- apply(dup.mat,1,function(x) { x[which(x %in% snp.excl)] <- NA; return(x) })
  }
  dup.mat <- narm(t(dup.mat))
  
  cal.mat <- matrix(numeric(),nrow=nrow(dup.mat),ncol=ncol(dup.mat))
  if(!any(colnames(snp.info) %in% col.name)){ stop("column name ",col.name," was not found in snp.info") }
  ind1 <- match(dup.mat[,1],rownames(snp.info))
  ind2 <- match(dup.mat[,2],rownames(snp.info))
  notfound <- c((dup.mat[,1][which(is.na(ind1))]),(dup.mat[,2][which(is.na(ind2))]))
  if(any(is.na(ind1)) | any(is.na(ind2))) { stop("some snps: ",comma(notfound)," in the duplicate matrix were not in the snp.info, add 'snp.excl' entries for such") }
  cal.mat[,1] <- snp.info[ind1,][[col.name]]
  cal.mat[,2] <- snp.info[ind2,][[col.name]]
  j <- apply(cal.mat,1,function(x) { head(which(x==min(x,na.rm=T)),1) })
  i <- 1:nrow(cal.mat)
  DUPPOS.EXCL <- dup.mat[cbind(i,j)]
  return(DUPPOS.EXCL)
}


# internal function used by sync.asnp.mats
#' Add blank SNPs to SnpMatrix
add.missing.snps <- function(X, snp.info) {
  if(!is.null(dim(snp.info))) {
    all.snps <- rownames(snp.info)
    if(is.null(all.snps)) { stop("snp.info must be a snp.info/snp.support object") }
  } else {
    stop("snp.info must be a snp.info/snp.support object")
  }
  to.add <- all.snps[!all.snps %in% colnames(X)]
  add.m <- matrix(raw(),nrow=nrow(X),ncol=length(to.add))
  rownames(add.m) <- rownames(X)
  colnames(add.m) <- to.add
  add.sm <- as(add.m,"SnpMatrix")
  add.asm <- new("aSnpMatrix",.Data=add.sm, snps=snp.info[to.add,],
                 samples=X@samples,phenotype=X@phenotype,alleles=X@alleles)
  m2 <- cbind2(X,add.asm)
  m2 <- m2[,all.snps] #reorder
  return(m2)
}


#internal
## for an aSnpMatrix, extract all snps from a given chromosome
get.chr <- function(x,chr) {
  select <- which(x@snps$chromosome==chr)
  if(length(select)>0) {
    return(x[,select])
  } else {
    return(NULL)
  }
}

# internal
## from a list of SnpMatrix objects with all SNPs but different samples,
# extract all samples for a given chromosome
Get.chr <- function(X,CHR) {
  print(CHR) # abysmal memory usage - big fail!
  each.chr <- lapply(X,get.chr,chr=CHR)
  comb.chr <- do.call("rbind3",args=each.chr)
  return(comb.chr)
}



#internal
# get the column means, using just a regular
# mean if x happens to be a vector rather
# than a matrix
colmeansel <- function(x,sel) { 
  if(is.null(dim(x))) {
    mean(x[sel]) 
  } else {
    colMeans(x[sel,])
  }
}


## LD-Prune a BigSnpMatrix to reduce the size
## internal
ld.prune.big <- function(X,stats="R.squared",thresh=.1,n.cores=1) {
  if(estimate.memory(ncol(X)^2)<0.1) { ld.prune.chr(X,stats=stats,thresh=thresh) }
  #kk <- proc.time()[3]
  mat <- big.ld(X,stats=stats[1],filename=paste0("tempLdMat",sample(10^6,1)),limit.ram=F,n.cores=n.cores)
  #jj <- proc.time()[3]; print(jj-kk)
  #  exc <- length(matt[matt>thresh])
  nsnp <- ncol(mat)
  #  cat(out.of(exc,nsnp^2)," in matrix exceed ld threshold\n")
 # kk <- proc.time()[3];
  #  mat[is.nan(mat)] <- runif(length(which(is.nan(mat))))/100
  selected <- rep(T,ncol(mat))
  winners <- NULL
  all.sums <- colsum(mat,na.rm=T)
  all.ns <- nrow(mat) - colna(mat)
  while(any(selected)) {
    selectedPre <- selected
    tot.scores <- abs(all.sums/all.ns)
    tot.scores[is.na(tot.scores)] <- 0
    winnerz <- which(tot.scores==max(tot.scores[selected],na.rm=T))
    winnerz <- winnerz[winnerz %in% which(selected)]
    if(length(winnerz)<1) {
      winnerz <- which(tot.scores==max(tot.scores[selected],na.rm=T))
      winnerz <- winnerz[winnerz %in% winners]
      cat("this happened")
    }
    winner <- winnerz[sample(length(winnerz),1)]
    if(!is.na(winner) & (!winner %in% winners)) {
      winners <- c(winners,winner)
      passz <- (mat[winner,] < thresh) # here is where we'd need to add in a window if desired
      passz[is.na(passz)] <- T
      cat("winner ",colnames(mat)[winner],"had",length(which(!passz)),"in ld. ",out.of(length(which(selected)),length(selected)),"\n")
      selected <- selected & passz
      if(any(is.na(selected))) { prv(selected,passz)}
    }
    selected[winner] <- F
    ii <- which(selectedPre & !selected)
    segm <- mat[ii, ,drop=F]
    tot.dif <- colSums(segm,na.rm=T)
    n.dif <- apply(segm,2,function(x) { length(which(is.na(x)))})
    all.sums <- all.sums - tot.dif
    all.ns <- all.ns - nrow(segm) + n.dif
  }
#  jj <- proc.time()[3]; print(jj-kk)
  cat("pruned",out.of((nsnp-length(winners)),nsnp),"SNPs\n")
  return(winners)
}



# internal 
# take a set of strings and extract any text common to all, or only the unique text.
# e.g, unique text can be used to find differing pre-fixes or suffixes from the same
# root file name
string.in.common <- function(instrings, inverse=FALSE)
{
  ## takes a set of strings and produces
  ## a single string which is made up of the characters common to
  ## each string. can be useful in creating filenames or titles to 
  ## avoid repetition of long strings
  numstrz <- length(instrings)
  if(inverse) { newCharz <- vector("list",numstrz) } 
  for (cc in 1:numstrz)
  {
    nextCharz <- strsplit(instrings[cc],"")[[1]]
    if (cc==1)
    {
      curCharz <- nextCharz
    } else {
      curCharz <- curCharz[!is.na(match(curCharz,nextCharz))]
    }
  }
  if(!inverse) {
    outstring <- paste(curCharz,collapse="")
  } else {
    for (cc in 1:numstrz)
    {
      nextCharz <- strsplit(instrings[cc],"")[[1]]
      ii <- match(nextCharz,curCharz)
      adup <- which(duplicated(ii))
      if(length((adup))>0) {
        for (dd in 1:length((adup))) {
          ii[adup[dd]] <- match(nextCharz[adup[dd]],curCharz[-ii[adup[dd]]])
        }
      }
      newCharz[[cc]] <- nextCharz[is.na(ii)]
    }
    outstring <- lapply(newCharz, function(X) { paste(X,collapse="") })
  }
  return(outstring)
}


# internal
#' apply a function to each element of a SnpMatrixList
#' 
#' Applys a function to each subset of a SnpMatrixList. Also works for normal SnpMatrix,
#' aSnpMatrix or XSnpMatrix objects.
#' @param snpMatLst a SnpMatrixList 
#' @param fun the function to apply
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @param ... further arguments to fun()
#' @return by default returns a list, with the result for each outputted element
#' corresponding to the each input element in the same order. 
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' is(sml)
#' fun.snpMatLst(sml,dim)
fun.snpMatLst <- function(snpMatLst,fun=nrow,fail=T,dir=NULL,n.cores=1,...)
{
  # ... further arguments to fun()
  # get dimensions of snpMatLst (SnpMatrix list) regardless
  # of whether it's a set of SnpMatrix objects or list of file locations
  .do1 <- function(lst,dir,fun) {
    TF <- is.file(paste(lst),dir)
    if(TF) {
      fnm <- find.file(paste(lst),dir)
      snpMat <- get.SnpMatrix.in.file(fnm)
      #snpMat <- get(paste(load())) 
      res <- fun(snpMat,...)
    } else {
      warning(paste("invalid snpMat file",lst))
      res <- NULL
    }
    return(res)
  }
  HD <- switch(snp.mat.list.type(snpMatLst,fail),snpmatrix="snpmatrix",memory=F,disk=T,error=NULL)
  if(is.null(dir)) { dir <- getwd() } # ; warning("no directory passed to function") }
  if(HD=="snpmatrix") { return(fun(snpMatLst,...)) }
  if(HD) {
    n.grp <- length(snpMatLst)
    list.spec <- vector("list",n.grp)
    if(n.cores>1) {
      list.spec <- mclapply(snpMatLst,.do1,dir=dir,fun=fun,mc.cores=n.cores)
    } else {
      list.spec <- lapply(snpMatLst,.do1,dir=dir,fun=fun)
    }
  } else {
    if(n.cores>1) {
      list.spec <- mclapply(snpMatLst,fun,mc.cores=n.cores,...)
    } else {
      list.spec <- lapply(snpMatLst,fun,...)
    }
  }
  return(list.spec)
}


#(load("/ipswich/data/gdxbase/ld/ld_4.18/RDATA/ThousandGenomes_05_11_CEU/genotypes_chr22_ThousandGenomes_CEU_r0511.RData"))
#snps$genotypes,support
# if(F) {
#   chr.list <- vector("list",22)
#   for (cc in 1:22) {
#     (load(paste0("/ipswich/data/gdxbase/ld/ld_4.18/RDATA/ThousandGenomes_05_11_CEU/genotypes_chr",cc,"_ThousandGenomes_CEU_r0511.RData")))
#     keepers <- get.ichip.1000(snps$support)
#     SNPS <- snps$genotypes[,rownames(keepers)]
#     colnames(SNPS) <- keepers[["ichip"]]
#     chr.list[[cc]] <- SNPS
#   }
#   sml.1000 <- chr.list
# }


# internal or specific
get.ichip.1000 <- function(support,build=37) {
  sup <- df.to.ranged(support,chr="chromosome",start="position",end="position",build=build)
  keepers <- subsetByOverlaps(sup,chip.support(build=build))
  ii <- match(start(keepers),start(chip.support()))
  keepers[["ichip"]] <- rownames(chip.support())[ii]
  return(keepers)
}                                                                                                                                                                               



# internal
# returns a mean correction vector for a SnpMatrix
correct.by <- function(X) {
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix")) { stop("X must be a SnpMatrix object") }
  cs <- suppressWarnings(col.summary(X))
  if(nrow(cs)!=ncol(X)) { stop("col.summary length (",nrow(cs),") did not match length of X (",ncol(X),")")}
  raf <- cs$RAF
  #print(length(which(is.na(raf))))
  #print(cs[which(is.na(raf)),])
  nn <- cs$Calls
  aa <- raf^2
  ab <- raf*(1-raf)*2
  bb <- (raf)^2
  mN <- (aa*1) + (ab*2) + (bb*3)
  ss <- ((bb*nn)*(3-mN)^2) + ((ab*nn)*(2-mN)^2) + ((aa*nn)*(1-mN)^2)
  #  sD <- sqrt((1/(nn-1))*ss)  
  sD <- ss/nn
  sD[is.na(sD)] <- 1
  sD[sD<.0001] <- 0.1
  return(list(mean=mN,sd=sD))
}


# internal
# remove snps missing from either of 2 snp matrices
# if 'return.mat' is false, just returns the snpids to exclude.
# if true, then if return.x is true, returns the cleaned 'X' snpMatrix, else the clean 'Y'
rmv.common.missing <- function(X,Y,return.mat=FALSE,return.X=TRUE) {
  if(is.null(dim(X)) | is.null(dim(Y))) { stop("X and Y must have 2 dimensions (e.g, SnpMatrix) and have SNP names as column names")}
  x.snps <- colnames(X)
  y.snps <- colnames(Y)
  rafX <- (col.summary(X)$RAF)
  rafY <- (col.summary(Y)$RAF)
  print(length(which(is.na(rafX))))
  print(length(which(is.na(rafY))))
  x.bad <- x.snps[is.na(rafX)]
  y.bad <- y.snps[is.na(rafY)]
  more.bad1 <- x.snps[!x.snps %in% y.snps]
  more.bad2 <- y.snps[!y.snps %in% x.snps]
  all.bad <- unique(c(x.bad,y.bad,more.bad1,more.bad2))
  if(return.mat) {
    if(return.X) {
      return(X[,which(!clean.snp.ids(colnames(X)) %in% clean.snp.ids(all.bad))])
    } else {
      return(Y[,which(!clean.snp.ids(colnames(Y)) %in% clean.snp.ids(all.bad))])
    }
  } else {
    return(all.bad)
  }
}




# internal 
# for a snpMatLst, sync with snp.info and sample.info to exclude snps/samples not in the support files
sync.snpmat.with.info <- function(snpMatLst,snp.info=NULL,sample.info=NULL,dir=NULL,n.cores=1)
{
  # auto det2ect whether snpMatLst is a list of SnpMatrix objects or a list of RData
  # file locations and act accordingly. autodetect whether snp and/or sample info inputted.
  #print(length(snpMatLst)); print(is(snpMatLst)); print(dim(snp.info)); print(dim(sample.info)); print(head(dir))
  reorder.fn <- function(snpMatPart,info,snp=T) {
    if(snp) { nms <- colnames; txt <- "columns" } else { nms <- rownames; txt <- "rows" }
    to.keep <- match(rownames(info),nms(snpMatPart)) ; ll1 <- length(to.keep)
    #prv(to.keep)
    to.keep <- to.keep[!is.na(to.keep)] ; ll2 <- length(to.keep)
    #prv(to.keep,snpMatPart,info)
    if((1-(ll2/(ll1)))>.25) {
      warning(" discarding ",round((1-(ll2/(ll1)))*100,2),
              "% of ",txt, " that don't match the annotation file") 
    }
    if(snp) { return(snpMatPart[,to.keep]) } else { return(snpMatPart[to.keep,]) }
  }
  if(n.cores>1) { 
    require(parallel); l_apply <- function(...) { parallel::mclapply(...,mc.cores=n.cores) } 
  } else { 
    l_apply <- function(...) { lapply(...) }
  }
  if(!is.null(snp.info)) {
    if (is(snp.info)[1]=="RangedData" & is(snpMatLst[[1]])[1] %in% c("XSnpMatrix","SnpMatrix"))
    {
      cat(" reordering/selecting markers from SnpMatrices 1 to",length(snpMatLst),"...")
      snpMatLst <- l_apply(X=snpMatLst,FUN=reorder.fn,info=snp.info,snp=T)
      cat("done\n")
      return(snpMatLst)
    } else {
      if(is.ch(snpMatLst)) {
        if(all(is.file(unlist(snpMatLst),dir$cr,dir))) {
          for (cc in 1:length(snpMatLst)) {
            cat(" reordering/selecting markers from SnpMatrix",cc,"...")
            fn <- find.file(paste(snpMatLst[[cc]]),dir$cr,dir)
            snpMat <- get(paste(load(fn)))
            to.keep <- match(rownames(snp.info),colnames(snpMat))
            to.keep <- to.keep[!is.na(to.keep)]
            snpMat <- snpMat[,to.keep]
            save(snpMat,file=fn)
            cat("done\n")
          }
          return(snpMatLst)
        } else {
          warning("invalid input parameters, could not find listed file: ",unlist(snpMatLst))
        }
      } else {
        warning("invalid input parameters, need 'RangedData' or list of file names")
      }
    }
  }
  if(!is.null(sample.info)) {
    if (is(sample.info)[1]=="data.frame" & is(snpMatLst[[1]])[1] %in% c("XSnpMatrix","SnpMatrix"))
    {
      cat(" reordering/selecting samples from SnpMatrix files 1 to",length(snpMatLst),"...")
      snpMatLst <- l_apply(X=snpMatLst,FUN=reorder.fn,info=sample.info,snp=F)
      cat("done\n")
      return(snpMatLst)
    } else {
      if(is(snpMatLst[[1]])[1]=="character") {
        if(all(is.file(unlist(snpMatLst),dir$cr,dir))) {
          for (cc in 1:length(snpMatLst)) {
            cat(" reordering/selecting samples from SnpMatrix file",cc,"...")
            fn <- find.file(paste(snpMatLst[[cc]]),dir$cr,dir)
            snpMat <- get(paste(load(fn)))
            to.keep <- match(rownames(sample.info),rownames(snpMat))
            to.keep <- to.keep[!is.na(to.keep)]
            snpMat <- snpMat[to.keep,]
            save(snpMat,file=fn)
            cat("done\n")
          }
          return(snpMatLst)
        } else {
          warning("invalid input parameters, need 'RangedData' and list of SnpMatrix objects")
        }
      } else {
        warning("invalid input parameters, need 'RangedData' and list of SnpMatrix objects")
      }
    }
  } else {
    if(!is.null(snp.info)) { warning("neither snp.info or sample.info synced to snpMatrix")  }
  }
}



# internal, make a sample.info object using either just names, or using a snpMatLst too
sample.info.from.sml <- function(snpMatLst, subIDs.actual=NULL, dir=NULL, plink=FALSE) {
  if(!is.null(dir)) {
    if(snp.mat.list.type(snpMatLst)=="disk") {
      snpMatLst <- lapply(snpMatLst,find.file,dir=dir)
    }
  } else {
    dir <- getwd()
  }
  if(is.null(subIDs.actual)) {
    if(!is.null(snpMatLst)) {
      subIDs.actual <- rownamesL(snpMatLst,list=FALSE,dir=dir)
    } else {
      stop("Must enter at least 'subIDs.actual' (list of sample ids), or snpMatLst (sample ids as rownames)")
    }
  }
  ## get sample call rates ##
  if(is.character(subIDs.actual)) {
    if(length(subIDs.actual)>0) { names(subIDs.actual) <- NULL } } # prevent names convolving with contents
  if(is(snpMatLst[[1]])[1] %in% c("XSnpMatrix","SnpMatrix","aSnpMatrix")) {
    nrT <- nrowL(snpMatLst,list=FALSE); nrS <- nrowL(snpMatLst,list=TRUE)
    if(nrT==sum(unlist(nrS))) {
      group.nums <- rep(c(1:length(snpMatLst)),nrS)
    } else {
      group.nums <- rep(1,length(nrT))
    }
    if(length(group.nums)!=length(subIDs.actual)) {
      # print(table(group.nums)); print(length(group.nums)); print(sapply(snpMatLst,dim)[1,])
      stop("Error: mismatch between snpMat objects and number of subject ids in 'subIDs.actual'\n")
    }
  } else {
    if(plink) {
      group.nums <- rep(1,times=length(subIDs.actual))
    } else {
      lsp <- get.snpMat.spec(snpMatLst,dir=dir) #; print(dim(lsp)); print(snpMatLst)
      if(!(min(dim(lsp))>1 & length(dim(lsp))>0 )) { stop("Error: 'get.snpMat.spec' did not return 2D array\n")}
      if(sum(lsp[1,])==length(subIDs.actual)) {
        group.nums <- rep(c(1:length(snpMatLst)),times=lsp[1,])
      } else {
        group.nums <- rep(1,times=length(subIDs.actual))
      }
    }
  }
  if(length(group.nums)!=length(subIDs.actual)) { stop(paste("Error: group was length",length(group.nums),
                                                             "but IDs had length",length(subIDs.actual),"so couldn't construct a data.frame"))}
  sample.info <- data.frame(grp=group.nums,row.names=subIDs.actual,QCfail=rep(0,times=length(group.nums)))
  return(sample.info)
}



# internal
# used to apply the thresholds at the last step of doSnpQC()
apply.snp.thresholds <- function(snp.info,callrate.snp.thr=0.95,hwe.thr=10^-5,
                                 grp.cr.thr=.001,grp.hwe.z.thr=4,maf=NA,proc=1) {
  ## now remove samples failing HWE + callrate 95 regardless of source
  # also optionally remove those failing group criteria if columns present in snp.info
  snp.info[["QCfail"]] <- rep(0,times=nrow(snp.info))
  snp.info$call.rate[is.na(snp.info$call.rate)] <- 0
  snp.info$P.hwe[is.na(snp.info$P.hwe)] <- 0
  snp.info$Z.hwe[is.na(snp.info$Z.hwe)] <- 0
  snp.info$maf[is.na(snp.info$maf)] <- 0
  cr.cond <- snp.info$call.rate < callrate.snp.thr
  call.rate.excl.snps <- row.names(snp.info)[cr.cond]
  HWE.cond <- snp.info$P.hwe < hwe.thr
  HWE.exclude.snps <- row.names(snp.info)[HWE.cond]
  if(!is.na(maf)) {
    if(maf==0) { 
      MAF.cond <- snp.info$maf == 0 
    } else {
      MAF.cond <- snp.info$maf < maf
    }
    MAF.excl.snps <- row.names(snp.info)[MAF.cond]
  } else { MAF.cond <- rep(F,nrow(snp.info)); MAF.excl.snps <- NULL }
  to.remove <- unique(c(call.rate.excl.snps,HWE.exclude.snps,MAF.excl.snps))
  out.list <- list(snp.info,call.rate.excl.snps,HWE.exclude.snps,MAF.excl.snps,cr.cond,HWE.cond,MAF.cond)
  names(out.list) <- c("SNP.INFO","CR.EXCL","HWE.EXCL","MAF.EXCL","CR.cond","HWE.cond","MAF.cond")
  ## now also do the group level stats if present ##
  if(all(c("grp.hwe.zmax","grp.miss.p") %in% colnames(snp.info))) {
    snp.info$grp.miss.p[is.na(snp.info$grp.miss.p)] <- 1
    snp.info$grp.hwe.zmax[is.na(snp.info$grp.hwe.zmax)] <- 0
    grp.cr.cond <- snp.info$grp.miss.p < grp.cr.thr
    grp.cr.exclude.snps <- row.names(snp.info)[grp.cr.cond]
    grp.hwe.cond <- snp.info$grp.hwe.zmax > grp.hwe.z.thr
    grp.hwe.exclude.snps <- row.names(snp.info)[grp.hwe.cond]
    to.remove <- unique(c(to.remove,grp.cr.exclude.snps,grp.hwe.exclude.snps))
    out.list <- list(snp.info,call.rate.excl.snps,HWE.exclude.snps,grp.cr.exclude.snps,
                     grp.hwe.exclude.snps,maf.excl.snps,cr.cond,HWE.cond,grp.cr.cond,grp.hwe.cond,MAF.cond)
    names(out.list) <- c("SNP.INFO","CR.EXCL","HWE.EXCL","GRP.CR.EXCL","GRP.HWE.EXCL","MAF.EXCL","CR.cond",
                         "HWE.cond","grp.CR.cond","grp.HWE.cond","MAF.cond")
  }
  if(length(to.remove)>0)  {
    ## mark those failing on any threshold
    snp.info[["QCfail"]][match(to.remove, rownames(snp.info))] <- proc  }
  snp.info$QCfail[is.na(snp.info$QCfail)] <- 1
  out.list$SNP.INFO <- snp.info
  return(out.list)
}



# internal

# detect the format of a snpMatLst (SnpMatrix list)
# returns snpSubGroups, sampleSubGroups, or singleEntry (else error)
#' @param allow.dif.dims logical, whether to allow lists with completely
#' different dimensions for samples and snps (only recommended for sample-wise purposes)
get.split.type <- function(snpMatLst,head.only=TRUE,dir=NULL,allow.dif.dims=FALSE) {
  if(is.null(dir)) { dir <- getwd() }
  if(head.only) { snpMatLst <- head(snpMatLst) } # reduces time to return type for long lists
  if(!is.list(snpMatLst)) {
    if(is(snpMatLst)[1] %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")) {
      return("snpmatrix")
    } else {
      stop("Error: not a valid snpMatLst (not even a list)") 
    }
  }
  typz <- snp.mat.list.type(snpMatLst,fail=TRUE)
  HD <- switch(typz,disk=TRUE,memory=FALSE,error=NA)
  #  typz <- sapply(lapply(snpMatLst,is),"[",1)
  dir <- validate.dir.for(dir,c("cr"))
  #  if(all(typz %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")))
  #  { HD <- F } else {
  #    if (all(typz=="character")) { HD <- T } else {
  #      stop("snpMatLst doesn't appear to be a list of SnpMatrix objects",
  #           "or a list of file locations, row.summary impossible")
  #    } 
  #  }
  dimz <- get.snpMat.spec(snpMatLst,dir=dir)
  #prv(dimz)
  if(ncol(dimz)==1) { return("singleEntry") }
  # get DIMs of separate SNP lists (e.g, chromosomes)
  zroSNP <- as.numeric(names(table(diff(dimz[2,]))))
  zroSMP <- as.numeric(names(table(diff(dimz[1,]))))
  #prv(dimz,zroSNP,zroSMP)
  same.num.snps <- same.num.samples <- TRUE
  if(length(zroSNP)>0) { 
    if (any(zroSNP!=0)) { 
      same.num.snps <- FALSE
    } 
  } else {
    stop("empty dims list") 
  }
  if(length(zroSMP)>0) { 
    if (any(zroSMP!=0)) { 
      same.num.samples <- FALSE
    } 
  } else {
    stop("empty dims list") 
  }
  if(same.num.snps & same.num.samples) {
    ## need to investigate further
    rows.equal <- cols.equal <- FALSE
    dn <- fun.snpMatLst(snpMatLst,fun=dimnames,fail=T,dir=dir)
    LL <- length(dn)
    rn.first <- dn[[1]][[1]]; rn.last <- dn[[LL]][[1]]
    cn.first <- dn[[1]][[2]]; cn.last <- dn[[LL]][[2]] 
    if(all(rn.first==rn.last)) { rows.equal <- TRUE  }
    if(all(cn.first==cn.last)) { cols.equal <- TRUE  }
    if(rows.equal & !cols.equal) { return("snpSubGroups") }
    if(!rows.equal & cols.equal) { return("sampleSubGroups") }
    if(!rows.equal & !cols.equal) { stop("inconsistent ids in snpMatLst object") }
    if(rows.equal & cols.equal) {
      if(length(snpMatLst)>1) { 
        stop("duplicate entries in snpMatLst") 
      } else {
        return("singleEntry")
      }
    }                
  } else {
    if(same.num.snps & !same.num.samples) { return("sampleSubGroups") }
    if(!same.num.snps & same.num.samples) { return("snpSubGroups") }
    if(!same.num.snps & !same.num.samples & !allow.dif.dims) {
      stop("invalid/incomplete snpMatLst object") 
    }
  }
}



# internal
# roughly random vector to create a grouping variable for
# n elements into 'grps' roughly equally sized groups.
split.vec <- function(n,grps=10) {
  init <- as.integer(cut(runif(n),breaks=sort(runif(grps))))
  init[is.na(init)] <- 0
  init <- sort(init+1)
  return(init)
}



# internal
# detect whether object is aSnpMatrix type, even if within a snpMatLst
#' @param trust.first logical, if TRUE, will only test the first two
#' files (if a list) in order to save time
is.aSnpMatrix <- function(X,any=FALSE,dir=NULL,trust.first=FALSE) {
  if(is(X)[1]=="aSnpMatrix") { return(TRUE) }
  if(is.list(X)) {
    if(trust.first) { X <- head(X,2) } # reduce time taken
    typ <- suppressWarnings(snp.mat.list.type(X,fail=F))
    if(typ=="error") { return(FALSE) }
    if(typ=="disk") { 
      mis <- smlapply(X,FUN=function(x) { is(x)[1] },c.fun=unlist,dir=dir) 
    } else {
      mis <- sapply(lapply(X,is),"[",1)
    }
    if(any) {
      if(any(mis %in% c("aXSnpMatrix","aSnpMatrix"))) { return(TRUE) } else { return(FALSE) }
    } else {
      if(all(mis %in% c("aSnpMatrix","aXSnpMatrix"))) { return(TRUE) } else { return(FALSE) }
    }
  } else {
    return(FALSE)
  }
}



# internal
# Create a [partial] sample support data.frame from an aSnpMatrix object
samp.from.annot <- function(aSnpMat) {
  if(!is(aSnpMat)[1]=="aSnpMatrix") { stop("aSnpMat must be of type 'aSnpMatrix'") }
  df <- aSnpMat@samples
  pheno.col <- which(colnames(df)==aSnpMat@phenotype)
  sex.col <- c(grep("sex",tolower(colnames(df))),grep("gender",tolower(colnames(df))))[1]
  if(is.na(sex.col)) { sex.col <- NULL }
  df <- df[,c(pheno.col,pheno.col,sex.col,pheno.col)]
  colnames(df)[1:2] <- c("grp","phenotype")
  if(ncol(df)==4) { colnames(df)[3:4] <- c("sex","QCfail") } else { colnames(df)[3] <- "QCfail" }
  df$grp <- rep(1,nrow(df)); df$QCfail <- rep(0,nrow(df))
  return(df)
}


##### MAIN FUNCTIONS ######


#' Synchronise SNPs amongst slightly differing datasets
#' 
#' When conducting a study (e.g, meta-analysis) and you are pooling
#' genotype datasets you may find discrepancies in the SNPs included in each
#' even if taken from the same platform, for example the study may only
#' make available post-QC datasets, and so the different datasets all have
#' slightly different subsets of the same overall SNP-set. This function 
#' aids such analyses by adding in blank entries for each dataset for any
#' SNPs found in other groups that are missing, to the effect that
#' each dataset will have the same SNP IDs.
#' @param myMat a list of SnpMatrix objects with different samples but mostly the same SNPs,
#'  but different omitted, this helps synchronise them by inserting missing vectors of snps 
#'  not missing in at least one dataset into the datasets where they are not present
#' @return a sample-subset SnpMatrixList where blanks are inserted to synchronise SNP labels for
#' each comprising SnpMatrix object
#' @examples
#' # too hard to simulate?
#' #' snp.set <- sample(rownames(chip.support()),20) # select a random set of 10 SNPs
#' data(exAnnotSnp)
#' sm1 <- exAnnotSnp[1:10,sample(20,14)]
#' sm2 <- exAnnotSnp[11:30,sample(20,14)]
#' sm3 <- exAnnotSnp[31:50,sample(20,14)]
#' print(cbind(colnames(sm1),colnames(sm2),colnames(sm3))) # colnames mismatch
#' sml <- sync.asnp.mats(sm1,sm2,sm3) # create a SnpMatrixList
#' print(cbind(colnames(sml[[1]]),colnames(sml[[2]]),colnames(sml[[3]]))) # colnames match
sync.asnp.mats <- function(m1,m2,...) {
  myMat <- list(m1,m2)
  extr <- list(...)
  if(length(extr)>0) {
    myMat <- c(myMat,extr)
  }
  typz <- sapply(myMat,function(X) { is(X)[1] })
  if(!all(typz %in% c("aSnpMatrix","aXSnpMatrix"))) {
    stop("all arguments must be aSnpMatrix objects") }
  nn <- length(myMat)
  sMat <- vector("list",nn)
  for (cc in 1:nn) {
    sMat[[cc]] <- myMat[[cc]]@snps
  }
  snpso <- do.call(rbind,args=sMat)
  snpso <- snpso[!duplicated(snpso$snp.name),]
  for (cc in 1:nn) {
    myMat[[cc]] <- add.missing.snps(myMat[[cc]],snpso)
  }
  return(myMat)
}



#' Flip the strand for allele codes
#' 
#' When strands are being converted sometimes strand flipping is necessary,
#' whereby C goes to G, A goes to T, etc. This function takes care of this
#' in the most intuitive way possible, and can also handle IUPAC ambiguity 
#' codes, and separators within the text.
#' @param acgt character, allele codes, a series of letters, may include
#' IUPAC ambiguity codes, "R","Y","M","K","S","W","H","B","V","D","N", in
#' addition to the standard A, C, G, T codes. Can also involve a separator
#' in the text, e.g, C/T, A/T, etc, as long as the seperator being used
#' is passed to the 'sep' parameter.
#' @param sep optional separator (ignore if entering single letters), that
#' goes between pairs of alleles. The function will work if more that 2
#' alleles are used as long as multiple separators are used
#' @export
#' @return returns a result in the same character format as that entered,
#' including the same separator if any, but with the strands flipped
#' @examples
#' flip.strand(c("A","C","G","M"))
#' flip.strand(c("A/C","C/T","G/A","M/B"))
#' flip.strand(c("A|C|T|G","C|T","G|A|C","M|B"),sep="|")
flip.strand <- function(acgt,sep="/") {
  new <- paste(acgt)
  nn <- length(acgt)
  sepz <- grep(sep,acgt,fixed=TRUE)
  if(length(sepz)> 0.5*nn) {
    two.sets <- strsplit(acgt,split=sep,fixed=TRUE)
    lenz <- sapply(two.sets, length)
    #if(any(lenz>2)) {
    # warning("A maximum of 2 allele codes should be entered with a separator, e.g, A/C",
    #          ". If more types are possible, use IUPAC ambiguity codes") }
    result <- sapply(lapply(two.sets, flip.strand),paste,sep="",collapse=sep)
    return(result)
  } else {
    if(length(which(nchar(acgt)>1)) > 0.5*nn) { 
      warning("allele codes should only be 1 character long, the majority of 'acgt' entries were longer than this")
    } else {
      ## valid, and continue
    }
  }
  if(is(new)[1]!="character") { stop("invalid input, needs to be coercible to characters") }
  # define recoding characters, including uncertain codes #
  pre <- c("A","T","G","C","R","Y","M","K","S","W","H","B","V","D","N")
  post <- c("T","A","C","G","Y","R","K","M","S","W","D","V","B","H","N")
  pre <- c(toupper(pre),tolower(pre)) # create upper and lower case versions
  post <- c(toupper(post),tolower(post)) # create upper and lower case versions
  if(!any(acgt %in% pre)) { warning("no characters in 'acgt' were valid allele codes") }
  if(length(pre)!=length(post)) { stop("internal function error, lists invalid, please report bug") }
  for (cc in 1:length(pre)) {
    new[acgt==pre[cc]] <- post[cc]
  }
  return(new)
}

#' Bind column-wise for more than 2 SnpMatrix objects
#' 
#' This is a cbind function for more than 2 [a]SnpMatrix objects at once
#' @param ... any number of SnpMatrix or aSnpMatrix objects with the same
#' set of samples, but different SNPs
#' @param silent logical, whether to display information on progress
#' @return a larger SnpMatrix with the same samples but more SNPs
#' @seealso rbind3
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' sm1 <- sml[[1]] # chromosome 1 SNPs from sml
#' sm2 <- sml[[2]] # chromosome 2 SNPs from sml
#' sm3 <- sml[[3]] # etc
#' sm1_3 <- cbind3(sm1,sm2,sm3) # yields an object with chr's 1 thru 3
cbind3 <- function(...,silent=TRUE) {
  objlist <- list(...)
  Dimm <- function(X) { paste(dim(X),collapse=",") }
  isz <- unlist(lapply(lapply(objlist,is),"[",1))
  if(all(isz=="SnpMatrix")) { return(cbind(...)) }
  if(all(isz=="XSnpMatrix")) { return(cbind(...)) }
  if(!all(isz=="aSnpMatrix")) { stop("all arguments must have type aSnpMatrix, from the annotSnpStats package") }
  lo <- length(objlist)
  if(lo==1) { return(...) }
  if(lo==2) { return(cbind2(...)) }
  if(!silent) { cat("object 1 has dim: ",Dimm(objlist[[1]]),"\nappending object 2/",lo," [dim:",Dimm(objlist[[2]]),"]\n",sep="") }
  outlist <- cbind2(objlist[[1]],objlist[[2]])
  for (cc in 3:lo) {
    if(!silent) { cat("appending object ",cc,"/",lo," [dim:",Dimm(objlist[[cc]]),"]\n",sep="") }
    outlist <- cbind2(outlist,objlist[[cc]])
  }
  if(!silent) { cat("produced combined object with dim:",Dimm(outlist),"\n") }
  return(outlist)
}


#' Bind row-wise for more than 2 SnpMatrixList objects
#' 
#' This is an rbind function for more than 2 SnpMatrixList objects at once, but
#' if 'mismatching.snps' are used, must all be aSnpMatrix based lists.
#' @param ... any number of SnpMatrixList objects with the same
#' set of SNPs [or overlapping], but different samples
#' @param mismatching.snps logical, whether to assume that the SNP columnnames
#' in each SnpMatrixList might be different (faster to set FALSE if you know they
#' are the same).
#' @param out.str character, root for the filenames of the resulting SnpMatrixList, can 
#' include a directory path.
#' @param dir directory (character) for the list if not using absolute paths
#' @return a larger SnpMatrix with the same SNPs but more samples
#' @seealso rbind3
#' @examples
#' # THIS FUNCTION IS CURRENTLY UNTESTED!!! #
rbind.sml <- function(snpMatLst,...,out.str="samplesCombined",mismatching.snps=TRUE,dir=NULL) {
  smlz <- list(...)
  smlz <- c(list(snpMatLst),smlz)
  ls <- length(smlz)
  if(ls==1) { return(snpMatLst) }
  ll <- sapply(smlz,length)
  dir <- dirname(out.str)
  fn <- basename(out.str)
  if(length(unique(ll))>1) { stop("snpMatLst and all ... lists must have the same length") }
  for (cc in 1:ls) {
    typ <- snp.mat.list.type(smlz[[cc]],fail=TRUE) #disk or memory
    typ2 <- get.split.type(smlz[[cc]],dir=dir)  # singleEntry snpSubGroups sampleSubGroups
    if(typ2=="snpmatrix") { stop("at least one argument was not a SnpMatrixList") }
    if(typ2=="sampleSubGroups") { stop("use cbind.sml for sample subgroups") }
    if(typ2!="snpSubGroups") { stop("at least one argument was an improper snpMatLst parameter") }    
  }
  #out <- vector("list",ll[1])
  SMLZ <- vector("list",ll[1])
  for (dd in 1:(ll[1])) {
    if(mismatching.snps) {
      argz <- lapply(smlz,function(X) { get.SnpMatrix.in.file(X[[dd]]) })
      prv(argz[[1]],argz[[2]],argz[[3]])
      obj <- do.call("sync.asnp.mats",args=argz)
    } else {
      obj <- lapply(smlz[[dd]],function(X) { get.SnpMatrix.in.file(X) })
    }
    obj2 <- do.call("rbind3",args=obj)
    if(typ=="disk") {
      fnm <- cat.path(dir,fn=fn,suf=dd,ext="RData")
      save(obj2,file=fnm)
      SMLZ[[dd]] <- fnm
    } else {
      SMLZ[[dd]] <- obj
    }
  }
  return(SMLZ)
}

  

#' Bind row-wise for more than 2 SnpMatrix objects
#' 
#' This is an rbind function for more than 2 [a]SnpMatrix objects at once
#' @param ... any number of SnpMatrix or aSnpMatrix objects with the same
#' set of SNPs, but different samples
#' @param silent logical, whether to display information on progress
#' @return a larger SnpMatrix with the same SNPs but more samples
#' @seealso cbind3
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sm1 <- snp.mat[1:10,]
#' sm2 <- snp.mat[11:20,]
#' sm3 <- snp.mat[21:30,]
#' sm4 <- snp.mat[31:40,]
#' smCMB <- rbind3(sm1,sm3,sm4) # yields an object with samples 1:10, 21:40
#' prv(smCMB)
rbind3 <- function(...,silent=TRUE){
  objlist <- list(...)
  Dimm <- function(X) { paste(Dim(X),collapse=",") }
  isz <- unlist(lapply(lapply(objlist,is),"[",1))
  if(all(isz=="SnpMatrix")) { return(rbind(...)) }
  if(all(isz=="XSnpMatrix")) { return(rbind(...)) }
  if(!all(isz=="aSnpMatrix")) { stop("all arguments must have type aSnpMatrix, from the annotSnpStats package") }
  lo <- length(objlist)
  if(lo==1) { return(...) }
  if(lo==2) { return(rbind2(...)) }
  if(!silent) { cat("object 1 has dim: ",Dimm(objlist[[1]]),"\nappending object 2/",lo," [dim:",Dimm(objlist[[2]]),"]\n",sep="") }
  outlist <- rbind2(objlist[[1]],objlist[[2]])
  for (cc in 3:lo) {
    if(!silent) { cat("appending object ",cc,"/",lo," [dim:",Dimm(objlist[[cc]]),"]\n",sep="") }
    outlist <- rbind2(outlist,objlist[[cc]])
  }
  if(!silent) { cat("produced combined object with dim:",Dimm(outlist),"\n") }
  return(outlist)
}


#' Show information on what type of SnpMatrixList is in use
#' 
#' Determine whether a snpMatLst is a list of filenames, SnpMatrix objects, 
#' or just a [a]SnpMatrix. SnpMatrixList objects can be a list stored in 'memory' or 
#' a list stored on 'disk'. This function returns a string to indicate which, or 'fail' 
#' if not a compatible object, or 'snpmatrix' if it is just a regular SnpMatrix.
#' @param snpMatLst a SnpMatrixList object
#' @param fail logical, whether to throw an error instead of a warning when 
#' snpMatLst is not a SnpMatrix or SnpMatrixList
#' @return character string, either "memory", "disk", "snpmatrix" or "fail".
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' snp.mat.list.type(sml)
#' data(exSnpMat)
#' snp.mat.list.type(exSnpMat)
#' data(exAnnotSnp)
#' snp.mat.list.type(exAnnotSnp)
#' snp.mat.list.type("nonsense")
snp.mat.list.type <- function(snpMatLst,fail=FALSE)
{
  # plumbCNV snp-qc can be run on a list of snpMatrix objects in memory ('memory') or
  # on a list of RData files on 'disk' to save RAM, this fn detects which type is in use
  if(!is.list(snpMatLst)) {
    if(is(snpMatLst)[1] %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")) {
      return("snpmatrix")
    } else { 
      HDt <- "error"
      if(fail) {  stop("Error: not a valid snpMatLst (not even a list)") } 
      return(HDt)
    }
  }
  typz <- sapply(lapply(snpMatLst,is),"[",1)
  if(all(typz %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")))
  { HDt <- "memory" } else {
    if (all(typz=="character")) { HDt <- "disk" } else {
      HDt <- "error"
      if(fail) {  stop("Error: not a valid snpMatLst!") 
      } else { warning("couldn't classify type as was not a valid snpMatLst")  }
    }
  }
  return(HDt)
}


#' Trim samples in a list of files to a common set
#' 
#' To make a SnpMatrixList from a set of files, the files
#' must have a common sample set, or a common SNP set. This function
#' helps to trim down to a common set if the sample/SNP-list differs
#' a little between files (e.g, from different sources, different exclusions, etc)
#' @param sml a list of file names containing [a]SnpMatrix objects (almost a SnpMatrixList)
#' but the requirement for a common SNP or sample set is not fulfilled.
#' @param trim default is 'common' which will look at all row/column names (depending)
#' on 'snps' and select the set common to all listed objects. Otherwise you can
#' select a custom set of SNPs/samples by inputting them here as a string vector
#' @param snps logical, whether to trim to a common/selected set of SNPs, or if
#' set to FALSE, then to set of samples
#' @export
#' @seealso snpSel, sampSel
#' @return return a SnpMatrixList trimmed appropriately
#' # file.set.trim(new.sml,dir="/chiswick/data/ncooper/imputation/NEW/")
file.set.trim <- function(sml,trim="common",snps=FALSE,dir=NULL,tracker=TRUE) {
  # even up number of samples #
  ll <- length(sml)
  if(!is.list(sml)) { stop("sml must be a list") }
  if(!all(sapply(sml,file.exists))) { stop("some elements of sml did not point to existing files") }
  if(all(trim=="common")) {
    rn <- lapply(sml,function(x) { i <- reader(x); return(if(!snps) {rownames(i) } else {colnames(i)}) })
    tt <- table(unlist(rn))
    my.rn <- names(tt)[tt==ll]
    if(length(my.rn)<1) { stop("There were no common",if(snps) {"SNPs"} else {"samples"},"across files in 'sml'. Try setting 'snps' to ",!snps,"\n") }
  } else {
    my.rn <- trim
  }
  fnz <- character(ll)
  if(is.character(dir)) { if(!file.exists(dir)) { dir <- NULL } }
  for (j in 1:ll) { 
    if(is.null(dir)) { dirn <- dirname(sml[[1]]) } else { dirn <- dir }
    X <- reader(sml[[j]]); 
    if(!snps) { 
      if(!any(rownames(X) %in% my.rn)) { 
        warning("no samples from 'trim' list found in file",j,"\n") 
        X <- X[0,]
      } else {
        X <- X[my.rn[my.rn %in% rownames(X)],] 
      } 
    } else { 
      if(!any(colnames(X) %in% my.rn)) { 
        warning("no SNPs from 'trim' list found in file",j,"\n") 
        X <- X[,0]
      } else {
        X <- X[,my.rn[my.rn %in% colnames(X)]] 
      }
    } 
    fnz[j] <- cat.path(dirn,basename(sml[[j]]))
    #prv(fnz[j])
    save(X,file=fnz[j]); gc() 
    if(tracker) { loop.tracker(j,ll) }
  }
  new.sml <- as.list(fnz)
  if(tracker) { cat("wrote SnpMatrixList files to",unique(dirname(fnz)),"\n") }
  return(new.sml)
}


#' Select all SNPs for a subset of samples in a snpMatLst
#' 
#' Subsetting method for a SnpMatrixList, but designed for
#' smaller subsets, and returns a SnpMatrix rather than a 
#' SnpMatrixList. 
#' @param snpMatLst a SnpMatrixList
#' @param character vector, a list of samples (indexes/labels), recommended less than 1,000, but depends on
#' on the number of SNPs
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param SnpMatrix logical, if TRUE then this function outputs a SnpMatrix object
#' containing the selected SNPs. If FALSE, outputs a SnpMatrixList in the same
#' location as the original unless 'dir' is provided
#' @return returns a SnpMatrix containing the specified SNP
#' subset in the specified order, and the full set of samples from snpMatLst.
#' The 'snps' selection must imply a dataset that can fit into a regular SnpMatrix.
#' @seealso sampSel
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' sub.samps <- rownamesL(sml)[c(2,5,9)]
#' subs <- sampSel(sml,sub.samps)
#' prv(subs)
sampSel <- function(snpMatLst,samples,dir=NULL,SnpMatrix=TRUE) {
  typ <- snp.mat.list.type(snpMatLst)
  if(typ=="error") { stop("snpMatLst was not a valid SnpMatrixList object") }
  if(!SnpMatrix) {
    return(file.set.trim(snpMatLst,trim=samples,snps=FALSE,dir=dir,tracker=FALSE))
  } 
  nc <- ncolL(snpMatLst,list=FALSE); nr <- length(samples)
  res <- check.snpmat.size(nr,nc)
  if(!res) { stop("selected SNP-set would result in a SnpMatrix too large to be created") }
  if(length(samples)>1000) { warning("This function is designed for small lists of 'samples', likely to fail for large datasets")}
  sampselfun <- function(snpMat,samples) {
    if(is.character(samples)) { samples <- narm(match(samples,rownames(snpMat))) }
    if(length(samples)>0) { return(snpMat[samples,]) } else { return(NULL) }
  }
  newlist <- fun.snpMatLst(snpMatLst,fun=sampselfun,dir=dir,samples=samples)
  nulz <- sapply(newlist,is.null)
  newlist <- newlist[!nulz] # remove NULLs from this list or else rbind doesn't work properly
  typ <- get.split.type(snpMatLst,dir=dir)
  if(typ=="snpSubGroups") {
    if(is(newlist[[1]])[1] %in% c("aSnpMatrix","aXSnpMatrix")) { fun <- "cbind3" } else { fun <- "cbind" }
  } else {
    if(is(newlist[[1]])[1] %in% c("aSnpMatrix","aXSnpMatrix")) { fun <- "rbind3" } else { fun <- "rbind" }
  }
  outobj <- do.call(fun,args=newlist)  
  if(is.character(samples)) { outobj <- outobj[samples,] }
  if(is.character(samples)) { outobj <- outobj[narm(match(samples,rownames(outobj))),] }
  return(outobj)
}

#' Select all samples for a subset of SNPs in a snpMatLst
#' 
#' Subsetting method for a SnpMatrixList, but designed for
#' smaller subsets, and returns a SnpMatrix rather than a 
#' SnpMatrixList. 
#' @param snpMatLst a SnpMatrixList
#' @param character vector, a list of snps (names/labels), recommended less than 50,000, but depends on
#' on the number of samples.
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param SnpMatrix logical, if TRUE then this function outputs a SnpMatrix object
#' containing the selected SNPs. If FALSE, outputs a SnpMatrixList in the same
#' location as the original unless 'dir' is provided
#' @return returns a SnpMatrix containing the specified SNP
#' subset in the specified order, and the full set of samples from snpMatLst.
#' The 'snps' selection must imply a dataset that can fit into a regular SnpMatrix.
#' @seealso sampSel
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' sub.snps <- snp.set[c(1,4)]
#' subs <- snpSel(sml,sub.snps)
#' prv(subs)
snpSel <- function(snpMatLst,snps,dir=NULL,SnpMatrix=TRUE) {
  typ <- snp.mat.list.type(snpMatLst)
  if(typ=="error") { stop("snpMatLst was not a valid SnpMatrixList object") }
  if(!SnpMatrix) {
    return(file.set.trim(snpMatLst,trim=snps,snps=TRUE,dir=dir,tracker=TRUE)) # set FALSE
  } 
  nc <- length(snps); nr <- nrowL(snpMatLst,list=FALSE)
  res <- check.snpmat.size(nr,nc)
  if(!res) { stop("selected SNP-set would result in a SnpMatrix too large to be created") }
  if(length(snps)>50000) { warning("This function is designed for small lists of 'snps', likely to fail for large datasets")}
  if(!is.character(snps)) { warning("This function is designed for character indexing, results may not be what you expect") }
  snpselfun <- function(snpMat,snps) {
    if(is.factor(snps)) { snps <- paste(snps) }
    if(is.character(snps)) { 
      snps <- clean.snp.ids(snps)
      snps <- narm(match(snps,clean.snp.ids(colnames(snpMat)))) 
    }
    if(is.numeric(snps)) {
      snps <- snps[snps %in% 1:ncol(snpMat)]
    } 
    if(length(snps)>0) { 
      #prv(snps,snpMat)
      #return(list(snps,snpMat))
      #print(length(which(is.na(snps)))); print(summary(snps))
      out <- snpMat[,snps]
      #prv(out)
      return(out) 
    } else { return(NULL) }
  }
  newlist <- fun.snpMatLst(snpMatLst,fun=snpselfun,dir=dir,snps=snps)
  nulz <- sapply(newlist,is.null)
  newlist <- newlist[!nulz] # remove NULLs from this list or else cbind doesn't work properly
 # prv(newlist)
  typ <- get.split.type(snpMatLst,dir=dir)
 
  if(length(newlist)>0) {  typ2 <- is(newlist[[1]])[1] } else { 
    typ2 <- ""; warning("result of selection seems empty") 
    return(NULL)
  }
  if(typ=="snpSubGroups") {
    if(typ2 %in% c("aSnpMatrix","aXSnpMatrix")) { fun <- "cbind3" } else { fun <- "cbind" }
  } else {
    if(typ2 %in% c("aSnpMatrix","aXSnpMatrix")) { fun <- "rbind3" } else { fun <- "rbind" }
  }
  #prv(newlist)
  if(length(newlist)>1) {
    outobj <- do.call(fun,args=newlist) 
  } else {
    outobj <- newlist[[1]]
  }
  colnames(outobj) <- clean.snp.ids(colnames(outobj))
  if(is.character(snps)) { snps <- clean.snp.ids(snps) ; outobj <- outobj[,narm(match(snps,colnames(outobj)))] }
  return(outobj)
}


# MY EXAMPLE of READING IMPUTE FILES ##
# if(F) {
# 
# WTCCC2.dir <- "/ipswich/data/chrisw/MS-sawcer/WTCCC2"
# subdir <- "CCC2_Cases_UKC/"
# fulldir <- cat.path(WTCCC2.dir,fn="",pref=subdir)
# dat.fn <- list.files(fulldir)
# dat.fn <- cat.path(fulldir,dat.fn[grep(".gen.gz",dat.fn)])
# samp.fn <- cat.path("/ipswich/data/chrisw/MS-sawcer/WTCCC2","MS_UKC_illumina.sample",pref=subdir)
# ms.dr <- "/chiswick/data/ncooper/imputation/MS/"
# sml <- read.impute.list(dat.fn=dat.fn, samp.fn=samp.fn, combine=TRUE,
#                             snpcol=1, out.pref="TEMP", out.dir=ms.dr, verbose=TRUE)
# }


#' Read impute file(s) into a SnpMatrix[List] format
#' 
#' IMPUTE2 is a popular imputation method. This function facilitates
#' reading the output into R to use with snpStats. When there are 
#' multiple files assumed to be same samples, different SNP sets (e.g, chromsomes).
#' A wrapper for the snpStats 'read.impute' function
#' @param dat.fn file name(s) of the impute files to read in, may be gzipped
#' @param sample.fn file name of the *.sample file
#' @param HD logical, whether to be a disk-type SnpMatrixList or a memory type
#' @param combine logical, whether to attempt to combine each chromosome into
#' a single large SnpMatrix object, or to leave as separate objects
#' @param snpcol passed to snpStats:read.impute, Which column of the input
#'  will be used as the SNP name. Default is column 2
#' @param out.pref character, prefix for output files
#' @param out.dir directory to store data files when HD=TRUE
#' @param verbose logical, TRUE to show extra information on progress
#' @return a SnpMatrixList (or SnpMatrix if combine=TRUE and object can fit in memory)
#' with the imported data
#' @export
#' @examples
#' ## not run ## (owndata) ?
read.impute.list <- function(dat.fn, samp.fn=NULL,
                            HD=TRUE, combine=FALSE, snpcol=2, out.pref="SnpMatrix", 
                            out.dir=getwd(), verbose=FALSE) {
  samp <- reader(samp.fn)
  if(length(rownames(samp)[1])==1) { samp <- samp[-1,] }
  samp.ids <- rownames(samp)
  ii <- which(tolower(colnames(samp)) %in% "sex")
  if(length(ii)>0) { sex <- samp[[ii]] } else { sex <- NULL }
  ii <- which(substr(tolower(colnames(samp)),1,4) %in% c("case","phen"))
  if(length(ii)>0) { pheno <- samp[[ii]] } else { pheno <- NULL }
  if(length(samp.ids)<2000 & combine) { combine <- TRUE } else { combine <- FALSE }
  snpMatLst <- vector("list",length(dat.fn))
  names(snpMatLst) <- basename(dat.fn)
  if(!combine) {
    nms <- string.in.common(basename(dat.fn),inverse=TRUE)
  }
  ofn <- character(length(dat.fn))
  for (cc in 1:length(dat.fn)) {
    snpMat <- read.impute(dat.fn[cc],rownames = samp.ids, snpcol = snpcol)
    if(HD & !combine) { 
      ofn[cc] <- cat.path(out.dir, fn=out.pref,suf=nms[cc],ext="RData")
      save(snpMat,file=ofn[cc]) 
      if(verbose) { cat("saved to",ofn[cc],"\n") }
      snpMatLst[[cc]] <- ofn[cc]
    } else {
      snpMatLst[[cc]] <- snpMat
    }
  }
  if(combine) {
    if(verbose) { cat("combining",length(snpMatLst),"objects with same samples, different regions\n") }
    snpMat <- do.call("cbind",args=snpMatLst)
    if(HD) {
      ofn <- cat.path(out.dir, pref="ALL_",string.in.common(basename(dat.fn)),ext="RData")
      save(snpMat, file=ofn)
      if(verbose) { cat("saved combined object to",ofn,"\n") }
      return(ofn)
    } else {
      return(snpMat)
    }
  } else {
    return(snpMatLst)
  }
}



                                                                                                                                                                   
  


# specific to me?
# reads snp matrices from RData binary file. if it finds multiple, will attempt to join them together
# can also handle an XSnpMatrix
get.SnpMatrix.in.file <- function(file,dir=NULL,warn=FALSE){
  if(!is.null(dir)) { file <- find.file(file,dir) }
  #print(file); print(dir)
  obj.nm <- paste(load(file))
  ## NOW MAKE SURE WE HAVE EXACTLY ONE SNPMATRIX OBJECT FROM THIS FILE ##
  if(length(obj.nm)>0) {
    if(length(obj.nm)==1) { 
      typz <- is(get(obj.nm))[1] 
    } else {
      typz <- sapply(obj.nm,function(X) { is(get(X))[1] })  
    }
    vld <- which(typz %in% c("XSnpMatrix","SnpMatrix","aSnpMatrix","aXSnpMatrix"))
    if(length(vld)<1) { stop("no SnpMatrix objects found in file:",file)}
    if(length(vld)>1) {
      if(warn) { warning("found multiple SnpMatrix objects in datafile, attempting to rbind() them") }
      concat.snp.matrix <- NULL;
      vld <- vld[order(names(typz)[vld])] # alphabetical order should ensure consistency across multiple data files
      try(concat.snp.matrix <- do.call("rbind3",args=lapply(obj.nm[vld],function(X) get(X))))
      if(is.null(concat.snp.matrix)) { 
        warning("SnpMatrix objects had different numbers of SNPs [cols], could not rbind")
        try(concat.snp.matrix <- do.call("cbind3",args=lapply(obj.nm[vld],function(X) get(X))))
      }
      if(is.null(concat.snp.matrix)) { stop("SnpMatrix objects had different numbers of Samples [rows], could not cbind either") }
      obj.nm <- "concat.snp.matrix"
    } else {
      obj.nm <- obj.nm[vld]
    }
  }
  ret.out <- get(obj.nm[1])
  if(is(ret.out)[1]=="XSnpMatrix") { if(warn) { warning("read XSnpMatrix from ",file) } }
  if(is(ret.out)[1]=="aSnpMatrix") { if(warn) { warning("read aSnpMatrix from ",file) } }
  return(ret.out)
}

#' Check whether a specific number of samples and snps will fit in a SnpMatrix object
#' 
#' SnpMatrix objects have a memory limit size of 2^31. This function can check a hypothetical
#' number of samples and SNPs (rows and columns) and returns whether this size can be represented
#' in memory. If not, suggest using BigSnpMatrix or SnpMatrixList objects to store your
#' data.
#' @param nsamp integer, proposed number of samples
#' @param nsnp integer, proposed number of SNPs
#' @return logical TRUE or FALSE, and if FALSE, also a warning message
#' @export
#' @examples
#' check.snpmat.size(2000,500000) # this dimensionality is fine
#' check.snpmat.size(15000,200000) # too big!
check.snpmat.size <- function(nsamp,nsnp) {
  maxsize <- ((2^31)-2)
  #prv(nsamp,nsnp,maxsize)
  if(!is.numeric(nsamp) | !is.numeric(nsnp)) { stop("Error: check.data.size takes numeric args")}
  size <- as.numeric(as.numeric(nsamp)*as.numeric(nsnp)); pc.max <- round(((size/2^31)*100),1)
  valid <- size<maxsize; if(is.na(valid)) { valid <- FALSE }
  if(!valid) { warning(nsamp," x ",nsnp," is ",pc.max,"% of the allowed object size for a SnpMatrix.\n",
          "Suggest splitting the files into smaller chunks (e.g, by sample subset or chromosomes)") }
  return(valid)
}

#' coerce a SnpMatrixList to a single SnpMatrix
#' 
#' also available as as(snpMatLst,"SnpMatrix")
#' @param snpMatLst a SnpMatrixList
#' @param verbose logical, whether to display information on progress
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @return returns the dataset rbind or cbind -ed together into a single SnpMatrix
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' snp.mat2 <- snpMatLst.collate(sml) # should have recreated snp.mat (except maybe not for column order)
#' snp.mat <- snp.mat[,match(colnames(snp.mat2), colnames(snp.mat))] # rearrange order of columns
#' all.equal(snp.mat,snp.mat2) # the objects are now identical
snpMatLst.collate <- function(snpMatLst, verbose=TRUE, dir=getwd()) {
  typ <- snp.mat.list.type(snpMatLst,fail=TRUE) #disk or memory
  typ2 <- get.split.type(snpMatLst,dir=dir)  # singleEntry snpSubGroups sampleSubGroups
  if(typ2=="snpmatrix") { return(snpMatLst) } # was already a SnpMatrix
  if(typ2=="singleEntry") { 
    # just 1 element, no binding required
    if(typ=="disk") {
      return(get.SnpMatrix.in.file(snpMatLst[[1]],dir=dir))
    } else {
      return(snpMatLst[[1]])
    }
  }
  spec <- get.snpMat.spec(snpMatLst,dir=dir)
  nr <- nrowL(snpMatLst,list=F)
  nc <- ncolL(snpMatLst,list=F)
  if(!check.snpmat.size(nr,nc)) { stop("Proposed SnpMatrix would be",nr,"rows by",nc,"columns, which exceeds the maximum possible size of a standard R object, suggest leaving this object as a SnpMatrixList") }
  if((typ2 %in% c("sampleSubGroups","snpSubGroups")) & typ=="disk") {
    if(verbose) { cat("Extracting",length(snpMatLst),if(typ2=="sampleSubGroups") {"sample"} else {"SNP"},"subsets from files\n") }
    snpMatLst <- lapply(snpMatLst,get.SnpMatrix.in.file,dir=dir)
  }
  if(verbose) { cat("Binding subsets to create a new SnpMatrix with dimensions:",nr,"samples x",nc,"snps\n") }
  if(typ2=="sampleSubGroups") {
    try(concat.snp.matrix <- do.call("rbind3",args=snpMatLst))
    if(is.null(concat.snp.matrix)) { 
      warning("SnpMatrix objects had different numbers of SNPs [cols], could not rbind")
      return(NULL)
    }
  }
  if(typ2=="snpSubGroups") {
    try(concat.snp.matrix <- do.call("cbind3",args=snpMatLst))
    if(is.null(concat.snp.matrix)) { 
      warning("SnpMatrix objects had different numbers of Samples [rows], could not cbind")
      return(NULL)
    }
  }
  return(concat.snp.matrix)
}


#' Get the dimensions of all of the SnpMatrix objects in a SnpMatrixList
#' 
#' Returns a table with number of samples and snps in each element
#' of a SnpMatrixList
#' @param snpMatLst a SnpMatrixList
#' @param fail logical, if TRUE, will stop() if snpMatLst is not a valid
#' SnpMatrixList, otherwise will throw a warning and try to proceed
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param print logical, if TRUE, will print the result to the console, as well
#' as returning the object
#' @return Returns a table with number of samples and snps in each element
#' of a SnpMatrixList, the first row is the sample counts, the second are SNP counts
#' @seealso nrowL, ncolL, rownamesL
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' get.snpMat.spec(sml)
#' sml <- list(snp.mat[1:20,],snp.mat[21:35],snp.mat[35:50,])
#' get.snpMat.spec(sml)
get.snpMat.spec <- function(snpMatLst,fail=T,dir=NULL,print=FALSE)
{
  # get dimensions of snpMatLst (SnpMatrix list) regardless
  # of whether it's a set of SnpMatrix objects or list of file locations
  HD <- switch(snp.mat.list.type(snpMatLst,fail),snpmatrix="snpmatrix",memory=F,disk=T,error=NULL)
  #if(is.null(dir)) { dir <- getwd() } #; warning("no directory passed to function") }
  if(HD=="snpmatrix") { return(matrix(dim(snpMatLst))) }
  if(HD) {
    n.grp <- length(snpMatLst)
    list.spec <- matrix(integer(),ncol=n.grp,nrow=2)
    for (cc in 1:n.grp)
    {
      if(!is.null(dir)) { TF <- is.file(paste(snpMatLst[[cc]]),dir) } else { TF <- file.exists(snpMatLst[[cc]]) }
      if(TF) { 
        if(!is.null(dir)) { fnm <- find.file(paste(snpMatLst[[cc]]),dir) } else { fnm <- paste(snpMatLst[[cc]]) }
        if(print) { print(fnm) }
        snpMat <- get.SnpMatrix.in.file(fnm,dir=dir)
        #snpMat <- get(paste(load())) 
        list.spec[,cc] <- dim(snpMat)
      } else {
        warning(paste("invalid snpMat file",cc)) 
      }
    }
  } else {
    list.spec <- sapply(snpMatLst,dim)
  }
  return(list.spec)
}


#' colnames function for a SnpMatrixList
#' 
#' Gets the colnames in a SnpMatrixList, which if subsets are sample-subsets
#' will return a single vector if list=FALSE, or if the subsets are based on SNPs, e.g,
#' groups/cohorts, then by default will return a list giving the rownames for
#' each subset of the list (each likely to be different)
#' @param snpMatLst a SnpMatrixList
#' @param list logical, whether to preferentially return results as a list for each
#' subset, or to sum all subset counts and give a total number
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @return returns the colnames in a SnpMatrixList, which if subsets are sample-subsets
#' will return a single vector if list=FALSE, or if the subsets are based on SNPs, e.g,
#' groups/cohorts, then by default will return a list giving the rownames for
#' each subset of the list (each likely to be different)
#' @seealso nrowL, ncolL, rownamesL
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' cn <- colnamesL(sml)
#' prv(cn)
#' cn <- colnamesL(sml,list=TRUE)
#' headl(cn,skip = 4, n = 30)
colnamesL <- function(snpMatLst,list=FALSE,dir=NULL) {
  alln <- smlapply(snpMatLst,FUN=colnames,dir=dir)
  if(list) { return(alln) }
  allN <- unlist(alln); if(length(snpMatLst)==1) { return(allN) }
  if(anyDuplicated(allN)) { 
    if(length(alln[[2]])==length(alln[[1]])) {
      if(all(alln[[1]]==alln[[2]])) { return(unlist(alln[[1]])) } # are all elements the same snps
      warning("duplicate (column) SNP names found across SnpMatrix objects in SnpMatrixList")   
    }
  }
  return(allN)
}

#' rownames function for a SnpMatrixList
#' 
#' Gets the rownames in a SnpMatrixList, which if subsets are SNP-subsets
#' will return a single vector if list=FALSE, or if the subsets are based on samples, e.g,
#' groups/cohorts, then by default will return a list giving the rownames for
#' each subset of the list (each likely to be different)
#' @param snpMatLst a SnpMatrixList
#' @param list logical, whether to preferentially return results as a list for each
#' subset, or to sum all subset counts and give a total number
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @return returns the rownames in a SnpMatrixList, which if subsets are SNP-subsets
#' will return a single vector if list=FALSE, or if the subsets are based on samples, e.g,
#' groups/cohorts, then by default will return a list giving the rownames for
#' each subset of the list (each likely to be different)
#' @seealso nrowL, ncolL, colnamesL
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- list(snp.mat[1:20,],snp.mat[21:35],snp.mat[35:50,])
#' rownamesL(sml)
#' rownamesL(sml,list=TRUE)
rownamesL <- function(snpMatLst,list=FALSE,dir=NULL) {
  alln <- smlapply(snpMatLst,FUN=rownames,dir=dir)
  if(list) { return(alln) }
  allN <- unlist(alln); if(length(snpMatLst)==1) { return(allN) }
  if(anyDuplicated(allN)) { 
    if(length(alln[[2]])==length(alln[[1]])) {
      if(all(alln[[1]]==alln[[2]])) { return(unlist(alln[[1]])) } # are all elements the same samples
      warning("duplicate (row) sample names found across SnpMatrix objects in SnpMatrixList") 
    }
  }
  return(allN)
}

#' nrow function for a SnpMatrixList
#' 
#' Gets the number of rows in a SnpMatrixList, which if subsets are SNP-subsets
#' will return a single number if list=FALSE, or if the subsets are based on samples, e.g,
#' groups/cohorts, then by default will return a list giving the number of rows for
#' each subset of the list (each likely to be different)
#' @param snpMatLst a SnpMatrixList
#' @param list logical, whether to preferentially return results as a list for each
#' subset, or to sum all subset counts and give a total number
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @return returns the number of rows in a SnpMatrixList, which if subsets are SNP-subsets
#' will be a single number if list=FALSE, or if the subsets are based on samples, e.g,
#' groups/cohorts, then by default will return a list giving the number of rows for
#' each subset of the list (each likely to be different)
#' @seealso ncolL, colnamesL, rownamesL
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' nrowL(sml)
#' nrowL(sml,list=FALSE)
nrowL <- function(snpMatLst,list=TRUE,dir=NULL) {
  out <- smlapply(snpMatLst,FUN=nrow,c.fun=unlist,dir=dir)
  if(!list)  { 
    oo <- unique(out)
    if(length(oo)==1) { out <- oo } else  { out <- sum(out) }
  }
  return(out)
}

#' ncol function for a SnpMatrixList
#' 
#' Gets the number of columns in a SnpMatrixList, which if subsets are sample-subsets
#' will return a single number if list=FALSE, or if the subsets are based on SNPs, e.g,
#' chromosomes, then by default will return a list giving the number of columns for
#' each subset of the list (each likely to be different)
#' @param snpMatLst a SnpMatrixList
#' @param list logical, whether to preferentially return results as a list for each
#' subset, or to sum all subset counts and give a total number
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @return returns the number of columns in a SnpMatrixList, which if subsets are sample-subsets
#' will be a single number if list=FALSE, or if the subsets are based on SNPs, e.g,
#' chromosomes, then by default will return a list giving the number of columns for
#' each subset of the list (each likely to be different)
#' @seealso nrowL, colnamesL, rownamesL
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 100 SNPs
#' snp.mat <- rSnpMatrix(50,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' ncolL(sml)
#' ncolL(sml,list=FALSE)
ncolL <- function(snpMatLst,list=TRUE,dir=NULL) {
  out <- smlapply(snpMatLst,FUN=ncol,c.fun=unlist,dir=dir)
  if(!list)  { 
    oo <- unique(out)
    if(length(oo)==1) { out <- oo } else  { out <- sum(out) }
  }
  return(out)
}


#' 'lapply' function for SnpMatrixList
#' 
#' Applys a function to each subset of a SnpMatrixList. Also works for normal SnpMatrix,
#' aSnpMatrix or XSnpMatrix objects.
#' @param snpMatLst a SnpMatrixList 
#' @param FUN the function to apply
#' @param c.fun function to combine the result, noting that the raw form of the result is a list, so
#' if the result should be a vector, then the appropriate combination function would be 'unlist'.
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @param ... further arguments to FUN()
#' @return by default returns a list, with the result for each outputted element
#' corresponding to the each input element in the same order. However if a 'c.fun'
#' is used, then result should correspond to the output format of 'c.fun'.
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),10) # select a random set of 10 SNPs
#' snp.mat <- rSnpMatrix(50,10); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' is(sml)
#' smlapply(sml,dim)
smlapply <- function(snpMatLst,FUN=nrow,c.fun=NULL,dir=NULL,n.cores=1,...) {
  if(is(snpMatLst)[1] %in% c("XSnpMatrix","SnpMatrix","aSnpMatrix")) { return(FUN(snpMatLst)) }
  typ <- snp.mat.list.type(snpMatLst,fail=TRUE) #disk or memory
  typ2 <- get.split.type(snpMatLst,dir=dir)  # singleEntry snpSubGroups sampleSubGroups
  if(typ2=="snpmatrix") { return(FUN(snpMatLst)) }
  if(typ2=="singleEntry") { return(FUN(snpMatLst[[1]])) }
  if(typ2=="sampleSubGroups" | typ2=="snpSubGroups") {
    #alln <- lapply(snpMatLst,rownames)
    alln <- fun.snpMatLst(snpMatLst,fun=FUN,dir=dir,n.cores=n.cores,...)
    if(is.function(c.fun)) {
      return(c.fun(alln))
    } else {
      return(alln) 
    }
  }
}

## HERE!! ###

#' Run sample quality control on a SnpMatrix or SnpMatrixList
#'
#' This function calculates summary statistics which it then uses to derive 
#' lists of samples to exclude based on the chosen thresholds which are the 
#' parameters of this function.
#' @param snpMatLst a SnpMatrixList
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param plink logical, if TRUE then will attempt to use extract the QC summary statistics
#' from plink QC files. This is mainly an option for plumbCNV.
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @param call.rate numeric, exclusion threshold for the sample callrate
#' @param het.lo numeric, lower bound threshold for exclusion based on heterozygosity,
#' default is 0.1.
#' @param het.hi numeric, upper bound threshold for exclusion based on heterozygosity,
#' default is 0.4.
#' @param subIDs.plink a parameter used by plumbCNV, leave as NULL.
#' @param proc integer, primarily for use by plumbCNV, controls the code to use for QCfail's.
#' @return a list containing SAMPLE.INFO: a data.frame of results, one row for each sample with
#' statistics like call rate, and heterozygosity, see snpStats:row.summary; elements: 
#' CR.EXCL, HZ.EXCL with vectors of sample-ids to exclude based on callrate and heterozygosity
#' thresholds.
#' @seealso doSnpQC, list.rowsummary, list.colsummary
#' @export
#' @examples
#' data(exSnpMat)
#' result <- doSampQC(exSnpMat)
#' headl(result)
#' data(exAnnotSnp)
#' result <- doSampQC(exAnnotSnp)
#' headl(result)
doSampQC <- function(snpMatLst=NULL, dir=getwd(), subIDs.plink=NULL, plink=FALSE,
                     call.rate=.95, het.lo=0.1, het.hi=0.4, n.cores=1, proc=1) 
{
  # Do sample callrate quality control on a snpMatrix object (or just use plink files
  # if the QC has already been performed in plink)
  # allow NULL for subIDs.actual if can use rownames of snpMatLst
  dir <- validate.dir.for(dir,c("ano","cr.plk"),warn=F)
  subIDs.actual <- subIDs.plink
  callrate.samp.thr <- call.rate
  if(!is.list(snpMatLst)) {  snpMatLst <- list(snpMatLst) }  # if a normal SnpMatrix is entered, this should save it                             
  sample.info <- sample.info.from.sml(snpMatLst, subIDs.actual=subIDs.actual, dir=dir, plink=plink)
  #print(head(sample.info)=dir)
  if (plink | is.null(snpMatLst)){
    if(!plink) { cat(" snpMatLst is NULL but plink=FALSE, trying to see whether plink QC data is available\n") }
    if (!"snpdataout.irem" %in% list.files(dir$cr.plk)) {
      stop(paste("Error: expecting plink samples-removed-by-QC file: snpdataout.irem, in",dir$cr.plk)) }
    if (!"snpdataout.imiss" %in% list.files(dir$cr.plk)) {
      stop(paste("Error: expecting plink samples-QC output file: snpdataout.imiss, in",dir$cr.plk)) }
    irem.fn <- cat.path(dir$cr.plk,"snpdataout.irem")
    check.fl.rows <- file.nrow(irem.fn)
    if(check.fl.rows>0) { 
      call.rate.excl.grp <- paste(read.table(irem.fn,header=F,stringsAsFactors=FALSE)[,2]) 
      sample.info[call.rate.excl.grp,"QCfail"] <- proc
    } else { 
      call.rate.excl.grp <- NULL 
    }
    cat(paste(" plink sample QC removed",length(call.rate.excl.grp),"samples\n"))
    imiss.fn <- cat.path(dir$cr.plk,"snpdataout.imiss")
    check.fl.rows <- file.nrow(imiss.fn)
    if(check.fl.rows>0) { imisser <- read.table(imiss.fn,header=T,stringsAsFactors=FALSE) } else { imisser <- data.frame(IID="",F_MISS="") }
    sample.info[["call.rate"]] <- 1-(imisser[match(rownames(sample.info),imisser$IID),"F_MISS"])
    sample.info$call.rate[is.na(sample.info$call.rate)] <- 0
    more.fail <- which(sample.info$call.rate[!is.na(sample.info$call.rate)]<callrate.samp.thr)
    if(length(more.fail)>0)
    {
      warning(paste("threshold entered (",callrate.samp.thr,") may differ from plink threshold (MIND=...) ",sep=""),
              paste("because",length(more.fail),"additional samples failed the call rate threshold. Removing additional samples.\n"),
              "Note that this will affect the SNP-QC as some snps with missing calls for these ",
              "bad samples might have failed unnecessarily.")
      sample.info$QCfail[more.fail] <- proc
      call.rate.excl.grp <- c(call.rate.excl.grp,rownames(sample.info)[more.fail])
    }
    # using plink QC, no heterozygosity exclusion?
    warning("SNP-QC with plink means samples not excluded outside heterozygosity limits")
    het.excl.grp <- paste(NULL)
  } else {
    sample.qc <- list.rowsummary(snpMatLst,dir=dir,n.cores=n.cores)
    #save(sample.info,dir,snpMatLst,sample.qc,file="fuckedupSQC.RData")
    #prv(sample.qc,sample.info)
    sample.info[["call.rate"]] <- sample.qc[rownames(sample.info),"Call.rate"]
    sample.info$call.rate[is.na(sample.info$call.rate)] <- 0
    sample.info[["heterozygosity"]] <- sample.qc[rownames(sample.info),"Heterozygosity"]
    sample.info$heterozygosity[is.na(sample.info$heterozygosity)] <- 0
    plot(density(sample.info$heterozygosity))
    call.rate.excl.grp <- paste(rownames(sample.info)[sample.info$call.rate<callrate.samp.thr])
    call.rate.excl.grp <- narm(call.rate.excl.grp)
    het.excl.grp <- paste(rownames(sample.info)[sample.info$heterozygosity<het.lo | sample.info$heterozygosity>het.hi])
    het.excl.grp <- narm(het.excl.grp)
    if(length(call.rate.excl.grp)>0) {
      sample.info[call.rate.excl.grp,"QCfail"] <- proc
    }
    if(length(het.excl.grp)>0) {
      sample.info[het.excl.grp,"QCfail"] <- proc
    }
  }
  out.list <- list(call.rate.excl.grp,het.excl.grp,sample.info)
  names(out.list) <- c("CR.EXCL","HZ.EXCL","SAMPLE.INFO")
  return(out.list)
}


#' Run SNP quality control on a SnpMatrix or SnpMatrixList
#'
#' This function calculates summary statistics which it then uses to derive 
#' lists of SNPs to exclude based on the chosen thresholds which are the 
#' parameters of this function.
#' @param snpMatLst a SnpMatrixList
#' @param snp.info If using aSnpMatrix or a SnpMatrixList comprising aSnpMatrix objects, or
#' if working with the annotation in chip.support(), then this paramter should be set to NULL.
#' Otherwise this parameter should be' a 'snp.info' annotation object, type RangedData,
#'  ChipInfo or GRanges, specifying the chromosome and position of each SNP, with the SNP labels
#'  as appearing in the SnpMatrix columns as the rownames of this object.
#' @param sample.info Not required unless you wish to consider groups of samples for outlier 
#' detection, e.g, looking for missingness differences between groups, etc. If so, there should
#' be a column 'grp' in the sample.info object, coded as integer, specifying the groups to analyse
#' separately. If using aSnpMatrix or a SnpMatrixList comprising aSnpMatrix objects then you might 
#' be able to leave this parameter as NULL. Otherwise this parameter should be' a 'sample.info' 
#' annotation object, type data.frame, specifying the phenotype, sex and other characteristics of
#' samples, with the sample labels as appearing in the SnpMatrix rownames as the rownames of this
#' object. 
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param plink logical, if TRUE then will attempt to use extract the QC summary statistics
#' from plink QC files. This is mainly an option for plumbCNV.
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @param call.rate numeric, exclusion threshold for the SNP callrate
#' @param hwe.p.thr numeric, threshold for exclusion based on Hardy Weinberg Equilibrium violation,
#' as a p-value.
#' @param grp.hwe.z.thr numeric, threshold for the worst HWE violation amongst any subgroup
#' for a SNP to be excluded, as a Z-score.
#' @param grp.cr.thr numeric, p value threshold for excluding SNPs with differences between
#' groups, if grp.miss=TRUE.
#' @param grp.miss logical, whether to exclude SNPs where there is a significant difference in 
#' missingness between cohorts/groups from different centres/datasources.
#' @param maf.thr, numeric, threshold for maf... uses 'less-than', except for zero, which
#' uses 'equal-to'. If you don't want to exclude any SNPs based on maf, use 'NA'.
#' @param autosomes.only logical, if TRUE, only analyse chromosomes 1-22
#' @param subIDs.plink a parameter used by plumbCNV, leave as NULL.
#' @param proc integer, primarily for use by plumbCNV, controls the code to use for QCfail's.
#' @param build the UCSC build to use for annotation (leave NULL if using getOption("ucsc"))
#' @return a list containing SNP.INFO: a data.frame of results, one row for each SNP with
#' statistics like call rate, and HWE, see snpStats:col.summary; elements: 
#' CR.EXCL, HWE.EXCL, and MAF.EXCL with vectors of snp-ids to exclude based on the respectively
#' callrate, hardy weinburg and MAF thresholds; elements: CR.cond, HWE.cond, and MAF.cond are
#' logical vectors for all snps to showing, where TRUE means a failure on that criterion,
#' respectively, callrate, hardy weinburg and MAF.
#' @seealso doSampQC, list.rowsummary, list.colsummary
#' @export
#' @examples
#' data(exSnpMat) # a plain SnpMatrix object
#' result <- doSnpQC(exSnpMat)
#' headl(result).
#' data(exAnnotSnp) # an aSnpMatrix object
#' result <- doSnpQC(exAnnotSnp)
#' headl(result)
#' snp.set <- sample(rownames(chip.support()),200) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(50,200); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' Dim(sml) # a aSnpMatrixList of data for ~22 chromosomes, for 50 samples, across 200 SNPS.
#' result <- doSnpQC(exAnnotSnp)
#' headl(result)
doSnpQC <- function(snpMatLst=NULL, snp.info=NULL, sample.info=NULL, dir=getwd(), plink=FALSE, n.cores=1,
                    call.rate=.95, hwe.p.thr=(10^-7), grp.hwe.z.thr=5.32, grp.cr.thr=.001, group.miss=F,
                    maf.thr=0.005, autosomes.only=T, subIDs.plink=NULL, proc=1, build=NULL) 
{
  # Do SNP callrate/HWE quality control on a snpMatrix object (or just use plink files
  # if the QC has already been performed in plink)
  # only need sample info for the group statistics
  # if entering a reasonable p-value for hwe it will auto convert to Z
  # need snp.info every time, unless snpMatLst contains aSnpMatrix objects to create one with
  dir <- validate.dir.for(dir,c("ano","cr.plk","cr"),warn=F)
  if(grp.hwe.z.thr<0.1) { grp.hwe.z.thr <- p.to.Z(grp.hwe.z.thr) }
  if(hwe.p.thr>1) { hwe.thr <- Z.to.p(hwe.p.thr) } else { hwe.thr <- hwe.p.thr }
  subIDs.actual <- subIDs.plink
  callrate.snp.thr <- call.rate
  if(is.null(snp.info)) {
    if(is.aSnpMatrix(snpMatLst,dir=dir)) {
      snp.info <- snp.info.from.annot(snpMatLst,dir=dir)
    } else { 
      if(is.null(build)) { build <- getOption("ucsc") }
      build <- ucsc.sanitizer(build)
      if(is.null(snp.info)) { snp.info <- chip.support(build=build) }
      cnX <- colnamesL(snpMatLst,list=FALSE)
      if(any(!cnX %in% rownames(snp.info))) {
        stop("all SNPs from snpMatLst must have support in the aSnpMatrix, snp.info, or the result of chip.support(build=build)")}
      snp.info <- snp.info[cnX,]
      if(!is(snp.info)[1] %in% c("GRanges","RangedData","ChipInfo")) { stop("snp.info must be GRanges, RangedData or ChipInfo") }
      snp.info <- as(snp.info,"RangedData")
    }
  }  
  if(is.null(sample.info)) {
    if(is.aSnpMatrix(snpMatLst,dir=dir)) { 
      sample.info <- sample.info.from.annot(snpMatLst,dir=dir)
    } else {
      if(snp.mat.list.type(snpMatLst) %in% c("memory","disk")) {
        sample.info <- sample.info.from.sml(snpMatLst, subIDs.actual=subIDs.actual, dir=dir, plink=plink)
      } else {
        sample.info <- NULL
      }
    }
  }
  ## get snp call rates ##
  HD <- switch(snp.mat.list.type(snpMatLst,fail=F),memory=FALSE,disk=TRUE,snpmatrix=FALSE,error=NULL)
  if(is(snpMatLst[[1]])[1] %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")) { HD <- F } else { HD <- T }
  if (plink | is.null(snpMatLst)){
    if (!"snpdataout.lmiss" %in% list.files(dir$cr.plk)) {
      stop(paste("Error: expecting plink snp-QC output file: snpdataout.lmiss, in",dir$cr.plk)) }
    if (!"snpdataout.hwe" %in% list.files(dir$cr.plk)) {
      stop(paste("Error: expecting plink snp-HWE output file: snpdataout.hwe, in",dir$cr.plk)) }
    lmisser <- read.table(paste(dir$cr.plk,"snpdataout.lmiss",sep=""),header=T,stringsAsFactors=FALSE)
    snp.info[["call.rate"]] <- 1-(lmisser[match(rownames(snp.info),lmisser$SNP),"F_MISS"])
    hweer <- read.table(paste(dir$cr.plk,"snpdataout.hwe",sep=""),header=T,stringsAsFactors=FALSE)
    pp <- snp.info[["P.hwe"]] <- hweer[match(rownames(snp.info),hweer$SNP),"P"]
    counts <- hweer[match(rownames(snp.info),hweer$SNP),"GENO"]
    three.counts <- strsplit(counts, "/",fixed=T)
    three.counts <- lapply(three.counts,function(X) { Y <- as.numeric(X[1:3]); Y[is.na(Y[1:3])] <- 0; Y })
    monoms <- sapply(three.counts,function(X){ length(which(X==0))})
    #prv(three.counts)
    mafz <- sapply(three.counts,function(X){ if(X[1]>=X[3]) { j <- X[3] } else { j <- X[1] } ; (2*j+X[2])/(sum(X)*2) })
    hetz <- sapply(three.counts,function(X){ X[2]/(sum(X)) })
    P.BB <- sapply(three.counts,function(X){ X[3]/(sum(X)) })
    snp.info[["maf"]] <- mafz
    snp.info[["het"]] <- hetz
    snp.info[["BAF"]] <- (.5*hetz + P.BB)
    cond <- snp.info$maf<.0005
    cond[is.na(cond)] <- F
    mmz <- rownames(snp.info)[cond]
    ofn <- cat.path(dir$cr,"monomorphic.txt"); cat("~wrote",length(mmz),"suspected monomorphic snps to:\n ",ofn,"\n")
    mafz <- cbind(mmz,round(snp.info[["maf"]][cond],7)); colnames(mafz) <- c("marker.id","minor.allele.frequency")
    write.table(mafz,file=ofn,quote=F,row.names=F,sep="\t")
    zz <- qnorm(1-(pp/2)) #only positive z values
    snp.info[["Z.hwe"]] <- zz
  } else {
    if(autosomes.only) { snp.info <- select.autosomes(snp.info) }
    if(group.miss & !is.null(sample.info)) {
      ## if this is done, 2 new columns will be added, if so, thresholds will be applied by:
      # apply.snp.thresholds() automatically upon detecting the colnames
      cat(" calculating SNP callrates separately for each cohort\n")
      snp.info <- get.grp.level.snp.stats(snpMatLst=snpMatLst,snp.info=snp.info,
                                          sample.info=sample.info,n.cores=n.cores,
                                          plot.fn=cat.path(dir$cr,"grpwiseCallrates.pdf"))
    }
    exp.snp <- length(chrNums(select.autosomes(snp.info))) #,22)
    # no need to convert i think as fine to recombine either way??? #
    if(FALSE & length(snpMatLst)!=exp.snp) {
      snpMatLst <- convert.smp.to.chr22(snpMatLst,snp.info=snp.info,dir=dir,n.cores=n.cores) }
    snp.qc <- list.colsummary(snpMatLst,dir=dir,n.cores=n.cores)
    ii <- match(rownames(snp.info),rownames(snp.qc))
    snp.info[["call.rate"]][!is.na(ii)] <- snp.qc[narm(ii),"Call.rate"]
    snp.info[["maf"]][!is.na(ii)] <- snp.qc[narm(ii),"MAF"]
    snp.info[["het"]][!is.na(ii)] <- snp.qc[narm(ii),"P.AB"]
    snp.info[["BAF"]][!is.na(ii)] <- ((.5*snp.qc[narm(ii),"P.AB"]) + snp.qc[narm(ii),"P.BB"])
    cond <- snp.info$maf<.0005
    cond[is.na(cond)] <- F
    mmz <- rownames(snp.info)[cond]
    ofn <- cat.path(dir$cr,"monomorphic.txt"); cat("~wrote",length(mmz),"suspected monomorphic snps to:\n ",ofn,"\n")
    mafz <- cbind(mmz,round(snp.info[["maf"]][cond],7)); colnames(mafz) <- c("marker.id","minor.allele.frequency")
    write.table(mafz,file=ofn,quote=F,row.names=F,sep="\t")
    zz <- snp.qc[narm(ii),"z.HWE"]
    pp <- (1-pnorm(as.numeric(abs(zz))))*2  #two-tailed p value
    snp.info[["P.hwe"]][!is.na(ii)] <- pp
    snp.info[["Z.hwe"]][!is.na(ii)] <- zz
  }
  out.list <- apply.snp.thresholds(snp.info,callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,
                                   grp.cr.thr=grp.cr.thr,grp.hwe.z.thr=grp.hwe.z.thr,maf=maf.thr,proc=proc)
  return(out.list)
}




#internal
get.grp.level.snp.stats <- function(snpMatLst,snp.info,sample.info,n.cores=1,plot.fn=NULL) {
  ### do some snp-qc on each group separately and look for between-grp differences
  if(length(snpMatLst)==length(unique(sample.info$grp))) {
    if(length(snpMatLst)>1) {
      # if more than one group do snp qc separately for each to test for differences between grps
      if(n.cores>1) { 
        # get SNP-qc stats summary by grp
        snp.qc.grpwise <- parallel::mclapply(snpMatLst,colSummary,filt=rownames(snp.info),mc.cores=n.cores) 
      } else { 
        snp.qc.grpwise <- lapply(snpMatLst,colSummary,filt=rownames(snp.info)) 
      }
      grplenz <- sapply(snpMatLst,nrow)
      if(is.character(plot.fn)) {
        pdf(plot.fn)
        for (dd in 1:length(snpMatLst)) {
          snp.sub.info <- snp.info 
          snp.sub.info[["call.rate"]] <- snp.qc.grpwise[[dd]][rownames(snp.info),"Call.rate"]
          draw.density.plots(fn=NULL,sample.info=sample.info[sample.info$grp==dd,],
                             snp.info=snp.sub.info, callrate.samp.thr=.9, callrate.snp.thr=.95)
        }
        dev.off()
        cat("~wrote call rate plots for each cohort to:\n ",plot.fn,"\n")
      }
      # run chi sq for number of valid calls across grps
      not.missing.cols <- sapply(snp.qc.grpwise,function(x) { x$Calls} )
      pvals <- apply(not.missing.cols,1,function(x) { 
        if(all(x==0)) { NA } else { chisq.test(x,p=(grplenz/sum(grplenz)))$p.value } } )
      snp.info[["grp.miss.p"]] <- pvals
      ## find the largest difference in HWE z-scores between cohorts for each SNP (cut ~4)
      hardy.zs <- sapply(snp.qc.grpwise, function(x) { x$z.HWE } )
      zmax.dif <- function(z) { max(abs(diff(z))) }
      hwe.z.max.difs <- apply(hardy.zs, 1, zmax.dif)
      snp.info[["grp.hwe.zmax"]] <- hwe.z.max.difs
    }       
  } else {
    warning("snp missingness comparison between groups skipped as number of groups didn't match number of snpMatrices")
  }
  return(snp.info)
}


#' SNP-wise Summary for a SnpMatrixList (SNP-qc sumamry)
#' 
#' SnpMatrixList objects can split the dataset by sample subsets, or by SNP subsets (e.g, 
#' chromosomes). This function applies the snpStats function 'row.summary' for 
#' a SnpMatrixList.
#' @param snpMatLst a SnpMatrixList
#' @param mode when mode is 'row' (default), behaviour is as described above. If mode='col'
#' then will perform a 'col.summary', a summary by SNP
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param warn logical, set TRUE to see warnings about missing genotypes
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @param verbose logical, use TRUE to show progress of each element, FALSE to hide processing text.
#' Mostly worthwhile for very large datasets where you'd like some indication of progress.
#' @return a data.frame of results, one row for each sample with statistics like call rate,
#' and heterozygosity, see snpStats:row.summary. Or, if mode='col' the frame will contain
#' 1 row per SNP, with MAF, allele frequencies, HWE, etc.
#' @seealso list.colsummary
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(20,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' result <- list.rowsummary(sml)
#' head(result)
list.rowsummary <- function(snpMatLst,mode="row",dir=getwd(),warn=T,n.cores=1,verbose=FALSE)
{
  # performs 'row.summary' or 'col.summary' snpStats functions
  # on a list of snpMatrix objects
  fail <- F
  # salvage if a regular object (not a list) is entered
  if(!fail & !is.list(snpMatLst)) { 
    if(is(snpMatLst)[1] %in% c("aSnpMatrix","XSnpMatrix","SnpMatrix")) { snpMatLst <- list(snpMatLst) } else { stop("invalid snpMatLst") }
  }
  typz <- sapply(lapply(snpMatLst,is),"[",1)
  dir <- validate.dir.for(dir,c("cr"))
  if(all(typz %in% c("aSnpMatrix","SnpMatrix","XSnpMatrix")))
  { HD <- F } else {
    if (all(typz=="character")) { HD <- T } else {
      warning("snpMatLst doesn't appear to be a list of SnpMatrix objects",
              "or a list of file locations, row.summary impossible")
      return(NULL)
    } 
  }
  dimz <- NULL
  # this line allows this function to be either row.summary or col.summary (mode)
  if(mode!="col" & mode!="row") { mode <- "row" }
  must.use.package("snpStats",T)
  wgts <- numeric(length(snpMatLst)) # vector to store number of SNPs/samples in each sublist
  my.summary <- switch(mode,row=row.summary,col=col.summary)
  dimz <- get.snpMat.spec(snpMatLst,dir=dir)
  mat.typez <- c("sampleSubGroups","snpSubGroups","singleEntry")
  mat.type <- get.split.type(snpMatLst,dir=dir)
  if(mat.type==mat.typez[3] & !HD) { return(my.summary(snpMatLst[[1]])) } # a list length 1!
  rowsum.list <- vector("list",length(snpMatLst))
  if(!warn) { 
    #avoid alarming user given known issue with genotypes that are uniformly zero (empty) in column mode
    options(warn=-1) 
  }
  if(n.cores>1 & !HD) {
    must.use.package("parallel")
    if(verbose) { cat(" processing elements in",min(n.cores,length(snpMatLst)),"parallel cores ...")  }
    rowsum.list <- mclapply(snpMatLst,my.summary)
  } else {
    if(verbose) {  cat(" processing element ") }
    for (dd in 1:length(snpMatLst))
    {
      if(verbose) { cat(dd,"..",sep="") }
      if(HD) {
        #fl.nm <- cat.path(dir$cr,snpMatLst[[dd]])
        fl.nm <- find.file(snpMatLst[[dd]],dir$cr,dir)
        #print(fl.nm)
        the.matrix <- get.SnpMatrix.in.file(file=fl.nm)
        ##############
        #print(dim(get(obj.nm))) ; badness <- which(is.na(get(obj.nm)),arr.ind=T); print(head(badness)); print(dim(badness))
        rowsum.list[[dd]] <- my.summary(the.matrix)
      } else {
        rowsum.list[[dd]] <- my.summary(snpMatLst[[dd]])
        #print(dim(snpMatLst[[dd]]))
      }
    }
  }
  if(verbose) { cat("done\n") }
  if(!warn) { options(warn=0) }
  colz <- sapply(lapply(rowsum.list,dim),"[",2)
  main <- Mode(colz)
  if(any(colz!=main)) {
    com.nms <- colnames(rowsum.list[[which(colz==main)[1]]])
    culpz <- which(colz!=main)
    for (cc in 1:length(culpz)) {
      nxt <- rowsum.list[[culpz[cc]]]
      if(all(com.nms %in% colnames(nxt))) {
        rowsum.list[[culpz[cc]]] <- nxt[,com.nms]
        warning("element number ",culpz[cc]," had more columns, so some will be truncated, likely an XSnpMatrix")
      } else {
        stop("summary column names from element number ",culpz[cc]," could not be matched to other elements.")
      }
    }
  }
  if(mode=="row" & mat.type==mat.typez[1]) {
    #sample summary on sampleSubGrps
    jj <- unlist(sapply(as.list(rowsum.list),rownames))
    out <- do.call("rbind",args=(rowsum.list))
    if(!all(rownames(out) %in% jj)) { 
      warning("sample-rownames were corrupted, fixing:\n",
              (paste(rownames(out)[1:3],"==>",jj[1:3],"\n")),"...\n",sep="")
      rownames(out) <- jj
    } 
    return(out)
  }
  if(mode!="row" & mat.type==mat.typez[2]) {
    #snp summary on snpSubGrps
    jj <- unlist(sapply(as.list(rowsum.list),rownames))
    out <- do.call("rbind",args=(rowsum.list))
    if(!all(rownames(out) %in% jj)) { 
      warning("snp-rownames were corrupted, fixing:\n",
              (paste(rownames(out)[1:3],"==>",jj[1:3],"\n")),"...\n",sep="")
      rownames(out) <- jj
    } 
    return(out)
  }
  wgts <- NA
  if(mode=="row" & mat.type==mat.typez[2]) {
    #sample summary on snpSubGrps
    wgts <- as.numeric(dimz[2,]) ; snpOnSamp <- FALSE
  }
  if(mode!="row" & mat.type==mat.typez[1]) {
    #snp summary on sampleSubGrps
    wgts <- as.numeric(dimz[1,]) ; snpOnSamp <- TRUE
  }
  if(any(is.na(wgts))) { stop("there were NAs in the list of weights (incorrect SnpMatrix dims)") }
  wgtsR <- wgts/sum(wgts)
  # now have weights (size) of each listed set
  result <- rowsum.list[[1]]; result[is.na(result)] <- 0; result <- result*0
  sampsorsnps <- rownames(result)
  # combine the results using the weights
  for (cc in 1:length(rowsum.list)) {
    for (dd in 1:ncol(result)) {
      subanalysis <- (rowsum.list[[cc]][sampsorsnps,dd]*wgtsR[cc])
      if(snpOnSamp & dd==1) { subanalysis <- subanalysis/wgtsR[cc] } # (this is 'calls' which is a count!)
      if(any(is.na(subanalysis))) {
        ## often due to non-called snps, need to avoid adding NAs, or everything will get set to NA
        if(dd==ncol(result)) {
          subanalysis[is.na(subanalysis)] <- result[is.na(subanalysis),dd]
        } else {
          subanalysis[is.na(subanalysis)] <- 0
        }
      }
      result[,dd] <- result[,dd] + subanalysis
      # prv(result[,dd])
    }
  }
  return(result)
}

#' SNP-wise Summary for a SnpMatrixList (SNP-qc sumamry)
#' 
#' SnpMatrixList objects can split the dataset by sample subsets, or by SNP subsets (e.g, 
#' chromosomes). This function applies the snpStats function 'col.summary' for 
#' a SnpMatrixList.
#' @param snpChrLst a SnpMatrixList
#' @param mode when mode is 'col' (default), behaviour is as described above. If mode='row'
#' then will perform a 'row.summary', a summary by sample.
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param warn logical, set TRUE to see warnings about missing genotypes
#' @param n.cores number of cores to use for parallel processing of large datasets
#' @return a data.frame of results, one row for each SNP with statistics like call rate,
#' MAF, allele frequencies, HWE, etc, see snpStats:col.summary. Or, if mode='row' the
#' frame will contain 1 row per sample, with call rate and heterozygosity for each.
#' @seealso list.rowsummary
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(20,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' result <- list.colsummary(sml)
#' head(result)
list.colsummary <- function(snpChrLst,mode="col",dir=getwd(),warn=F,n.cores=1)
{
  # wrapper to make 'list.rowsummary' work for 'col.summary' too
  # warn=F helps to avoid alarming user given known issue with genotypes
  # that are uniformly zero (empty) - chip qc snps, etc
  if(mode!="col" & mode!="row") { mode <- "col" }
  return(list.rowsummary(snpChrLst,mode=mode,dir=dir,warn=F,n.cores=n.cores))
}

#' Convert a SnpMatrixList which is split by sample subsets to one split by chromosomes
#' 
#' SnpMatrixList objects can split the dataset by sample subsets, or by SNP subsets (e.g, 
#' chromosomes). This function converts from the sample subset way to chromosome way.
#' @param snpMatLst a SnpMatrixList that is stored by sample subsets (e.g, by case/control/region)
#' @param snp.info if snpMatLst is an aSnpMatrix type, or if all SNPs contained are in the
#' current chip.support() object, then leave this as NULL. Otherwise this parameter should be
#' a 'snp.info' annotation object, type RangedData, ChipInfo or GRanges, specifying the
#' chromosome and position of each SNP, with the SNP labels as appearing in the SnpMatrix columns
#' as the rownames of this object. 
#' @param dir the path to the files in snpMatLst, if full paths are not used therein
#' @param build the UCSC build to use for annotation (leave NULL if using getOption("ucsc"))
#' @param verbose logical, use TRUE to show progress of each step, FALSE to hide processing text.
#' @return a SnpMatrixList reformatted into chromosome subsets rather than sample subsets.
#' @seealso convert.chr22.to.smp
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),100) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(20,100); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' rand.grp.vec <- sort(sample(1:3,20,replace=T)) # create a random grouping vector to categorize each sample
#' Dim(sml) # shows same number of rows, different numbers of columns in each element
#' sml.sub.grps <- convert.chr22.to.smp(sml, rand.grp.vec)
#' sml2 <- convert.smp.to.chr22(sml.sub.grps)
#' all.equal(sml,sml2)
convert.smp.to.chr22 <- function(snpMatLst,snp.info=NULL,dir="",n.cores=1,build=NULL,verbose=FALSE)
{
  ## convert snp.matrix list separated into different sample sets
  ## with all chromosomes, into 22 chromosome lists with all samples
  ## HD=T saves chromosome list elements to disk (loc) instead of memory
  ## in a situation where data is very large or memory is limited
  ## In HD=T the output will be a list of disk locations rather than an array
  dir <- validate.dir.for(dir,"cr")
  #print(snp.mat.list.type(snpMatLst,fail=T))
  typ <- is(snpMatLst[[1]])[1]
  if(is.null(snp.info) & (typ %in% c("aSnpMatrix","aXSnpMatrix")))  { snp.info <- snp.info.from.annot(snpChrLst) }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(snp.info)) { snp.info <- chip.support(build=build) }
  cnX <- colnamesL(snpMatLst,list=FALSE)
  if(any(!cnX %in% rownames(snp.info))) {
    stop("all SNPs from snpMatLst must have support in the aSnpMatrix, snp.info, or the result of chip.support(build=build)")}
  snp.info <- snp.info[cnX,]
  if(!is(snp.info)[1] %in% c("GRanges","RangedData","ChipInfo")) { stop("snp.info must be GRanges, RangedData or ChipInfo") }
  snp.info <- as(snp.info,"RangedData")
  HD <- switch(snp.mat.list.type(snpMatLst,fail=T),memory=F,disk=T,error=NULL)
  chr.set <- chrNames(snp.info); n.chr <- length(chr.set)
  #print(chr.set)
  if(verbose) { cat("\nConverting from 'n' group list of all markers to",n.chr,"chromosome list of all samples\n") }
  snpChrLst <- vector("list", n.chr)
  #num.subs <- length(subIDs.actual)
  if (is(snp.info)[1]=="RangedData")
  {
    must.use.package("genoset",bioC=T)
  } else {
    stop("snp.info parameter must be a snp.info object")
  }
  must.use.package("snpStats",bioC=T)
  list.spec <- get.snpMat.spec(snpMatLst,dir=dir)
  #prv(list.spec)
  rownames(list.spec) <- c("Samples","SNPs"); 
  if(!is.null(colnames(list.spec))) { colnames(list.spec) <- basename(colnames(list.spec)) }
  if(verbose) { cat("\n"); print(list.spec); cat("\n") }
  #e.g, sanger    t1d    uva
  #[1,]   4537   6808   5461   --> samples
  #[2,] 196524 196524 196524   --> snps
  if(any(diff(range(list.spec[2,]))!=0)) { cat("Warning! SNP files have different numbers of markers\n") }
  st <- 1+c(0,cumsum(list.spec[1,])[1:(ncol(list.spec)-1)])
  en <- cumsum(list.spec[1,])
  num.subs <- sum(list.spec[1,],na.rm=T)
  # define a local function to do each chromosome #
  if(HD) { 
    rnm <- paste(1:num.subs)
  } else {
    rnm <- unlist(c(sapply(snpMatLst,rownames)))
  }
  #return(snp.info)
  one.chr <- function(cc) {
    # do all the procesing for one chromosome
    # snp.info is a RangedData object (IRanges)
    nch <- paste(chr.set)[cc]
    #if(is.numeric(nch)) { if(nch>22) { nch <- chrNames(snp.info)[cc] } }
    next.chr <- rownames(snp.info[nch]) 
    #print(nch) ; print(snp.info[next.chr,])
    next.mat <- matrix(raw(),nrow=num.subs,ncol=length(next.chr))
    colnames(next.mat) <- next.chr; rownames(next.mat) <- rnm 
    next.snpmat <- new("SnpMatrix", next.mat)
    rm(next.mat) #; gc() this was really slow!
    for (dd in 1:length(snpMatLst))
    {
      if(verbose) { cat(".") }
      if(HD){
        snpMat <- get(paste(load(paste(snpMatLst[[dd]]))))
        col.select <- next.chr %in% colnames(snpMat)
        next.snpmat[st[dd]:en[dd],col.select] <- snpMat[,next.chr[col.select]]
        rownames(next.snpmat)[st[dd]:en[dd]] <- rownames(snpMat)
        rm(snpMat)
      } else {
        col.select <- next.chr %in% colnames(snpMatLst[[dd]])
        next.snpmat[st[dd]:en[dd],col.select] <- snpMatLst[[dd]][,next.chr[col.select]]
        rm(col.select)
      }
    }
    if(HD) {
      fnm <- cat.path(dir$cr,"chrmat_",suf=chr.set[cc],ext="RData")
      save(next.snpmat,file=fnm)
      return(fnm)
    } else {
      return(next.snpmat)
    }
  }
  ## run for each chromosome in serial or parallel mode
  if(length(chrNames(snp.info))!=length(snpChrLst)) {
    stop("number of chromosomes in snp.info [",length(chrNames(snp.info)),"] did not match snpChrLst [",length(snpChrLst),"]") }
  if(n.cores>1) {
    if(verbose) { cat(" processing chromosomes",paste(range(chr.set),collapse="-")) }
    must.use.package("parallel"); snpChrLst <- mclapply(1:n.chr,one.chr,mc.cores=min(n.cores,4))
    if(verbose) { cat(" .. done\n") }
  } else {
    if(verbose) { cat(" processing chromosome ") }; snpChrLst <- vector("list",n.chr)
    for (j in 1:n.chr) {
      if(verbose) { cat(chr.set[j]) }; snpChrLst[[j]] <- one.chr(j)
    }
    if(verbose) { cat("..done\n") }
  }
  return(snpChrLst)
}



#' Convert a SnpMatrixList which is split by chromosome to one split by sample subsets
#' 
#' SnpMatrixList objects can split the dataset by sample subsets, or by SNP subsets (e.g, 
#' chromosomes). This function converts from the chromosome way to sample subset way.
#' @param snpChrLst a SnpMatrixList that is stored by snp subsets (usually by chromosome)
#' @param snp.info if snpChrList is an aSnpMatrix type, or if all SNPs contained are in the
#' current chip.support() object, then leave this as NULL. Otherwise this parameter should be
#' a 'snp.info' annotation object, type RangedData, ChipInfo or GRanges, specifying the
#' chromosome and position of each SNP, with the SNP labels as appearing in the SnpMatrix columns
#' as the rownames of this object. 
#' @param subIDs optional list of subject IDs to include in the new object, if NULL, the ids will
#' be take from the rownames of 'snpChrLst'. Only necessary if you want to select a subset,
#' to reorder the samples, or both.
#' @param group.nums vector of integers with the same length as 'subIDs' or else the same as
#' the number of samples in 'snpChrLst'. If left as NULL, a random vector, with the same
#' number of groupings as the number of chromosomes in snpChrLst will be created to roughly
#' maintain the size of the subsets. Generally most imagineable uses for this function should
#' provide some sensible grouping vector here, rather than the random default.
#' @param dir the path to the files in snpChrLst, if full paths are not used therein
#' @param build the UCSC build to use for annotation (leave NULL if using getOption("ucsc"))
#' @param verbose logical, use TRUE to show progress of each step, FALSE to hide processing text.
#' @return a SnpMatrixList reformatted into sample subsets rather than SNP (chromosome) subsets.
#' @seealso convert.smp.to.chr22
#' @export
#' @examples
#' snp.set <- sample(rownames(chip.support()),200) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(50,200); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,]) # create a SnpMatrixList
#' rand.grp.vec <- sort(sample(1:3,50,replace=T)) # create a random grouping vector to categorize each sample
#' Dim(sml) # shows same number of rows, different numbers of columns in each element
#' sml.sub.grps <- convert.chr22.to.smp(sml, rand.grp.vec)
#' Dim(sml.sub.grps) # shows same number of columns, different numbers of rows in each element
#' convert.chr22.to.smp(sml[1:5], group.nums=rand.grp.vec, verbose=TRUE) # show progress (limit to first 5 chr's)
convert.chr22.to.smp <- function(snpChrLst,group.nums=NULL,snp.info=NULL,subIDs=NULL,dir=NULL,build=NULL,verbose=FALSE)
{
  ## convert snp.matrix list separated into 22+ chromosome lists with all samples
  ## into different sample sets with all chromosomes. 
  ## HD=T saves chromosome list elements to disk (loc) instead of memory
  ## in a situation where data is very large or memory is limited
  ## In HD=T the output will be a list of disk locations rather than an array
  typ <- is(snpChrLst[[1]])[1]
  if(is.null(snp.info) & (typ %in% c("aSnpMatrix","aXSnpMatrix")))  { snp.info <- snp.info.from.annot(snpChrLst) }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(snp.info)) { snp.info <- chip.support(build=build) }
  cnX <- colnamesL(snpChrLst,list=FALSE)
  if(any(!cnX %in% rownames(snp.info))) {
    stop("all SNPs from snpChrLst must have support in the aSnpMatrix, snp.info, or the result of chip.support(build=build)")}
  snp.info <- snp.info[cnX,]
  if(!is(snp.info)[1] %in% c("GRanges","RangedData","ChipInfo")) { stop("snp.info must be GRanges, RangedData or ChipInfo") }
  snp.info <- as(snp.info,"RangedData")
  HD <- switch(snp.mat.list.type(snpChrLst,fail=T),memory=F,disk=T,error=NULL)
  if (is(snp.info)[1]=="RangedData")
  { must.use.package("genoset",bioC=T)
  } else { stop("snp.info object was not of type 'RangedData'")}
  chr.set <- chrNums(snp.info); n.chr <- length(chr.set)
  if(is.null(subIDs)) { subIDs <- rownamesL(snpChrLst,list=FALSE) }
  if(is.null(group.nums)) { group.nums <- split.vec(length(subIDs),n.chr) } # use the same number of sample grps as chr grps
  num.grps <- length(unique(group.nums))
  snpMatLst <- vector("list", num.grps)
  if(length(group.nums)!=length(subIDs)) { stop("'group.nums' must have the same length as subIDs") }
  num.subs <- length(subIDs)
  grp.sizes <- table(group.nums)
  if(verbose) { 
    cat("\nConverting from",n.chr,"chr list of all samples to",length(unique(group.nums)),"sub-group list of all markers\n")
  }
  must.use.package("snpStats",bioC=T)
  list.spec <- get.snpMat.spec(snpChrLst,dir=dir) # even with 1 chr, should still have 2 rows
  #e.g, chr1  chr2  chr3  chr4  chr5  chr6  chr7snp  , ...   ... 
  #[1,] 16806 16806 16806 16806 16806 16806 16806 <-- samples
  #[2,] 23314 21355 12307  6444 13768 22316  7751  <-- markers
  if((any(diff(range(list.spec[1,]))!=0))) { cat("Warning! SNP files have different numbers of samples\n") }
  st <- 1+c(0,cumsum(list.spec[2,])[1:(ncol(list.spec)-1)])
  en <- cumsum(list.spec[2,])
  num.snps <- rowSums(list.spec)[2]
  for (dd in 1:num.grps) {
    if(verbose) { cat(paste(" writing group",dd,"\n")) }
    next.mat <- matrix(raw(),nrow=grp.sizes[dd],ncol=num.snps)
    colnames(next.mat) <- rownames(snp.info)
    next.grp <- subIDs[group.nums==dd]
    rownames(next.mat) <- next.grp
    next.snpmat <- new("SnpMatrix", next.mat)
    for (cc in 1:n.chr)
    {
      if(HD){
        snpMat <- get(paste(load(paste(snpChrLst[[cc]]))))
        row.select <- next.grp %in% rownames(snpMat)
        next.snpmat[row.select,st[cc]:en[cc]] <- snpMat[next.grp[row.select],]
        snpMat <- NULL
      } else {
        if(verbose) {  cat(" processing chromosome ",chr.set[cc],"...",sep="") }
        snpChrLstCC <- snpChrLst[[cc]]
        row.select <- next.grp %in% rownames(snpChrLst[[cc]])
        next.snpmat[row.select,st[cc]:en[cc]] <- snpChrLst[[cc]][next.grp[row.select],]
        if(verbose) { cat("..done\n") }
      }
    }
    if(HD) {
      fnm <- cat.path(dir$cr,"snpmat_",suf=dd,ext="RData")
      save(next.snpmat,file=fnm)
      snpMatLst[[dd]] <- fnm
    } else {
      snpMatLst[[dd]] <- next.snpmat
    }
    if(verbose) { cat(" complete\n") }
  }
  return(snpMatLst)
}



#' Create a full sample.info object from an aSnpMatrix or SnpMatrixList
#' 
#' The package annotSnpStats simplifies and makes robust annotation 
#' that assists with several aspects of conducting a GWAS analysis in R. 
#' This function extracts a standardized 'sample.info' formatted data.frame used
#' by some functions in humarray and plumbCNV from an aSnpMatrix object, or a 
#' SnpMatrixList containing aSnpMatrix objects.
#' @param X an aSnpMatrix, or an SnpMatrixList containing aSnpMatrix objects
#' @param dir path to the SnpMatrixList files [soon redundant]
#' @return returns a data.frame based on the aSnpMatrix@samples slot, or
#' combines these slots for all sample subsets to produce a data.frame for
#' all sample sets (if there are multiple).
#' @export
#' @examples
#' data(exAnnotSnp)
#' si <- sample.info.from.annot(exAnnotSnp)
#' head(si) # preview the resulting annotation object
sample.info.from.annot <- function(aSnpMatLst,dir=getwd()) {
  if(!is.aSnpMatrix(aSnpMatLst,dir=dir)) { stop("not an aSnpMatrix object") }
  if(is.list(aSnpMatLst)) {
    typ <- get.split.type(aSnpMatLst,dir=dir,allow.dif.dims=TRUE)
    if(typ=="sampleSubGroups") {
      outlist <- smlapply(aSnpMatLst,FUN=samp.from.annot,dir=dir)
      for (cc in 1:length(outlist)) {
        outlist[[cc]][["grp"]] <- rep(cc,nrow(outlist[[cc]]))
      }
      out <- do.call("rbind",args=outlist)
      return(out)
    } else {
      return(samp.from.annot(snpMatLst.collate(aSnpMatLst[1],dir=dir)))
    }
  } else {
    return(samp.from.annot(aSnpMatLst))
  }
}

#' Create a full snp.info object from an aSnpMatrix or SnpMatrixList
#' 
#' The package annotSnpStats simplifies and makes robust annotation 
#' that assists with several aspects of conducting a GWAS analysis in R. 
#' This function extracts a standardized 'snp.info' formatted RangedData object used
#' by some functions in humarray and plumbCNV from an aSnpMatrix object, or from a 
#' SnpMatrixList containing aSnpMatrix objects.
#' @param X an aSnpMatrix, or an SnpMatrixList containing aSnpMatrix objects
#' @param dir path to the SnpMatrixList files [soon redundant]
#' @return returns a RangedData object based on the aSnpMatrix@snps slot, or
#' combines these slots for all chromosomes to produce combined annotation for
#' all sample sets (if there are multiple).
#' @export
#' @seealso aSnpMatrix.to.ChipInfo
#' @examples
#' data(exAnnotSnp)
#' si <- snp.info.from.annot(exAnnotSnp)
#' head(si) # preview the resulting annotation object
snp.info.from.annot <- function(aSnpMatLst,dir=getwd(),alleles=FALSE) {
  if(!is.aSnpMatrix(aSnpMatLst,dir=dir)) { stop("not an aSnpMatrix object") }
  if(is.list(aSnpMatLst)) {
    typ <- get.split.type(aSnpMatLst,dir=dir)
    if(typ=="snpSubGroups") {
      outlist <- smlapply(aSnpMatLst,FUN=snp.from.annot,dir=dir,alleles=alleles)
      out <- do.call("rbind",args=outlist)
      return(out)
    } else {
      return(snp.from.annot(snpMatLst.collate(aSnpMatLst[1],alleles=alleles)))
    }
  } else {
    return(snp.from.annot(aSnpMatLst,alleles=alleles))
  }
}

# create a (partial) snp.info object from a single aSnpMatrix object
snp.from.annot <- function(aSnpMat,alleles=FALSE) {
  if(!is(aSnpMat)[1]=="aSnpMatrix") { stop("aSnpMat must be of type 'aSnpMatrix'") }
  df <- aSnpMat@snps
  cn <- colnames(df)
  chr.col <- which(tolower(substr(cn,1,2)) %in% "ch")
  pos.col <- which(tolower(substr(cn,1,2)) %in% "po")
  DF <- df.to.ranged(df,start=cn[pos.col],end=cn[pos.col],chr=cn[chr.col],GRanges=FALSE)
  if(alleles) {
    al <- paste(aSnpMat@alleles)
    if(all(nchar(al)>0) & length(al)>1) {
      DF[["A1"]] <- df[[al[1]]]
      DF[["A2"]] <- df[[al[2]]]
    } else {
      warning("alleles option was set TRUE but did not find anything in aSnpMat@alleles")
    }
  }
  DF[["QCfail"]] <- rep(0,nrow(DF))
  return(DF)
}






#' Create an ChipInfo support object from an aSnpMatrix object
#' 
#' This function links the annotSnpStats package to the humarray
#' package.
#' The package annotSnpStats simplifies and makes robust annotation 
#' that assists with several aspects of conducting a GWAS analysis in R. 
#' This function allows addition of full humarray package annotation
#' based on the working annotSnpStats aSnpMatrix object (or list thereof).
#' @param X an aSnpMatrix, or list of aSnpMatrix objects
#' @param incorporate logical, if TRUE, then attempt to save the resulting
#' object and set the getOption('chip.info') value to this new file. If
#' FALSE, simply return the ChipInfo object produced.
#' @param location the path to save the derived object. Can be a file or
#' a directory (in which case a filename will be generated using today's date).
#' If NULL, a file will be saved in the tempdir() (temporary R session folder).
#' @return returns a ChipInfo object based on the aSnpMatrix@snps slot, or
#' if incorporate = TRUE, then saves this to a file and points humarray to 
#' this location for use as the default annotation for the current session.
#' @export
#' @examples
#' #ownexample
#' data(exAnnotSnp)
#' ci <- aSnpMatrix.to.ChipInfo(exAnnotSnp)
#' ci # preview the resulting annotation object
#' # aSnpMatrix.to.ChipInfo(exAnnotSnp,incorporate=TRUE) # would set the working chip.support to this result #
aSnpMatrix.to.ChipInfo <- function(X,incorporate=FALSE,location=NULL) {
  #if(!is(X)[1] %in% c("aSnpMatrix","aXSnpMatrix")) { stop("X must be an aSnpMatrix object (package: annotSnpStats)") }
  snp.info <- snp.info.from.annot(X,alleles=TRUE)
  ci <- as(snp.info,"ChipInfo")
  if(incorporate) {
    if(is.character(location)) {
      if(get.ext(location)=="" & file.exists(location)) {
        # assume a directory was passed, so need to add filename
        location <- cat.path(location,"aSnpMatChipInfo",suf=simple.date(),ext="RData")
      }
    } else {
      location <- cat.path(tempdir(),"aSnpMatChipInfo",suf=simple.date(),ext="RData")
    }
success <- T
success <- tryCatch(save(ci ,file=location),error=function(e) { F } )
if(!is.logical(success)) { success <- T }
if(success) {
  options(chip.info=location)
} else {
  stop("'location' for the chip.support() object was invalid, save() failed")
}
  } else {
    return(ci)
  }
}



#' Convert a SnpMatrix to a BigSnpMatrix
#'
#' BigSnpMatrix objects allow for larger datasets that can be stored
#' in a regular SnpMatrix. This function allows conversion of a SnpMatrix, or aSnpMatrix
#' into a BigSnpMatrix, which is highly suitable to running ancestry PCA.
#' @param X a SnpMatrix or aSnpMatrix object
#' @param replace.missing logical, whether to replace any missing values
#' before creating the BigSnpMatrix. This is highly recommended, especially for PCA
#' as the missing values are stored as zeros, which would result in a skewed
#' representation of the data for many functions, such as PCA.
#' @param ref.data optional reference SnpMatrix containing data with the desired
#' allele frequencies to impute to, or otherwise will use the frequencies of 
#' the current dataset in snpMatLst.
#' @param filename base file name prefix for the objects on disk that will store
#' the BigSnpMatrix object, include a path if not wanting to write to the working
#' directory.
#' @param tracker logical, whether to show a progress bar during processing
#' @param n.cores integer, this function can be run with a multicore setup, use
#' this parameter to specify the maximum number of cores to utilise for the time
#' consuming step of missing value replacement.
#' @param limit.ram logical, whether to regularly flush RAM during processing
#' of large objects (can create some long pauses during processing, but results
#' in more conservative, stable, memory management). Usually works ok when set FALSE.
#' @param mn.zero logical, set TRUE to mean correct the values in the dataset; if
#' ref.data is in use, will correct using the SNP means in that dataset, if not,
#' then will correct to zero (using the means in the snpMatLst dataset). For PCA
#' you would likely set this option to TRUE.
#' @param sd.hwe logical, set TRUE to SD correct the values in the dataset; if
#' ref.data is in use, will correct using the SNP means in that dataset, if not,
#' then will correct using the means in the snpMatLst dataset. Either way this
#' correction is done by using the expected standard deviation of the reference
#' allele frequency at each SNP if it were in perfect Hardy Weinburg Equilibrium.
#' For PCA you would often set this option to TRUE.
#' @return returns a BigSnpMatrix object which is essentially a big.matrix object
#' from the bigmemory() package, where genotypes are coded 1,2,3 (aa, ab, bb) where
#' b is the reference allele and zero is missing. However if mn.zero or sd.hwe are
#' TRUE, then the values will be post correction, so will not be such neat integer
#' values for each genotype. PCA can be performed on this object using the bigpca package.
#' @export
#' @seealso snp.from.annot, chip.support, bigpca
#' @examples 
#' data(exSnpMat)
#' big.mat <- bigSnpMatrix(exSnpMat,mn.zero=FALSE,sd.hwe=FALSE) # genotype values are 1,2,3
#' big.mat <- bigSnpMatrix(exSnpMat) # genotype values are now continuous, and ready for PCA
bigSnpMatrix <- function(X,filename="tempMatrix",tracker=TRUE, n.cores=1,
                         limit.ram=FALSE, mn.zero=TRUE, sd.hwe=TRUE, replace.missing=TRUE, ref.data=NULL) {
  aaa <- proc.time()[3]
  typ <- is(X)[1]
  typ2 <- is(ref.data)[1]
  max.gb <- 2.5
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) {
      stop("X must be a SnpMatrix or aSnpMatrix object")
  }
  if(!is.null(ref.data)) {
    if(!typ2 %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) {
      warning("ref.data should be a SnpMatrix or aSnpMatrix object, will ignore")
      ref.data <- NULL
    } else {
      if(all(clean.snp.ids(colnames(X)) %in% clean.snp.ids(colnames(ref.data)))) {
        ref.data <- ref.data[,match(clean.snp.ids(colnames(X)),clean.snp.ids(colnames(ref.data)))]
      } else {
        warning("ref.data should contain all column names in X, will ignore")
      }
    }
  }
  bbb <- proc.time()[3]; cat("ref data took ",round(bbb-aaa)," seconds\n")
  if(replace.missing) {
    X <- randomize.missing(X,verbose=TRUE,n.cores=n.cores) # was rm2()
    nmis <- length(which(X==as.raw("00")))
    if(nmis>0) {
      warning(out.of(nmis,do.call("*",args=as.list(Dim(X))))," data points are still missing from X")
    }
  }
  ccc <- proc.time()[3]; cat("missing data took ",round(ccc-bbb)," seconds\n")
  raw <- X@.Data
  if(mn.zero | sd.hwe) {
    if(!is.null(ref.data)) {
      oo <- correct.by(ref.data); mN <- oo$mean; sD <- oo$sd
    } else {
      oo <- correct.by(X); mN <- oo$mean; sD <- oo$sd
    }
  } 
  ddd <- proc.time()[3]; cat("mean/sd calc took ",round(ddd-ccc)," seconds\n")
  #return(oo)
  #raw[raw==0] <- NA
  nC <- ncol(X); nR <- nrow(X)
  if(!mn.zero) { mN <- rep(0,nC) }
  if(!sd.hwe) { sD <- rep(1,nC) }
  sd.er <- (length(which(is.na(sD)))) 
  mn.er <- (length(which(is.na(mN))))
  if(sd.er>0)  { warning("St. Dev. could not be calculated for ",sd.er,"SNPs due to missingness; suggest using 'rmv.common.missing()' with ref.data prior to analysis") }
  if(mn.er>0)  { warning("Mean offset could not be calculated for ",mn.er,"SNPs due to missingness; suggest using 'rmv.common.missing()' with ref.data prior to analysis") }
  #prv(mN,sD)
  dir <- dirname(filename); filename <- basename(filename)
  if(dir=="") { dir <- getwd() }
  if(tracker) { cat(" creating",nR,"x",nC,"target matrix,",cat.path(dir,filename),"...\n") }
  des <- paste(filename,"dsc",sep=".")
  bck <- paste(filename,"bck",sep=".")
  bigSnp <- big.matrix(nrow=nR,ncol=nC, backingfile=bck, dimnames=list(rownames(X),colnames(X)),
                         backingpath=dir, descriptorfile=des)
  eee <- proc.time()[3]; cat("creating big.matrix took ",round(eee-ddd)," seconds\n")
  split.to <- round(10*estimate.memory(bigSnp)) # split into .1GB chunks, save RAM without creating groups too small to process
  #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
  stepz <- round(seq(from=1,to=nC+1,length.out=round((split.to+1))))
  if((tail(stepz,1)) != nC+1) { stepz <- c(stepz,nC+1) }
  split.to <- length(stepz)-1
  prv(split.to)
  for (cc in 1:split.to)
  {
    # within submatrix cols
    c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
    # do the copying
    lilColRange <- c(c1:c2)
    if(tracker) {      loop.tracker(cc,split.to) }
    if(is.finite(sum(as.numeric(lilColRange)))) {
      #cat(range(lilColRange)); cat(dim(bigTrans)); cat(dim(bigMat))
      nxt.chunk <- raw[1:nR,lilColRange]
      nxt.chunk <- t(apply(nxt.chunk,2,as.numeric))
      new.val <- (nxt.chunk-mN[lilColRange])/sD[lilColRange]
      bigSnp[1:nR,lilColRange] <- t(new.val)
      #prv(new.val,nxt.chunk,lilColRange,c1,c2,cc)
    } else {
      cat(" Warning: empty interval ignored\n")
    }
    if(limit.ram & (cc %% (max.gb*10) == 0)) {
      # reset memory after every 'max.gb' 1GB chunks to prevent skyrocketing RAM use #
      fl.suc <- bigmemory::flush(bigSnp) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()  
      if(T) {
        RR <- describe(bigSnp); rm(bigSnp); bigSnp <- attach.big.matrix(RR,path=dir)
      }
    }
  }
  fff <- proc.time()[3]; cat("loop took ",round(fff-eee)," seconds\n")
  if(tracker) {
    prv(bigSnp)
  }
  return(bigSnp)
}




#' Calculate a full linkage disequilibrium matrix of large size
#'
#' The ld() function in SnpStats can only give a full matrix (not sparse)
#' as a result when the square of the number of SNPs is less than the
#' maximum size for a data.frame. For LD-pruning of a large dataset it
#' may be desirable to calculate LD for an entire chromosome. This function
#' runs the calculations in subsets and feeds the results to a big.matrix
#' object which has a much larger capacity, and is able to produce a full
#' LD matrix for very large SnpMatrix objects, which can be processed with
#' multiple cores (increases speed).
#' @param X a SnpMatrix/aSnpMatrix
#' @param stats passed to snpStats:ld(stats=), default is "R.squared", can
#' also calculate "D.prime", "R" and others. For this function only the first
#' statistic will be used.
#' @param n.cores maximum number of cores to use for parallel processing of
#' ld() calculations.
#' @param filename base file name prefix for the big.matrix objects on disk that will 
#' store the LD matrix result, include a path if not wanting to write to the working
#' directory.
#' @param tracker logical, whether to show a progress bar during processing
#' @param limit.ram logical, whether to regularly flush RAM during processing
#' of large objects (can create some long pauses during processing, but results
#' in more conservative, stable, memory management). Usually works ok when set FALSE.
#' @return returns a big.matrix object  from the bigmemory() package, 
#' containing R.squared (or other stat parameter) values. Can be used for pruning selection
#' and other purposes.
#' @export
#' @seealso ld.prune
#' @examples 
#' data(exSnpMat)
#' result <- big.ld(exSnpMat)
#' prv(result) # R.squared matrix
#' result <- big.ld(exSnpMat,stats="D.prime")
#' prv(result) # D.prime matrix
big.ld <- function(X,stats="R.squared",n.cores=1,filename="ldmat",tracker=F,limit.ram=FALSE) {
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { stop("X must be a SnpMatrix/aSnpMatrix") }
  n.snp <- ncol(X)
  filename <- paste0(filename,sample(10^6,1))
  dir <- dirname(filename)
  max.gb <- 4
  if(all(dir==".")) { dir <- getwd() }
  des <- paste(filename,"dsc",sep=".")
  bck <- paste(filename,"bck",sep=".")
  if(n.cores>1) { multi <- T } else { multi <- F }
  #prv(n.snp,bck,X,dir,des)
  bigLD <- big.matrix(nrow= n.snp,ncol= n.snp, backingfile=bck, dimnames=list(colnames(X),colnames(X)), backingpath=dir, descriptorfile=des)
  split.to <- round(100*estimate.memory(n.snp^2)) # split into .1GB chunks, save RAM without creating groups too small to process
  #prv(split.to)
 # if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
  stepz <- round(seq(from=1,to=n.snp+1,length.out=round((split.to+1))))
  if((tail(stepz,1)) != n.snp+1) { stepz <- c(stepz,n.snp+1) }
  split.to <- length(stepz)-1
  #prv(bigLD,stepz)
  if(multi) {
    job.count <- 0 
    cc.collect <- numeric(n.cores)
    runs <- lilColRange <- vector("list",n.cores)
    for (cc in 1:split.to)
    {
      job.count <- job.count + 1
      cc.collect[job.count] <- cc
      # within submatrix cols
      c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
      if(c2>ncol(X)) { prv(c2) } # remove this line, for debugging only
      if(c1<1) { prv(c1) } # remove this line, for debugging only
      if(is.na(c1) | is.na(c2)) { prv(c1,c2) } # remove this line, for debugging only
      # do the copying
      lilColRange[[cc]] <- c(c1:c2)
      if(tracker) { loop.tracker(cc,split.to) }
      runs[[job.count]] <- parallel::mcparallel(
        {
          if(is.finite(sum(lilColRange[[cc]]))) {
            ld(X,X[, lilColRange[[cc]]],stats=stats)
          } else {  
            NULL
          }
      })
      ## after a run of n.cores jobs, consolidate
      if(job.count>=n.cores | cc>=split.to) {
        ## collect previous #'n.cores' runs 
        list.set <- parallel::mccollect(runs[1:job.count])
        for (dd in (which(cc.collect>0))) {
          if(!is.null(list.set[[dd]])) {
            bigLD[,lilColRange[[cc.collect[dd]]]] <- list.set[[dd]]
          }
        }
        job.count <- 0; runs <- vector("list",n.cores); cc.collect <- numeric(n.cores)
        if(limit.ram & (cc %% (max.gb*10) == 0)) {
          # reset memory after every 'max.gb' 1GB chunks to prevent skyrocketing RAM use #
          fl.suc <- bigmemory::flush(bigLD) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()
          if(T) {
            RR <- describe(bigLD); rm(bigLD); bigLD <- attach.big.matrix(RR,path=dir)
          }
        }
      }
      ####
      loop.tracker(cc,split.to)
    }
  } else {
    for (cc in 1:split.to)
    {
      # within submatrix cols
      c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
      # do the copying
      lilColRange <- c(c1:c2)
      if(tracker) {      loop.tracker(cc,split.to) }
      if(is.finite(sum(as.numeric(lilColRange)))) {
        nxt.chunk <- ld(X,X[, lilColRange],stats=stats)
        bigLD[,lilColRange] <- nxt.chunk
      } else {
        cat(" Warning: empty interval ignored\n")
      }
      if(limit.ram & (cc %% (max.gb*10) == 0)) {
        # reset memory after every 'max.gb' 1GB chunks to prevent skyrocketing RAM use #
        fl.suc <- bigmemory::flush(bigLD) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()
        if(T) {
          RR <- describe(bigLD); rm(bigLD); bigLD <- attach.big.matrix(RR,path=dir)
        }
      }
    }
  }
  return(bigLD)
}


# MY personal EXAMPLE
# sml <- list.files("/chiswick/data/ncooper/iChipData/",pattern="temp.ichip")
# sml <- as.list(sml)
# sml <- sml[c(11,15:22,1:10,12:14,24,23)]
# dd <- "/chiswick/data/ncooper/iChipData/"
# dd <- ("/chiswick/data/ncooper/imputation/T1D/PQDATA/")
# sml <- list.files(dd)
# sml <- as.list(sml)
# sml <- sml[c(18:19,26:40,1:17,20:24)]

#' Convert a SnpMatrixList to a single BigSnpMatrix
#'
#' SnpMatrixList objects allow for larger datasets that can be stored
#' in a regular SnpMatrix. Sometimes however, certain functions or operations
#' require that the whole dataset be accessible at once, such as an ancestry
#' PCA. This function allows conversion of a SnpMatrixList into a BigSnpMatrix, 
#' which allows a SnpMatrix of almost any size to be created (at a substantial cost of speed).
#' @param X a SnpMatrixList
#' @param verbose logical, whether to supply information about progress
#' @param dir location of the SnpMatrixList elements, not required if the
#' elements are full path locations. Default is the working directory.
#' @param replace.missing logical, whether to replace any missing values
#' before creating the BigSnpMatrix. This is highly recommended, especially for PCA
#' as the missing values are stored as zeros, which would result in a skewed
#' representation of the data for many functions, such as PCA.
#' @param ref.data optional reference SnpMatrix containing data with the desired
#' allele frequencies to impute to, or otherwise will use the frequencies of 
#' the current dataset in snpMatLst.
#' @param filename base file name prefix for the objects on disk that will store
#' the BigSnpMatrix object, include a path if not wanting to write to the working
#' directory.
#' @param tracker logical, whether to show a progress bar during processing
#' @param limit.ram logical, whether to regularly flush RAM during processing
#' of large objects (can create some long pauses during processing, but results
#' in more conservative, stable, memory management). Usually works ok when set FALSE.
#' @param mn.zero logical, set TRUE to mean correct the values in the dataset; if
#' ref.data is in use, will correct using the SNP means in that dataset, if not,
#' then will correct to zero (using the means in the snpMatLst dataset). For PCA
#' you would likely set this option to TRUE.
#' @param sd.hwe logical, set TRUE to SD correct the values in the dataset; if
#' ref.data is in use, will correct using the SNP means in that dataset, if not,
#' then will correct using the means in the snpMatLst dataset. Either way this
#' correction is done by using the expected standard deviation of the reference
#' allele frequency at each SNP if it were in perfect Hardy Weinburg Equilibrium.
#' For PCA you would often set this option to TRUE.
#' @return returns a BigSnpMatrix object which is essentially a big.matrix object
#' from the bigmemory() package, where genotypes are coded 1,2,3 (aa, ab, bb) where
#' b is the reference allele and zero is missing. However if mn.zero or sd.hwe are
#' TRUE, then the values will be post correction, so will not be such neat integer
#' values for each genotype. PCA can be performed on this object using the bigpca package.
#' @export
#' @seealso snp.from.annot, chip.support, bigpca
#' @examples 
#' snp.set <- sample(rownames(chip.support()),200) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(50,200); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,])
#' big.mat <- bigSnpMatrixList(sml,mn.zero=FALSE,sd.hwe=FALSE) # genotype values are 1,2,3
#' big.mat <- bigSnpMatrixList(sml) # genotype values are now continuous, and ready for PCA
bigSnpMatrixList <- function(snpMatLst, verbose=TRUE, dir=getwd(),replace.missing=TRUE, ref.data=NULL,
                             filename="tempMatrix",tracker=TRUE, limit.ram=FALSE, mn.zero=TRUE, sd.hwe=TRUE) {
  typ <- snp.mat.list.type(snpMatLst,fail=TRUE) #disk or memory
  typ2 <- get.split.type(snpMatLst,dir=dir)  # singleEntry snpSubGroups sampleSubGroups
  max.gb <- NA
  if(typ2=="snpmatrix") { return(bigSnpMatrix(snpMatLst)) } # was already a SnpMatrix
  if(typ2=="singleEntry") { 
    # just 1 element, no binding required
    if(typ=="disk") {
      oo <- (get.SnpMatrix.in.file(snpMatLst[[1]],dir=dir))
    } else {
      oo <- (snpMatLst[[1]])
    }
    return(bigSnpMatrix(oo))
  }
  spec <- get.snpMat.spec(snpMatLst,dir=dir)
  nR <- nrowL(snpMatLst,dir=dir,list=F)
  nC <- ncolL(snpMatLst,dir=dir,list=F)
  cN <- colnamesL(snpMatLst,dir=dir,list=F)
  rN <- rownamesL(snpMatLst,dir=dir,list=F)
  #if(!check.snpmat.size(nr,nc)) { stop("Proposed bigSnpMatrix would be",nR,"rows by",nC,"columns, which exceeds the maximum possible size of a standard R object, suggest leaving this object as a SnpMatrixList") }
  dir2 <- dirname(filename); filename <- basename(filename)
  if(dir2=="") { dir2 <- getwd() }
  if(tracker) { cat(" creating",nR,"x",nC,"target matrix,",cat.path(dir,filename),"...\n") }
  des <- paste(filename,"dsc",sep=".")
  bck <- paste(filename,"bck",sep=".")
  bigSnp <- big.matrix(nrow=nR,ncol=nC, backingfile=bck, dimnames=list(rN,cN),
                       backingpath=dir2, descriptorfile=des)
  #typ2 <- "snpSubGroups" # why???
  if((typ2 %in% c("sampleSubGroups","snpSubGroups"))) {
    if(verbose) { cat("Extracting",length(snpMatLst),if(typ2=="sampleSubGroups") {"sample"} else {"SNP"},"subsets from files\n") }
    if(typ2=="snpSubGroups") { 
      stz <- c(1,1+cumsum(spec[2,])[-ncol(spec)]); enz <- cumsum(spec[2,])  # snp index
    } else {
      stz <- c(1,1+cumsum(spec[1,])[-ncol(spec)]); enz <- cumsum(spec[1,])  # sample index
    }
    #### HERE###
    n.chr <- length(snpMatLst)
    for (dd in 1:n.chr) {
     # print(snpMatLst[[dd]])
      if(typ=="disk") {
        next.chr <- (get.SnpMatrix.in.file(snpMatLst[[dd]],dir=dir))
      } else {
        next.chr <- snpMatLst[[dd]]
      }
      ### check type of next chromosome (or sample subset)
      typ3 <- is(next.chr)[1]
      if(typ3 %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) {
        if(replace.missing) { next.chr <- randomize.missing2(next.chr) }
        raw <- next.chr@.Data
      } else {
        stop("SnpMatLst must contain only SnpMatrix or aSnpMatrix objects", "; type was: ",typ3)
      }
      #### means and SDs for standardizing ####
      if(mn.zero | sd.hwe) {
        if(!is.null(ref.data)) {
          oo <- correct.by(ref.data[,colnames(next.chr)]); mN <- oo$mean; sD <- oo$sd
        } else {
          oo <- correct.by(next.chr); mN <- oo$mean; sD <- oo$sd
        }
      }
      if(!mn.zero) { mN <- 0 }
      if(!sd.hwe) { sD <- 1 }
      #########
      split.to <- round(10*estimate.memory(next.chr)) # split into .1GB chunks, save RAM without creating groups too small to process
      #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
      nr <- nrow(next.chr); nc <- ncol(next.chr)
      if(length(mN)==1) { mN <- rep(mN,nc) } # if there is a single mean correction factor, extend to number of SNPs
      if(length(sD)==1) { sD <- rep(sD,nc) } # same for SD correction factor
      stepz <- round(seq(from=1,to=nc+1,length.out=round((split.to+1))))
      if((tail(stepz,1)) != nc+1) { stepz <- c(stepz,nc+1) }
      split.to <- length(stepz)-1
      for (cc in 1:split.to)
      {
        # within submatrix cols
        c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
        # do the copying
        lilColRange <- c(c1:c2)
        if(is.finite(sum(as.numeric(lilColRange)))) {
          #cat(range(lilColRange)); cat(dim(bigTrans)); cat(dim(bigMat))
          nxt.chunk <- raw[1:nr,lilColRange,drop=F]
          nxt.chunk <- t(apply(nxt.chunk,2,as.numeric))
          if(typ2=="snpSubGroups") {
#            prv(stz[dd],enz[dd])
            #return(list(bigSnp=bigSnp,stz=stz,enz=enz,dd=dd,nr=nr,lilColRange=lilColRange,nxt.chunk=nxt.chunk,mN=mN,sD=sD))
            bigSnp[1:nr,(stz[dd]:enz[dd])[lilColRange]] <- t((nxt.chunk-mN[lilColRange])/sD[lilColRange]  )
          } else {
            bigSnp[stz[dd]:enz[dd],lilColRange] <- t((nxt.chunk-mN[lilColRange])/sD[lilColRange])
          }
        } else {
          cat(" Warning: empty interval ignored\n")
        }
        if(limit.ram & (cc %% (max.gb*10) == 0)) {
          # reset memory after every 'max.gb' 1GB chunks to prevent skyrocketing RAM use #
          fl.suc <- bigmemory::flush(bigSnp) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()  
          if(T) {
            RR <- describe(bigSnp); rm(bigSnp); bigSnp <- attach.big.matrix(RR,path=dir)
          }
        }
      }
      if(tracker) { loop.tracker(dd,n.chr) }
    }
  } else {
    stop("incorrect input caused failure to convert from SnpMatList to bigSnpMatrix format")
  }
  if(tracker) { prv(bigSnp) }
  return(bigSnp)
}

#ld2 <- function(X) {
#  mat <- SnpMatrix.to.data.frame(X,T)
#  kk <- proc.time()[3]; out <- cor(mat,use="pairwise.complete")^2 ;jj <- proc.time()[3]; print(jj-kk)
#  return(out)
#}



#' Convert a SnpMatrix to a SnpMatrixList split by chromosomes
#'
#' SnpMatrixList objects allow use of functions that apply to
#' SnpMatrix objects, but with storage of the dataset in a list
#' (by subsample or chromosome for instance) preventing either
#' reaching the maximum allowed size of a SnpMatrix object,
#' or reducing the RAM footprint. This function converts a regular
#' SnpMatrix object into a SnpMatrixList, split by chromosome.
#' @param X a SnpMatrix, XSnpMatrix, aSnpMatrix or aXSnpMatrix object
#' @param snp.info optional annotation object giving the chromosome and
#' position information for the SNPs in X, should contain at least all of
#' the colnames of X as rownames in snp.info. Should be GRanges, RangedData
#' or ChipInfo type. If X is an aSnpMatrix, then you should leave snp.info as
#' NULL and this function will extract the snp information from the aSnpMatrix.
#' Similarly if you are working with SNPs present in the current chip.support()
#' object, then positions will automatically be derived from this support object
#' when snp.info=NULL.
#' @param dir either 'NULL' to generate a memory based SnpMatrixList (full list is
#' held in the R workspace), or provide an existing directory to store the binary
#' files if you would like the SnpMatrixList to be a disk based list (data stored
#' on the hard drive and loaded only when required).
#' @param build the ucsc build in use, e.g, "hg19", however leave NULL if you are
#' working with the build = getOption("ucsc"). It is recommended when using the
#' humarray package to always set this option to the ucsc build you are working
#' with.
#' @param genomeOrder logical, TRUE to force the list to store the matrices
#' in genome order (e.g, sort the chromosomes), or FALSE leaves the order as inputted.
#' @return returns a SnpMatrixList object which if dir=NULL will be a large list of
#' SnpMatrix objects, or if dir is valid, will save the data to disk and the object
#' returned is compact, containing links to data stored on disk.
#' @export
#' @seealso snp.from.annot, chip.support
#' @examples 
#' snp.set <- sample(rownames(chip.support()),200) # select a random set of 200 SNPs
#' snp.mat <- rSnpMatrix(50,200); colnames(snp.mat) <- snp.set # create a random SnpMatrix
#' sml <- SnpMatrix.to.sml.by.chr(snp.mat,chip.support()[snp.set,])
SnpMatrix.to.sml.by.chr <- function(X,snp.info=NULL,dir=NULL,build=NULL,genomeOrder=TRUE) {
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { stop("X must be a SnpMatrix or equivalent") }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(snp.info) & (typ %in% c("aSnpMatrix","aXSnpMatrix")))  { snp.info <- snp.from.annot(X) }
  if(is.null(snp.info)) { snp.info <- chip.support(build=build) }
  if(any(!colnames(X) %in% rownames(snp.info))) { stop("all SNPs from X must have support in the aSnpMatrix, snp.info, or the result of chip.support(build=build)")}
  snp.info <- snp.info[colnames(X),]
  if(!is(snp.info)[1] %in% c("GRanges","RangedData","ChipInfo")) { stop("snp.info must be GRanges, RangedData or ChipInfo") }
  if(genomeOrder) {
    snp.info <- toGenomeOrder(snp.info)
    X <- X[,rownames(snp.info)]
  }
  ch <- paste(chrm(snp.info))
  chrz <- unique(paste(chrm(snp.info)))
  n.chr <- length(chrz)
  sml <- vector("list",n.chr)
  if(!is.null(dir)) { if(file.exists(dir)) { disk <- TRUE } else { stop("dir '",dir,"'' did not exist") } } else { disk <- FALSE }
  for (cc in 1:n.chr) {
    if(disk) { 
      ofn <- cat.path(dir,"CHR",suf=chrz[cc],ext="RData")
      var <- paste0("chr",chrz[cc])
      x <- X[,ch %in% chrz[cc]]
      save(x,file=ofn)
      sml[[cc]] <- ofn
    } else {
      sml[[cc]] <- X[,ch %in% chrz[cc]]
    }
  }
  return(sml)
}


# a vignette?
# to PCA 1000 genomes for available iChip Snps
# # assume already extracted 1000 genomes for ichip and aligned using annotSnpStats to our dataset (T1D-ichip)
# (load("/chiswick/data/ncooper/imputation/THOUSAND/aligned1000g.RData"))
# newsml <- SnpMatrix.to.sml.by.chr(asm.1000g.aligned)
# ld.pruned <- smlapply(newsml,ld.prune.chr,thresh=0.1,n.cores=12) # this step is slow (several hours?)
# pruned.1000g <- NULL
# pruned.1000g <- asm.1000g.aligned[,names(unlist(ld.pruned))]
# save(ld.pruned,pruned.1000g,file="pruned.list1000g.RData")
# (load("pruned.list1000g.RData"))
# (load("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData"))
# big1000.A <- bigSnpMatrix(pruned.1000g,"big1000gA")
# big1000.B <- bigSnpMatrix(pruned.1000g,"big1000gB",mn.zero=TRUE,sd.hwe=FALSE)
# result.quick.A <- big.PCA(big.t(big1000.A),return.loadings=TRUE)
# result.quick.B <- big.PCA(big.t(big1000.B),return.loadings=TRUE)
# anc <- factor(sample.info$ancestry[match(rownames(result.quick.A$PCs),rownames(sample.info))])
# pdf("PCA1000.A.pdf"); plot(result.quick.A$PCs[,"PC1"],result.quick.A$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
# legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()
# pdf("PCA1000.B.pdf"); plot(result.quick.B$PCs[,"PC1"],result.quick.B$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
# legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()
# use snpSel to get all the desired SNP data for the ancestry SNPs
# # top5pc.A <- (which(abs(result.quick.A$loadings[,1])>.01234))
# # divide 1 by 2 of vec: estimate.eig.vpcs(result.quick.A$Evalues,M=big.t(big1000.A))$variance.pcs[1:2]
# top5pc.A <- rev(order(abs(result.quick.A$loadings[,1])+1.47*abs(result.quick.A$loadings[,2])))[1:10000]
# # top5pc.B <- (which(abs(result.quick.B$loadings[,1])>.013295))
# top5pc.B <- rev(order(abs(result.quick.A$loadings[,1])+1.27*abs(result.quick.A$loadings[,2])))[1:10000]
# pca1.pred <- scale(big1000.A[,top5pc.A],center=T,scale=F) %*% result.quick.A$loadings[top5pc.A,1]
# pca2.pred <- scale(big1000.A[,top5pc.A],center=T,scale=F) %*% result.quick.A$loadings[top5pc.A,2]
# cor(pca1.pred,result.quick.A$PCs[,1],use="pairwise.complete")
# # [1,] 0.9981955
# cor(pca2.pred,result.quick.A$PCs[,2],use="pairwise.complete")
# # [1,] 0.9992125
# cor(pca2.pred,result.quick.A$PCs[,2],use="pairwise.complete")
# pdf("pca1.pdf"); plot(pca1.pred,result.quick.A$PCs[,1]); dev.off()
# top5pc.An <- colnames(big1000.A)[top5pc.A]
# top5pc.Bn <- colnames(big1000.B)[top5pc.B]
# save(top5pc.An,top5pc.Bn,result.quick.A,result.quick.B,anc,sample.info,pruned.1000g,file="pruned.list.PCA.1000g.RData")
### GET all T1D data ###
# dd <- ("/chiswick/data/ncooper/imputation/T1D/PQDATA/")
# sml <- list.files(dd)
# sml <- as.list(sml)
# sml <- sml[c(18:19,26:40,1:17,20:24)]
# SM <- snpSel(sml,snps=top5pc.An,dir=dd)
# BSM.A <- bigSnpMatrix(SM,"bigPCAa",ref.data=pruned.1000g)
# SM2 <- snpSel(sml,snps=top5pc.Bn,dir=dd)
# BSM.B <- bigSnpMatrix(SM2,"bigPCAb",mn.zero=TRUE,sd.hwe=FALSE,ref.data=pruned.1000g)
# pca.bsm.a.pred.1 <- scale(BSM.A[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.A,1]
# pca.bsm.a.pred.2 <- scale(BSM.A[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.A,2]
# pca.bsm.b.pred.1 <- scale(BSM.B[,top5pc.Bn],center=T,scale=F) %*% result.quick.B$loadings[top5pc.Bn,1]
# pca.bsm.b.pred.2 <- scale(BSM.B[,top5pc.Bn],center=T,scale=F) %*% result.quick.B$loadings[top5pc.Bn,2]
# pdf("PCAT1D.A.pdf"); plot(result.quick.A$PCs[,"PC1"],result.quick.A$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
# points(pca.bsm.a.pred.1,pca.bsm.a.pred.2,col="black")
# legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()
# pdf("PCAT1D.B.pdf"); plot(result.quick.B$PCs[,"PC1"],result.quick.B$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
# points(pca.bsm.b.pred.1,pca.bsm.b.pred.2,col="black")
# legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()


#' Select chromosome subset for aSnpMatrix
#' 
#' Returns the object filtered for specific chromosomes for a aSnpMatrix object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "aSnpMatrix", function(object,chr) {
  si <- snp.info.from.annot(object)
  subs <- rownames(humarray::chrSelect(si,chr))
  return(object[,subs])
})
 
#' Select chromosome subset for aSnpMatrix
#' 
#' Returns the object filtered for specific chromosomes for a aXSnpMatrix object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "aXSnpMatrix", function(object,chr) {
  si <- snp.info.from.annot(object)
  subs <- rownames(humarray::chrSelect(si,chr))
  return(object[,subs])
})


#' Linkage disequilibrium reduction of a SNP-set
#' 
#' For purposes of data reduction, or principle components
#' analysis for detecting ancestry, it is common practice to 'LD prune'
#' a microarray genotype SNP dataset from a larger to a smaller set.
#' The method involves calculating the matrix of linkage disequilibrium
#' weights and then selecting the strongest signals, removing the correlates
#' of each.
#' @param X X should be a SnpMatrix or aSnpMatrix object
#' @param stats select which kind of LD calculate to use, see arguments for the SnpStats function 
#' 'ld()'. Values can be: 'R.squared', 'D.prime', or 'R'
#' @param thresh numeric, correlation/statistical threshold to use as a cutoff, above which
#' all correlates of the strongest signal in each region will be discarded/
#' @return vector of integer indexes of the set of SNPs in X to keep, based on the LD statistic
#' and threshold selected
#' @export
#' @examples
#' # ownexample
#' data(exSnpMat)
#' res <- ld.prune.chr(exSnpMat)
#' small.dat <- exSnpMat[,res]
#' rest.dat <- exSnpMat[,-res]
#' LD <- ld(small.dat,small.dat,stats="R.squared")
#' print(LD)
#' LD <- ld(rest.dat,small.dat,stats="R.squared")
#' print(LD)
#' apply(LD,1,max) # each SNP excluded has a correlation >.1 with one of the two kept
ld.prune.chr <- function(X,stats="R.squared",thresh=.1) {
  if(estimate.memory(ncol(X)^2)>=0.1) { ld.prune.big(X,stats=stats,thresh=thresh) } 
  kk <- proc.time()[3]; mat <- ld(X,X,stats=stats[1]) ; jj <- proc.time()[3]; print(jj-kk)
  matt <- mat; matt[is.na(matt)] <- 0
  exc <- length(matt[matt>thresh])
  nsnp <- nrow(mat)
  cat(out.of(exc,nsnp^2)," in matrix exceed ld threshold\n")
  kk <- proc.time()[3];
  mat[is.nan(mat)] <- runif(length(which(is.nan(mat))))/100
  selected <- rep(T,nrow(mat))
  winners <- NULL
  while(any(selected)) {
    tot.scores <- abs(rowMeans(mat[,selected,drop=F],na.rm=T))
    tot.scores[is.na(tot.scores)] <- 0
    winnerz <- which(tot.scores==max(tot.scores[selected],na.rm=T))
    winnerz <- winnerz[winnerz %in% which(selected)]
    if(length(winnerz)<1) {
      winnerz <- which(tot.scores==max(tot.scores[selected],na.rm=T))
      winnerz <- winnerz[winnerz %in% winners]
      cat("this happened")
    }
    winner <- winnerz[sample(length(winnerz),1)]
    if(!is.na(winner) & (!winner %in% winners)) {
      winners <- c(winners,winner)
      passz <- (mat[,winner] < thresh)
      passz[is.na(passz)] <- T
      #cat("winner ",rownames(mat)[winner],"had",length(which(!passz)),"in ld. ",out.of(length(which(selected)),length(selected)),"\n")
      selected <- selected & passz
      if(any(is.na(selected))) { prv(selected,passz)}
    } 
    selected[winner] <- F
  }
  jj <- proc.time()[3]; print(jj-kk)
  cat("pruned",out.of((nsnp-length(winners)),nsnp),"SNPs\n")
  return(winners)
}


#' Split an aSnpMatrix into separate files for each chromosome
#'
#' When analysing large datasets it is common to split the datafile
#' into 22 separate files by chromosome to allow parallel processing
#' of smaller subsets. This function
#' takes an aSnpMatrix object from annotSnpStats and saves a file for
#' each chromosome, return a list of filenames (SnpMatrixList).
#' @param aSnpMatrix an aSnpMatrix object
#' @param build string, ucsc annotation build, default is 37/hg19
#' @param fname character, prefix to use for the name of the output .RData files.
#' Can include a directory.
#' @param verbose logical, whether to display information about progress
#' @param tracker logical, whether to display a progress bar (only if verbose=FALSE)
#' @return returns a list of file names for each chromosome.
#' @export
#' @examples
#' #ownexample
#' data(exAnnotSnp)
#' prv(exAnnotSnp)
#' outlist <- split.chr(exAnnotSnp)
#' prv(outlist)
split.chr <- function(aSnpMat,build=37,fname="SnpMat",verbose=FALSE,tracker=TRUE) {
  si <- snp.from.annot(aSnpMat)
  if(!is(si)[1]=="RangedData") { si <- as(si,"RangedData") }
  if(!is(si)[1]=="RangedData") { stop("couldn't extract snp.info from aSnpMat") }
  chrz <- chrNames(si)
  chrz2 <- gsub("chr","",paste(chrz),ignore.case = TRUE)
  if(any(rownames(si)!=colnames(aSnpMat))) { stop("mismatch between aSnpMat names and annotation") }
  sml <- vector("list",length(chrz))
  for(cc in 1:length(chrz)) {
    nxt.chr <- chrSel(si,paste(chrz[cc]))
    posz <- start(nxt.chr)
    indx <- match(rownames(nxt.chr),rownames(si))
    if(length(indx)<1) { warning("found no SNPs in chromosome ",chrz[cc])}
    X <- aSnpMat[,indx]
    sml[[cc]] <- ofn <- cat.path(dir="",paste(fname[1]),suf=chrz[cc],ext="RData")
    save(X,file=ofn)
    if(verbose) { cat("wrote chr",chrz[cc],"to",ofn,"\n") } else {
      if(tracker) {
        loop.tracker(cc,length(chrz))
      }
    }
  }
  return(sml)
}



#' Split an aSnpMatrix into a long and short chromosome arm
#'
#' When analysing large datasets it is common to split the datafile
#' into 22 separate files by chromosome to allow parallel processing
#' of smaller subsets. Sometimes even with this split a dataset may
#' be too large, particularly for the larger chromosomes. Another
#' step that can be taken is to split further by chromosome 'arm'.
#' Labelled as the long and short arms (arm 'p' and arm 'q'), these
#' are the portions of the chromosome separated by the centromere,
#' with biological justification for separation as interaction and LD
#' is far less likely across the centromere. This function
#' takes an aSnpMatrix object from annotSnpStats and splits into a
#' separate p and q arm. Note for many chips it's common to have
#' SNPs for only 1 arm, so allow for this in your code.
#' @param aSnpMatrix an aSnpMatrix object
#' @param build string, ucsc annotation build, default is 37/hg19
#' @param pqvec logical, if TRUE simply return a character vector of 
#' p's and q's indicating which arm for each SNP in order. If FALSE,
#' will return whole aSnpMatrix objects.
#' @param verbose logical, whether to display information about progress
#' @return If pqvec is TRUE simply return a character vector of 
#' p's and q's indicating which arm for each SNP in order. If FALSE,
#' will return two whole aSnpMatrix objects in a named list (arm p and q).
#' @export
#' @examples
#' #ownexample
#' data(exAnnotSnp)
#' prv(exAnnotSnp)
#' outlist <- split.pq(exAnnotSnp)
#' prv(outlist)
split.pq <- function(aSnpMat,build=37,pqvec=FALSE,verbose=TRUE) {
  si <- snp.from.annot(aSnpMat)
  if(!is(si)[1]=="RangedData") { si <- as(si,"RangedData") }
  if(!is(si)[1]=="RangedData") { stop("couldn't extract snp.info from aSnpMat") }
  cent <- get.centromere.locs(build=build)
  cnt.ch <- chrm(cent)
  chrz <- chrNames(si)
  chrz2 <- gsub("chr","",paste(chrz),ignore.case = TRUE)
  arm <- rep("p",nrow(si))
  for(cc in 1:length(chrz)) {
    ii <- which(cnt.ch==chrz2[cc])
    if(length(ii)>0) { 
      if(length(ii)>1) {
        warning("more than one centromere entry for a single chromosome - import failure likely")
      }
      nxt.chr <- chrSel(si,paste(chrz[cc]))
      cent.loc <- start(cent[ii[1],])
      if(verbose) { cat("centromere location:",cent.loc,"\n") }
      posz <- start(nxt.chr)
      indx <- match(rownames(nxt.chr),rownames(si))
      #prv(cent.loc,posz)
      arm[indx][posz > cent.loc] <- "q"
    } else {
      if(!chrz2[cc] %in% c("XY","MT","M")) {
        warning("centromere not found for chromosome labelled:",chrz2[cc])
      }
    }
  }
  #bb <- Band(si,build=37)
  #arm <- rep("p",nrow(si))
  #arm[grep("q",bb)] <- "q"
  if(pqvec) { return(arm) }
  if(length(unique(arm))==1) { warning("only arm values of ",unique(arm)," were found"); return(aSnpMat) }
  mat1 <- aSnpMat[,arm=="p"]
  mat2 <- aSnpMat[,arm!="p"]
  if(verbose) { cat("split aSnpMatrix into",ncol(mat1),"on arm 'p' and",ncol(mat2),"on arm 'q'\n")  }
  out <- list(arm.p=mat1,arm.q=mat2)
  return(out)
}

#' Create aSnpMatrix from SnpMatrix plus SNP/sample support files
#' 
#' @param snpMat a SnpMatrix or XSnpMatrix
#' @param snp.info annotation with rownames as snp labels, then chromosome, position and allele 
#' codes, can be RangedData, GRanges, ChipInfo or data.frame type
#' @param sample.info data.frame with rownames of sample ids, may also contain 'sex' and some
#' phenotype variable, such as 'pheno' or 'affected' or 'case' (autodetected), may also contain
#' pedigree information, such as a plink .fam/ped file
#' @param snp.excl text filename or character vector containing SNP ids to exclude
#' @param samp.excl text filename or character vector containing sample ids to exclude
#' @export
#' @return returns an annotSnpStats::aSnpMatrix object with @samples and @snps info taken
#' from the inputted support files. Will only include samples and SNPs occuring in both
#' the dataset and the support files.
#' @examples
#' ## my ownexample for now, only works @DIL ##
#' t1.dr <- "/chiswick/data/ncooper/imputation/T1D/"
#' setwd(t1.dr)
#' smp <- reader("/chiswick/data/ncooper/iChipData/sample.info.RData")
#' snp.excl <- "/chiswick/data/ncooper/imputation/T1D/snpsExcluded.txt"
#' smp.excl <- "/chiswick/data/ncooper/imputation/T1D/sampsExcluded.txt"
#' print(load("CHR20.RData"))
#' x <- gt.for.impute
#' si <- sup.for.impute
#' aSnpMat <- annot.sep.support(x,si,smp, snp.excl=snp.excl, samp.excl=smp.excl)
annot.sep.support <- function(snpMat, snp.info, sample.info, snp.excl=NULL, 
                              samp.excl=NULL, warn.int=TRUE, genome.order=TRUE, verbose=TRUE) {
  #### check for valid input types ####
  if(!is(snpMat)[1] %in% c("XSnpMatrix","SnpMatrix")) { 
    if(snp.mat.list.type(snpMat)!="error") {
      snpMat <- snpMatLst.collate(snpMat)
    } else {
      stop("snpMat must be of type 'SnpMatrix' or 'XSnpMatrix' or a 'snpMatrixList'") 
    }
  }
  if(is(snp.info)[1] %in% c("RangedData","GRanges","ChipInfo")) { 
    snp.info <- as(snp.info,"GRanges") 
    ii <- grep("elementMetadata.",colnames(mcols(snp.info)))
    if(length(ii)>0) {
      colnames(mcols(snp.info)) <- gsub("elementMetadata.","",colnames(mcols(snp.info)))
    }
  } else { 
    snp.info <- as(snp.info,"data.frame")
    snp.info <- column.salvage(snp.info,"pos",c("Start","Pos","pos","POS","Pos37","pos37","Pos36","pos36","position320"),ignore.case=FALSE)
    snp.info <- column.salvage(snp.info,"chr",c("chr","chromosome","chrom","chromo","space","seqnames"),ignore.case=TRUE)
    #snp.info <- column.salvage(snp.info,"end",c("end"),ignore.case=TRUE)
    print(head(snp.info))
    snp.info <- df.to.ranged(snp.info,GRanges=TRUE)
  }
  if(!is(snp.info)[1] %in% c("GRanges")) { stop("snp.info must be coercible to a data.frame") }
  if(genome.order) { snp.info <- toGenomeOrder(snp.info) }
  sample.info <- force.frame(sample.info)
  if(!is(sample.info)[1] %in% c("data.frame")) { stop("sample.info must be coercible to a data.frame") }
  ###############################
  # find the right columns in the sample file
  sample.info <- column.salvage(sample.info,"affected",c("affected","phenotype","pheno","case","cases","phenot","phen","aff","ph","controls","ctrls","ctrl","control"),ignore.case=TRUE)
  # extract optional fields
  suppressWarnings({
    sample.info <- column.salvage(sample.info,"sex",c("sex","gender","mf"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"father",c("father","paternal","pat"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"mother",c("mother","maternal","mat"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"pedigree",c("pedigree","family","fam"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"id",c("id","sampleid","caseid","member"),ignore.case=TRUE)
    sample.info <- column.salvage(sample.info,"plate",c("plate","array.plate","plateid","plate.id"),ignore.case=TRUE)
  })
  if(!is.null(snp.excl)) { snp.excl <- force.vec(snp.excl) }
  if(!is.null(samp.excl)) { samp.excl <- force.vec(samp.excl) }
  ## now can assume valid input parameters have been checked ##
  samps.in.dat <- rownames(snpMat)
  snps.in.dat <- colnames(snpMat)
  samps.in.support <- rownames(sample.info)
  snps.in.support <- rownames(snp.info)
  if(warn.int) {
    # warn when sample/snp support doesn't perfectly match the SnpMatrix row/colnames #
    if(!all(samps.in.dat %in% samps.in.support)) { warning("sample list in support differs from data, will take intersection")}
    if(!all(snps.in.dat %in% snps.in.support)) { warning("SNP list in support differs from data, will take intersection")}
  }
  final.samps <- samps.in.dat[(samps.in.dat %in% samps.in.support) & (!samps.in.dat %in% samp.excl)]
  final.snps <- snps.in.dat[(snps.in.dat %in% snps.in.support) & (!snps.in.dat %in% snp.excl)]
  if(length(final.snps)<1) { warning("after filtering, no SNPs remain") ; return(NULL) }
  if(length(final.samps)<1) { warning("after filtering, no samples remain") ; return(NULL) }
  if(genome.order) { snp.ord <- rownames(snp.info); final.snps <- snp.ord[snp.ord %in% final.snps] }
  pre <- Dim(snpMat)
  snpMat <- snpMat[final.samps,final.snps]
  post <- Dim(snpMat)
  if(verbose) { cat("snpMat original size: ",comma(pre),"; new size:",comma(post),"\n") }
  pre <- Dim(sample.info) ; sample.info <- sample.info[final.samps,] ;  post <- Dim(sample.info)
  if(verbose) { cat("sample.info original size: ",comma(pre),"; new size:",comma(post),"\n") }
  #prv(snp.info)
  pre <- Dim(snp.info) ; snp.info <- snp.info[final.snps,] ;  post <- Dim(snp.info)
  if(verbose) { cat("snp.info original size: ",comma(pre),"; new size:",comma(post),"\n") }
  #################################
  ##### CREATE THE @SNPS SLOT #####
  # find alleles, and other support in the snp.info object
  df.info <- ranged.to.data.frame(snp.info,include.cols=TRUE)
  #
  df.info <- column.salvage(df.info,"allele.1",c("allele.1","allele_1","allele-1","a1","allele1",
                                                 "allele-A","alleleA","allele.a","allele_a","forward","fwd","allele.f","allele.fwd"),ignore.case=TRUE)
  df.info <- column.salvage(df.info,"allele.2",c("allele.2","allele_2","allele-2","a2","allele2",
                                                 "allele-B","alleleB","allele.b","allele_b","backward","bk","bck","allele.bck"),ignore.case=TRUE)
  # extract optional fields
  suppressWarnings({
    df.info <- column.salvage(df.info,"cM",c("CM","distance","map"),ignore.case=TRUE);
    df.info <- column.salvage(df.info,"snp.name",c("snp.name","snpid","snp.id","dbSNP","id","label","rs.id","rsid"),ignore.case=TRUE)
  })
  snp.rn <- rownames(snp.info)
  snp.ch <- chrm(snp.info)
  if("snp.name" %in% colnames(df.info)) { snp.nm <- df.info[,"snp.name"] } else { snp.nm <- snp.rn }
  if("cM" %in% colnames(df.info)) { snp.cm <- df.info[,"cM"] } else { snp.cm <- rep(NA,length(snp.rn)) }
  snp.ps <- start(snp.info)
  if("allele.1" %in% colnames(df.info)) { snp.a1 <- df.info[,"allele.1"] } else { snp.a1 <- rep(NA,length(snp.rn)) }
  if("allele.2" %in% colnames(df.info)) { snp.a2 <- df.info[,"allele.2"] } else { snp.a2 <- rep(NA,length(snp.rn)) }
  if(!all(length(snp.rn)==c(length(snp.ch),length(snp.nm),length(snp.cm),length(snp.ps),length(snp.a1),length(snp.a2)))) {
    stop("Import of snp.info failed, columns had different lengths") # this should be almost impossible to happen HERE
  } else {
    if(verbose) { cat("@snps slot created successfully\n") }
  }
  snps.slot <- data.frame(chromosome=snp.ch,snp.name=snp.nm,cM=snp.cm,position=snp.ps,allele.1=snp.a1,allele.2=snp.a2)
  rownames(snps.slot) <- snp.rn
  allele.slot <- c("allele.1","allele.2")
  #################################
  #### CREATE THE @SAMPLES SLOT ###
  smp.rn <- rownames(sample.info)
  if("pedigree" %in% colnames(sample.info)) { smp.pd <- sample.info[,"pedigree"] } else { smp.pd <- smp.rn }
  if("member" %in% colnames(sample.info)) { smp.mb <- sample.info[,"member"] } else { smp.mb <- smp.rn }
  if("father" %in% colnames(sample.info)) { smp.ft <- sample.info[,"father"] } else { smp.ft <- rep(NA,length(smp.rn)) }
  if("mother" %in% colnames(sample.info)) { smp.mt <- sample.info[,"mother"] } else { smp.mt <- rep(NA,length(smp.rn)) }
  if("sex" %in% colnames(sample.info)) { smp.sx <- sample.info[,"sex"] } else { smp.sx <- rep(NA,length(smp.rn)) }
  if("plate" %in% colnames(sample.info)) { smp.pl <- sample.info[,"plate"] } else { smp.pl <- rep(NA,length(smp.rn)) }
  smp.ph <- sample.info[,"affected"]
  if(!all(length(smp.rn)==c(length(smp.pd),length(smp.mb),length(smp.ft),length(smp.mt),length(smp.sx),length(smp.ph),length(smp.pl)))) {
    stop("Import of sample.info failed, columns had different lengths") # this should be almost impossible to happen HERE
  } else {
    if(verbose) { cat("@samples slot created successfully\n") }
  }
  samps.slot <- data.frame(pedigree=smp.pd,member=smp.mb,father=smp.ft,mother=smp.mt,sex=smp.sx,affected=smp.ph) #,plate=snp.pl)
  rownames(samps.slot) <- smp.rn
  pheno.slot <- c("affected")
  #################################
  
  ## aSnpMatrix format SPECIFICATION: ##
  #@.Data a SnpMatrix
  #@snps  
  #              chromosome     snp.name cM position allele.1 allele.2
  # imm_1_898835          1 imm_1_898835 NA   898835     <NA>        A
  #@samples
  #          pedigree             member father mother sex affected
  # 2447     2447 120434_A01_BLOOD320736     NA     NA   1        1
  #@phenotype character() [phenotype colname]
  #@alleles = c("allele.1", "allele.2") [allele colnames]
  #######################################
  aSnpMat <- new("aSnpMatrix", .Data = snpMat, snps = snps.slot, 
                 samples = samps.slot, alleles = allele.slot, 
                 phenotype = pheno.slot)
  return(aSnpMat)
}



#########################
## SnpMatrix Functions ##
#########################


#' Calculate Lambda inflation factors for SNP dataset
#' 
#' This function calculates SNP-wise or overall Lambda and Lambda_1000 statistics for inflation due
#' to population structure. It works on a SnpMatrix object or dataframe coded 0,1,2,NA (autodetects
#'  which).
#' @param X SnpMatrix or data.frame coded 0,1,2,NA containing genotypes
#' @param pheno vector of coding for phenotypes, using (0,1) or (1,2) for controls, cases, should 
#' match nrow(X)
#' @param checks logical, whether to perform some sanity checks that will slightly slown down 
#' performance
#' @param output character vector, containing any of the following "all", "lambda","l1000", "both";
#' if snp.wise is false, the scalar value(s) specified by 'output', or otherwise a matrix
#' of parameters for each SNP, optionally limited to just the lambda column(s) by the value of 
#' 'output'.
#' @param snp.wise logical, if TRUE, return lambda statistics separately for each SNP (rather than
#'  an overall)
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @seealso \code{\link{lambda_nm}}
#' @examples
#' # http://en.wikipedia.org/wiki/Population_stratification
#' # Note that Y2 ~ L*X2, where allele counts are symbolized:
#' # Case     r0  r1  r2  R
#' # Control  s0  s1  s2  S
#' # Total    n0  n1  n2  N
#' require(snpStats) ; data(testdata)
#' pheno <- rep(0,nrow(Autosomes))
#' pheno[subject.data$cc=="case"] <- 1
#' lambdas(Autosomes,pheno)
#' lambdas(Autosomes[,1:5],pheno,snp.wise=TRUE) # list everything snp-wise for the first 5 SNPs
#' lambdas(Autosomes[,1:5],pheno) # just for the first 5 SNPs
lambdas <- function(X, pheno, checks=TRUE, 
                    output=c("all","lamba","l1000","both"),snp.wise=FALSE) {
  ## i don't think this should ever be used: ?? cc1000=TRUE, ??
  cc1000 <- FALSE
  ## workhorse internal function ##
  do.lambda <- function(x,ph,snpmat=NULL) {
    if(!is.null(snpmat)) {
      if(!snpmat) {
        tt <- table(round(as.numeric(x)),ph)
      }
    }
    if(is.null(snpmat)) {
      tt <- table(round(as.numeric(x)),ph)
      if("3" %in% rownames(tt)) { snpmat <- T } else { snpmat <- F }
    }
    if(snpmat) {
      xx <- round(as.numeric(x))-1; xx[xx==-1] <- NA
      tt <- table(xx,ph)
    }
    #tt <- narm(tt)
    if(checks) {
      # checks will slow down running so can optionally turn them off
      if(length(tt)<1) { warning("all allele counts were empty!"); return(NA) }
      if(!all(rownames(tt) %in% paste(c(0,1,2)))) { 
        warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA) 
      }
      if(!all(colnames(tt) %in% paste(c(0,1,2)))) {
        warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  
      }
    }
    
    if("0" %in% rownames(tt)) { Ctrl0 <- tt["0","0"]; Case0 <- tt["0","1"] } else { Ctrl0 <- Case0 <- 0 }
    if("1" %in% rownames(tt)) { Ctrl1 <- tt["1","0"]; Case1 <- tt["1","1"] } else { Ctrl1 <- Case1 <- 0 }
    if("2" %in% rownames(tt)) { Ctrl2 <- tt["2","0"]; Case2 <- tt["2","1"] } else { Ctrl2 <- Case2 <- 0 }
    Ctrl0[is.na(Ctrl0)] <- 0; Ctrl1[is.na(Ctrl1)] <- 0; Ctrl2[is.na(Ctrl2)] <- 0
    Case0[is.na(Case0)] <- 0; Case1[is.na(Case1)] <- 0; Case2[is.na(Case2)] <- 0
    
    A0 <- Case0+Ctrl0; A1 <- Case1+Ctrl1; A2 <- Case2+Ctrl2
    a0 <- (A0*2)+A1; a2 <- (A2*2)+A1
    if(a0 > a2) { type <- "minor" }
    if(a2 > a0) { type <- "major" }
    if(a0==a2) { type <- "neutral"  }
    if( length(which(c(A0,A1,A2)==0))==2 ) { type <- "monomorph" }
    #cat(type,"\n")
    if(type=="minor") {
      #Case
      r0 <- Case2  ; r1 <- Case1  ; r2 <- Case0  ; R <- Case0+Case1+Case2
      #Control  s0  s1  s2  S
      s0 <- Ctrl2  ; s1 <- Ctrl1  ; s2 <- Ctrl0  ; S <- Ctrl0+Ctrl1+Ctrl2
    } else {
      #Case
      r0 <- Case2  ; r1 <- Case1  ; r2 <- Case0  ; R <- Case0+Case1+Case2
      #Control  s0  s1  s2  S
      s0 <- Ctrl2  ; s1 <- Ctrl1  ; s2 <- Ctrl0  ; S <- Ctrl0+Ctrl1+Ctrl2
    }
    return(c(r0,r1,r2,R,s0,s1,s2,S))  #,n0,n1,n2,N,xx2,yy2,LL,L1000))
    #return(c(LL,L1000))
  }
  ## main code ##
  if(!snp.wise) { cc1000 <- FALSE }
  if(!max(Dim(pheno)) %in% Dim(X)) { warning("Phenotype data different size to dataset X"); return(NA)}
  if(any(is.na(pheno))) {
    bdz <- which(is.na(pheno))
    message("removed ",length(bdz)," samples from X due to phenotype = NA")
    pheno <- pheno[-bdz]
    X <- X[-bdz,]
  }
  if(all(pheno %in% c(1,2))) { pheno <- pheno-1 }
  if(!all(pheno %in% c(0,1))) { warning("Phenotype must be coded as controls,cases=0,1; or =1,2"); return(NA) }
  if(length(Dim(X))!=2) {
    if(length(Dim(X))==1) { return(do.lambda(as.numeric(X),ph=pheno)) } else {
      warning("invalid object for lambda inflation calculation"); return(NA)
    }
  }
  snpmat <- F
  if(is(X)[1] %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { snpmat <- T } else {
    tt.temp <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt.temp)) { snpmat <- T }
  }
  if(snpmat) {
    cnts0 <- get.allele.counts(X[pheno==0,],cc1000=cc1000)
    cnts1 <- get.allele.counts(X[pheno==1,],cc1000=cc1000)
    cnts <- cbind(cnts1,cnts0)
  } else {
    cnts <- apply(X,2,do.lambda,ph=pheno,snpmat=snpmat)
    cnts <- t(cnts)
  }
  #Total
  colnames(cnts) <- c("r0","r1","r2","R","s0","s1","s2","S")
  if(cc1000) { cnts <- round(cnts) }
  #  r0 <- s0 <- r1 <- s1 <- r2 <- s2 <- R <- S <- 0  # initialize these values to zero
  #  tryCatch(detach("cnts"),error=function(e) {NULL}) # in case erroneously attached from earlier
  #  attach(cnts) # should replace all zeros above with values
  n0 <- cnts$r0 + cnts$s0 ; n1 <- cnts$r1 + cnts$s1 ; n2 <- cnts$r2 + cnts$s2 ; N <-cnts$ R + cnts$S
  xx2 <- X_2(cnts$r1,cnts$r2,n1,n2,N,cnts$R)
  yy2 <- Y_2(cnts$r1,cnts$r2,n1,n2,N,cnts$R)
  LL <- yy2/xx2
  L1000 <- lambda_nm(Lnm=LL,n=1000,m=1000,nr=cnts$R,mr=cnts$S)
  #  detach(cnts)
  all.res <- cbind(cnts,xx2,yy2,LL,L1000)
  colnames(all.res)[9:12] <- c("X2","Y2","Lambda","L1000")
  output <- tolower(paste(output[1]))
  output <- gsub("lamda","lambda",output); output <- gsub("lamba","lambda",output)
  if(!output %in% c("lambda","l1000","both")) { output <- "all" }
  if(!snp.wise) {
    # return overall scalar result(s) across all SNPs (median based)
    lam <- medianna(all.res[,"Y2"])/.456
    if(output=="all") { output <- "both" }
    if(output %in% c("both","l1000")) { 
      lam1000 <- lambda_nm(lam,1000,1000,medianna(all.res$R),medianna(all.res$S))
    }
    out <- switch(output,lambda=lam,l1000=lam1000, both=c(Lambda=lam,L1000=lam1000))
  } else {
    # return separate result(s) for each SNP
    out <- switch(output,all=all.res,lambda=all.res[["Lambda"]],
                  l1000=all.res[["L1000"]], both=all.res[,c("Lambda","L1000")])
  }
  return(out)
}


#' Calculate Lambda inflation factors for a SnpMatrixList
#' 
#' This function calculates SNP-wise or overall Lambda and Lambda_1000 statistics for inflation due
#' to population structure. It works on a SnpMatrix object or dataframe coded 0,1,2,NA (autodetects
#'  which).
#' @param X SnpMatrix or data.frame coded 0,1,2,NA containing genotypes
#' @param pheno vector of coding for phenotypes, using (0,1) or (1,2) for controls, cases, should 
#' match nrow(X)
#' @param checks logical, whether to perform some sanity checks that will slightly slown down 
#' performance
#' @param output character vector, containing any of the following "all", "lambda","l1000", "both";
#' if snp.wise is false, the scalar value(s) specified by 'output', or otherwise a matrix
#' of parameters for each SNP, optionally limited to just the lambda column(s) by the value of 
#' 'output'.
#' @param snp.wise logical, if TRUE, return lambda statistics separately for each SNP (rather than
#'  an overall)
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @seealso \code{\link{lambda_nm}}
#' @examples
#' # http://en.wikipedia.org/wiki/Population_stratification
#' # Note that Y2 ~ L*X2, where allele counts are symbolized:
#' # Case     r0  r1  r2  R
#' # Control  s0  s1  s2  S
#' # Total    n0  n1  n2  N
#' require(snpStats) ; data(testdata)
#' pheno <- rep(0,nrow(Autosomes))
#' pheno[subject.data$cc=="case"] <- 1
#' lambdas(Autosomes,pheno)
#' lambdas(Autosomes[,1:5],pheno,snp.wise=TRUE) # list everything snp-wise for the first 5 SNPs
#' lambdas(Autosomes[,1:5],pheno) # just for the first 5 SNPs
lambda.sml <- function(sml, pheno, checks=TRUE, output=c("all","lamba","l1000","both"),snp.wise=FALSE) {
  output <- tolower(paste(output[1]))
  output <- gsub("lamda","lambda",output); output <- gsub("lamba","lambda",output)
  if(!output %in% c("lambda","l1000","both")) { output <- "all" }
  outl <- smlapply(sml,lambdas, pheno=pheno, checks=checks, output="all",snp.wise=TRUE)
  all.res <- do.call("rbind",args=outl)
  #colnames == "r0"     "r1"     "r2"     "R"      "s0"     "s1"     "s2"     "S"      "X2"     "Y2"     "Lambda" "L1000"
  if(!snp.wise) {
   # return overall scalar result(s) across all SNPs (median based)
   lam <- medianna(all.res[,"Y2"])/.456
   if(output=="all") { output <- "both" }
   if(output %in% c("both","l1000")) { 
     lam1000 <- lambda_nm(lam,1000,1000,medianna(all.res$R),medianna(all.res$S))
   }
   out <- switch(output,lambda=lam,l1000=lam1000, both=c(Lambda=lam,L1000=lam1000))
  } else {
   # return separate result(s) for each SNP
   out <- switch(output,all=all.res,lambda=all.res[["Lambda"]],
                 l1000=all.res[["L1000"]], both=all.res[,c("Lambda","L1000")])
  }
  return(out)
}

#(load("~/barrett/OUTPUT/impute-22-17000001-18000001SnpMatrix.RData"))

# si <- sample.info[rownames(X),]


# POINTY #

#' Convert a snpStats SnpMatrix object to a dataframe
#' 
#' Converts a snpStats::SnpMatrix object to a dataframe where coding becomes 0,1,2,NA,
#' which represents genotypes as the number of copies of the reference allele.
#' @param SnpMat a snpStats::SnpMatrix object
#' @param matrix logical, whether to convert to a normal matrix (TRUE), or a data.frame (default, FALSE)
#' @export
#' @return data.frame object with genotype coding where 0,1,2 are the number of copies
#' of the reference allele (by default the letter latest in the alphabet), and NA is
#' missing data.
#' @seealso \code{\link{df.to.SnpMatrix}}
#' @examples
#' require(snpStats)
#' samp <- matrix(sample(c(0:3),18,replace=TRUE),
#'   ncol=3,dimnames=list(paste0("ID0",1:6),c("rs123","rs232","rs433")))
#' samp # note that 0's for SnpMatrix objects are missing values, and other
#' # values are the number of copies of the reference allele plus 1. The reference
#' # allele is by default the letter that is latest in the alphabet, e.g for C/G,
#' # G would be the reference, for T/A, T would be the reference
#' test.mat <- new("SnpMatrix", samp)
#' test.mat # preview does not show data, so use as.data.frame(test.mat)
#' is(as.data.frame(test.mat)[[1]]) # note that using 'as.data.frame' the type is 'raw'
#' SnpMatrix.to.data.frame(test.mat) 
#' # ^ now missing cells are NA, and genotypes are coded as number of copies of reference
SnpMatrix.to.data.frame <- function(SnpMat,matrix=FALSE) {
  #if(is(SnpMat)[1]=="snp.matrix") { SnpMat <- as(SnpMat,"SnpMatrix") }
  if(is(SnpMat)[1]=="aSnpMatrix") { SnpMat <- SnpMat@.Data }
  if(is(SnpMat)[1]!="SnpMatrix") { stop("SnpMat must be a SnpMatrix/aSnpMatrix object") }
  cov.data <- as.data.frame(SnpMat)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
  if(matrix) { cov.data <- as.matrix(cov.data) }
  return(cov.data)
}


#' Convert a data.frame to a snpStats SnpMatrix object
#' 
#' Converts a dataframe to a snpStats::SnpMatrix object where the object contains
#' genotypes coded as number of copies of the reference allele: 0,1,2, and missing=NA.
#' This is an alternative to using new("SnpMatrix",data.frame()). Using 'new' the required 
#' format for the 'data.frame' argument is not as intuitive, as NA's are not allowed, and
#'  2 copies of the reference allele must be coded as 3, and 1 copy as 2, 0 copies as 1.
#' Note that this function will also accept data.frames/matrices coded in that way, and
#' will detect the coding automatically.
#' @param X a data.frame or matrix object containing allele codes 0,1,2 with missing=NA
#' @export
#' @seealso \code{\link{SnpMatrix.to.data.frame}}
#' @return SnpMatrix object
#' @examples
#' require(snpStats)
#' test.frame <- matrix(sample(c(0:2),18,replace=TRUE),
#'   ncol=3,dimnames=list(paste0("ID0",1:6),c("rs123","rs232","rs433")))
#' test.frame[2,2] <- NA # set one genotype to missing
#' test.frame 
#' test.mat <- df.to.SnpMatrix(test.frame) 
#' test.mat
#' snp.mat <- new("SnpMatrix",test.frame) # note this does not handle the NAs/0s nicely
#' all.equal(print(as.data.frame(snp.mat)),print(as.data.frame(test.mat))) # shows offset by 1
#' test.mat
#' row.summary(test.mat) # call rate analysis by sample using snpStats function
#' col.summary(test.mat) # analysis by SNP
df.to.SnpMatrix <- function(X){
  if(is.data.frame(X)) {
    if(any(sapply(lapply(X,is),"[",1) %in% c("character","factor"))) {
      for(cc in 1:ncol(X)) {
        X[[cc]] <- as.numeric(X[[cc]])
      }
    }
  } else { 
    if(!is.matrix(X)) { warning("X should be a matrix or data.frame, ",
                                "conversion is likely to fail if X is not sufficiently matrix-like")}
  }
  mxx <- maxna(X)
  if(mxx>3) { warning("Dataframe does not appear to contain allele codes") }
  X <- round(X)
  if(mxx==3) { X <- X-1 ; X[X<0] <- NA }
  NN <- as.matrix(X)
  #NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}



#' Determine major or minor allele status for a set of SNPs
#' 
#' For a snpStats object or data.frame containing values 0,1,2, NA representing genotypes AA, AB, 
#' BB and no-call. Determines whether the reference allele is the major or minor allele 
#' Where the homozygous genotype coded as highest value = reference, e.g, if AA=0, AB=1, BB=2, 
#' then B is considered the reference here, and by the snpStats package. Combines this with
#' frequencies of the alleles to evaluate whether 'BB' is major or minor. Note that default
#' behaviour for a SnpMatrix is to code alleles alphabetically, so usually the reference allele
#' is the letter later in the alphabet, e.g, it is never an 'A' allele.
#' @param X a SnpMatrix object (see snpStats), or a data.frame coded by reference allele copies,
#' 0,1,2, with missing as NA
#' @param checks logical, whether to perform additional checks for valid values during conversion,
#' setting FALSE will give a slight increase in speed, but is not recommended.
#' @param pheno integer vector, length=nrow(X), must be coded as controls=0 and cases=1, or alternatively
#' using controls=1, cases=2 will be automatically detected and recoded to 0,1.
#' @param tag.mono logical, whether to append a prefix of 'mono' to the major/minor factor code
#' for monomorphic SNPs (minor allele frequency = zero)
#' @param tag.neutral logical, whether to call SNPs with exactly equal A and B allele frequencies
#' 'neutral'; if this option is false, the default is to call neutral SNPs 'minor'.
#' @return returns a factor vector of the same length as the number of SNPs (columns) in the
#' SnpMatrix object/data.frame, indicating for each SNP whether the reference 'B' allele is the 
#' major or minor allele. If the SnpMatrix was creating using a snpStats import function the 
#' reference allele should be the nucleotide letter that is latest in the alphabet. So, for 
#' instance if a SNP is either T/A or A/T, then the reference (B) allele will be 'T'. This 
#' function indicates whether the 'B' allele is the major or minor allele (the major  allele
#' has the greatest frequency). This function can also code 'neutral' if both alleles have 
#' equal frequency when 'tag.neutral' is TRUE, and can add the prefix 'mono.' when 
#' 'tag.mono' is TRUE and one allele has 100% frequency, i.e, is monomorphic. When pheno is given
#' will calculate major/minor alleles with respect to the control MAF only, ignoring cases, where
#' controls as coded as zero.
#' @seealso \code{\link{caseway}}
#' @export
#' @examples
#' require(snpStats)
#' cn <- c("rs123","rs232","rs433","rs234")
#' rn <- paste0("ID0",1:6)
#' samp <- rbind(c(0,2,0,2),c(0,0,0,2),c(0,1,1,2),c(0,0,2,NA),c(0,2,1,1),c(0,2,2,0))
#' dimnames(samp) <- list(rn,cn)
#' test.mat <- df.to.SnpMatrix(samp)
#' majmin(test.mat)
#' col.summary(test.mat) # show call rates, MAF, allele frequencies
#' samp
#' majmin(samp) # also works on a matrix/data.frame with same structure as a Snp.Matrix
#' majmin(test.mat,tag.neutral=TRUE) # rs433 has equal A & B frequencies, and can be tagged neutral
#' majmin(test.mat,tag.mono=TRUE) # rs123 has zero B allele frequency, and can be tagged 'mono'
#' samp[2,2] <- 99 # insert invalid value into matrix
#' majmin(samp,TRUE,TRUE,TRUE) # warning for invalid value
#' majmin(samp,checks=FALSE) # invalid value is converted to NA without warning
majmin <- function(X,pheno=NULL,checks=TRUE,tag.mono=FALSE,tag.neutral=FALSE) {
  ## workhorse internal function ##
  #print(is(X)[1])
  allow.missing.pheno <- TRUE
  if(!is.null(pheno)) {
    if(!max(Dim(pheno)) %in% Dim(X)) { warning("Phenotype data different size to dataset X, returning NA"); return(NA)}
    if(any(is.na(pheno))) {
      if(allow.missing.pheno) {
        X <- X[-which(is.na(pheno)),]
        pheno <- pheno[-which(is.na(pheno))]
        warning("There were NA phenotype values, which will be ignored for MAF calculations") 
      } else {
        warning("There were NA phenotype values, returning NA"); return(NA) 
      }  
    }
    if(all(pheno %in% c(1,2))) { pheno <- pheno-1 }
    if(!all(pheno %in% c(0,1))) { warning("Phenotype must be coded as controls,cases=0,1; or =1,2, returning NA"); return(NA) }
    if(any(pheno==1)) {
      X <- X[which(pheno==0),] # remove cases from matrix
    }
  }
  do.mm <- function(x,snpmat=NULL) { 
    if(!is.null(snpmat)) { 
      if(!snpmat) { 
        tt <- table(round(as.numeric(x))) 
      }
    }
    if(is.null(snpmat)) {  
      tt <- table(round(as.numeric(x)))
      if("3" %in% names(tt)) { snpmat <- T } else { snpmat <- F }
    }
    if(snpmat) {
      xx <- round(as.numeric(x))-1; xx[xx==-1] <- NA
      tt <- table(xx)
    }
    if(checks) {
      # checks will slow down running so can optionally turn them off
      if(length(tt)<1) { warning("all allele counts were empty!"); return(NA) }
      if(!all(names(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
    }
    type <- "unknown"
    if("0" %in% names(tt)) { A0 <- tt[["0"]] } else { A0 <- 0 }
    if("1" %in% names(tt)) { A1 <- tt[["1"]] } else { A1 <- 0 }
    if("2" %in% names(tt)) { A2 <- tt[["2"]] } else { A2 <- 0 } 
    a0 <- (A0*2)+A1; a2 <- (A2*2)+A1
    if(a0 > a2) { type <- "minor" }
    if(a2 > a0) { type <- "major" }
    if(a0==a2) { if(tag.neutral) { type <- "neutral"  } else { type <- "minor" } }
    if( length(which(c(A0,A1,A2)==0))==2 ) { if(tag.mono) { type <- paste("mono",type,sep=".") } }
    return(type)
  }
  ## main code ##
  if(length(Dim(X))!=2) { 
    if(length(Dim(X))==1) { return(do.mm(as.numeric(X))) } else {
      warning("invalid object for major/minor allele testing"); return(NA)
    }
  }
  snpmat <- snp.mat <- F
  if(is(X)[1] %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { snpmat <- snp.mat <- T} else {
    tt1 <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt1)) { snpmat <- T }
  }
  if(snp.mat) {  
    ii <-  col.summary(X)
    #print(ii)
    all.typz <- c("minor","major")[as.numeric(round(ii$RAF,3)!=round(ii$MAF,3))+1]
    if(tag.neutral) { all.typz[round(ii$P.AA,6)==round(ii$P.BB,6)] <- "neutral" }
    if(tag.mono) { all.typz[round(ii$MAF,6)==0] <- paste("mono",all.typz[round(ii$MAF,6)==0],sep=".") }
  } else {
    all.typz <- apply(X,2,do.mm,snpmat=snpmat)
  }
  fl <- c("major","minor")
  if(tag.mono) { fl <- c(fl,"mono.major","mono.minor") }
  if(tag.neutral) { fl <- c(fl,"neutral") }
  return(factor(all.typz,levels=fl))
}



#' Find the direction of GWAS effects between cases and controls
#' 
#' After conducting an case-control association analysis on a SnpMatrix, e.g, 
#' GWAS using for example snp.rhs.tests from the snpStats package, it is not
#' always trivial to determine the direction of the effects with respect
#' to the reference allele. This function calculates which way around the
#' data is for case vs control (pheno) and will indicate with respect to 
#' cases whether they have more copies of the reference allele, or less, or also 
#' can highlight whether the heterozygous is the affected genotype. This function
#' works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which). Note that
#' using a SnpMatrix, with het.effects=FALSE can be much faster (5-100x) than 
#' using a data.frame and/or setting het.effects=TRUE.
#' @param X a SnpMatrix or data.frame containing genotypes for a set of samples
#' @param pheno integer, must be coded as controls=0 and cases=1, or alternatively
#' using controls=1, cases=2 will be automatically detected and recoded to 0,1.
#' @param checks logical, whether to perform additional checks for valid values 
#' during conversion, setting FALSE will give a slight increase in speed, but is 
#' not recommended.
#' @param long logical, whether to use longer more explicit text in the resulting
#' categories describing effect directions, or to use abbreviated codes. For instance,
#' if long==TRUE, then SNPs where cases have more of the reference allele will produce
#' a resulting factor of "cases have more 2, less 0", whereas if long==FALSE, the result
#' would be "CasesRef+".
#' @aliases allow.missing.pheno logical, whether to allow missing values (NA) in the pheno
#' parameter. If FALSE and there are missing values, the function will simply return NA for
#' the whole dataset, with a warning.
#' @param het.effects logical, whether to allow for the possibility that the risk
#' effect is seen for the heterozygous genotype rather than either homozygous allele,
#' if het.effects is TRUE, this will add to two additional possible categories to
#' the output vector, for and increase or decrease in the het allele for cases. If
#' het.effects is FALSE, then the borderline result will default to one side,
#' depending on how the allele is coded in the SnpMatrix object. Note that such
#' instances are unlikely to matter as the most commonly used 1df tests should 
#' be non-significant. Whereas if a 2df test is used, het.effects should be set TRUE.
#' @export
#' @return a factor with the same length as the number of SNPs in X. If het.effects
#' is FALSE and long is FALSE, this will take values 'CasesRef+' indicating cases
#' had a higher frequency of the reference allele than controls, or 'CasesRef-', 
#' indicating that cases had a lower frequency of the reference allele than controls.
#' If 'het.effects'=TRUE, then two more categories are added, 'CasesHet+' and 'CasesHet-',
#' indicating that cases had more of the heterozygous allele, or less, respectively, 
#' compared to controls. Or if 'long' is TRUE, then these four categories are respectively
#' changed to: 'cases have more 2, less 0', 'cases have more 0, less 2', and for het,
#' 'cases have more 1, less 0,2', 'cases have less 1, more 0,2'.
#' Note that these categories are based on absolute counts, and do not necessarily
#' reflect statistical differences in frequency, these should just be used to interpret
#' the direction of your findings, not to perform the analysis.
#' If there are any errors in coding a single NA value will be returned for the whole run,
#' and a warning given as to the reason.
#' @seealso \code{\link{majmin}}
#' @examples
#' require(snpStats)
#' cn <- c("rs123","rs232","rs433","rs234","rs222")
#' rn <- paste0("ID0",1:10)
#' samp <- rbind(c(0,2,0,2,0),c(0,0,0,2,0),c(0,1,1,2,0),c(0,0,2,NA,1),c(0,2,1,1,0),
#'               c(0,2,2,0,2),c(0,2,0,2,2),c(0,1,NA,2,1),c(0,1,0,0,2),c(0,2,2,2,2))
#' dimnames(samp) <- list(rn,cn)
#' test.mat <- df.to.SnpMatrix(samp)
#' phenotype <- c(rep(0,5),rep(1,5))
#' cbind(samp,phenotype) # show example data as a matrix with phenotype
#' caseway(test.mat,phenotype,long=TRUE) # long version
#' caseway(samp,phenotype,het.effects=TRUE) # also works using a data.frame
#' # example conducting association analysis, then adding interpretation of results
#' result <- data.frame(p.value=p.value(snp.rhs.tests(formula=phenotype~1,snp.data=test.mat)))
#' result[["direction"]] <- caseway(test.mat,phenotype)
#' result[["referenceIs"]] <- majmin(test.mat)
#' result # note that only rs222 is significant
caseway <- function(X, pheno, checks=TRUE, long=FALSE, het.effects=FALSE, allow.missing.pheno=TRUE) {
  # coding of output based on long=T/F
  if(long) { r1 <- "cases have more 1, less 0,2" } else { r1 <- "CasesHet+" }
  if(long) { r2 <- "cases have less 1, more 0,2" } else { r2 <- "CasesHet-" }
  if(long) { r3 <- "cases have more 0, less 2" } else { r3 <- "CasesRef-" }
  if(long) { r4 <- "cases have more 2, less 0" } else { r4 <- "CasesRef+" }
  ## workhorse internal function ##
  do.cw <- function(x,ph,snpmat=NULL,r1,r2,r3,r4) { 
    if(!is.null(snpmat)) { 
      if(!snpmat) { 
        tt <- table(round(as.numeric(x)),ph) 
      }
    }
    if(is.null(snpmat)) {  
      tt <- table(round(as.numeric(x)),ph)
      if("3" %in% rownames(tt)) { snpmat <- T } else { snpmat <- F }
    }
    if(snpmat) {
      xx <- round(as.numeric(x))-1; xx[xx==-1] <- NA
      tt <- table(xx,ph)
    }
    #tt <- narm(tt)
    if(checks) {
      # checks will slow down running so can optionally turn them off
      if(length(tt)<1) { warning("all allele counts were empty!"); return(NA) }
      if(!all(rownames(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
      if(!all(colnames(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
    }
    if("0" %in% rownames(tt)) { Ctrl0 <- tt["0","0"]; Case0 <- tt["0","1"] } else { Ctrl0 <- Case0 <- 0 }
    if("1" %in% rownames(tt)) { Ctrl1 <- tt["1","0"]; Case1 <- tt["1","1"] } else { Ctrl1 <- Case1 <- 0 }
    if("2" %in% rownames(tt)) { Ctrl2 <- tt["2","0"]; Case2 <- tt["2","1"] } else { Ctrl2 <- Case2 <- 0 } 
    Ctrl0[is.na(Ctrl0)] <- 0; Ctrl1[is.na(Ctrl1)] <- 0; Ctrl2[is.na(Ctrl2)] <- 0
    Case0[is.na(Case0)] <- 0; Case1[is.na(Case1)] <- 0; Case2[is.na(Case2)] <- 0
    ctrl.pc0 <- Ctrl0/sum(Ctrl0,Ctrl1,Ctrl2); case.pc0 <- Case0/sum(Case0,Case1,Case2)
    ctrl.pc1 <- Ctrl1/sum(Ctrl0,Ctrl1,Ctrl2); case.pc1 <- Case1/sum(Case0,Case1,Case2)
    ctrl.pc2 <- Ctrl2/sum(Ctrl0,Ctrl1,Ctrl2); case.pc2 <- Case2/sum(Case0,Case1,Case2)
    #prv(ctrl.pc0,ctrl.pc1,ctrl.pc2,case.pc0,case.pc1,case.pc2)
    if(long) { res <- "unclear results" } else { res <- "???" }
    if((case.pc0 < ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r1  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r2  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r4  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r3  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r4  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <- r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <- r4  }
    return(res)
  }
  ## main code ##
  typ <- is(X)[1]
  if(typ %in% c("aSnpMatrix","aXSnpMatrix")) {
    if(typ=="aSnpMatrix") { X <- as(X,"SnpMatrix") }
    if(typ=="aXSnpMatrix") { X <- as(X,"XSnpMatrix") }
  }
  if(!max(Dim(pheno)) %in% Dim(X)) { warning("Phenotype data different size to dataset X, returning NA"); return(NA)}
  if(any(is.na(pheno))) {
    if(allow.missing.pheno) {
      X <- X[-which(is.na(pheno)),]
      pheno <- pheno[-which(is.na(pheno))]
    } else {
      warning("allow.missing.pheno was set FALSE, but there were NA phenotype values, returning NA"); return(NA) 
    }  
  }
  if(all(pheno %in% c(1,2))) { pheno <- pheno-1 }
  if(!all(pheno %in% c(0,1))) { warning("Phenotype must be coded as controls,cases=0,1; or =1,2, returning NA"); return(NA) }
  if(length(Dim(X))!=2) { 
    if(length(Dim(X))==1) { return(do.cw(as.numeric(X),ph=pheno)) } else {
      warning("invalid object for case/control effect direction testing, returning NA"); return(NA)
    }
  }
  snpmat <- F
  typ <- is(X)[1]
  if(typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { 
    snpmat <- T
    if(!het.effects) {
      SSTS <- snpStats::single.snp.tests(pheno, snp.data=X, score=T)
      direc <- effect.sign(SSTS)
      if(long) { r3 <- "cases have more of the allele coded '0'" } else { r3 <- "CasesRef-" }
      if(long) { r4 <- "cases have more of the allele coded '2'" } else { r4 <- "CasesRef+" }
      all.res <- rep("???",length(direc))
      all.res[direc==1] <- r4
      all.res[direc==-1] <- r3
      return(factor(all.res))
    }
  } else {
    tt.temp <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt.temp)) { snpmat <- T }
  }
  all.res <- apply(X,2,do.cw,ph=pheno,snpmat=snpmat,r1=r1,r2=r2,r3=r3,r4=r4)
  return(factor(all.res))
}



#' Multicore randomised replacement of missing genotypes
#' 
#' snpStats imputation only works if there are correlated SNPs with non-missing values
#' that can be used to interpolate missing SNPs. If any correlated SNPs are missing
#' 'impute.missing' will leave these blank. This function mops up the remainder
#' by randomly inserting values consistent with the minor allele frequency of each SNP.
#' This can be run using multiple cores to speed up progress for large matrices.
#' @param X a SnpMatrix object, presumably with missing genotypes
#' @param verbose logical, whether to report on progress, or suppress extraneous output
#' @param n.cores integer, maximum number of processing cores to use (limited by system)
#' @return returns a SnpMatrix object with all missing values replaced
#' @export
#' @examples
#' require(snpStats)
#' test <- rSnpMatrix(nsamp=1000,call.rate=.95) # create a random SnpMatrix
#' prv(test) 
#' col.summary(test) # shows some missing values
#' test2 <- randomize.missing(test,verbose=TRUE)
#' col.summary(test2) # shows none missing
randomize.missing <- function(X,verbose=FALSE,n.cores=1) {
  miss1 <- function(x) { 
    x <- as.numeric(x)
    TX <- table(c(round(x),0,1,2,3))-c(1,1,1,1) # to force zero counts to be in the table
    naz <- which(x==0)
    if(length(naz)>0 & length(TX)>0) {
      x[naz] <- sample(as.numeric(names(TX)),size=length(naz),
                       replace=T,prob=as.numeric(TX))
    }
    return(as.raw(x))
  }
  miss2 <- function(sel,aa,ab,bb) { 
    x <- X@.Data[,sel]
    naz <- which(x==0)
    p.vec <- c(aa,ab,bb)
    if(any(is.na(p.vec))) {
      if(all(is.na(p.vec))) {
        x[naz] <- as.raw(rsnp(length(naz),cr=1))
        return(x)
      } else {
        stop(p.vec)
        #Pleft <- 1-sum(p.vec,na.rm=T)  
      }
    }
    x[naz] <- sample(as.raw(1:3),size=length(naz),replace=T,prob=p.vec)
    return(x)
  }
  # randomly generate replacements for missing values using current distribution for each column of X
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { stop("X must be a SnpMatrix object or similar") }
  cs <- col.summary(X)
  FF <- nrow(X)
  nmiss <- FF - cs$Calls
  select <- nmiss>0
  if(length(which(select))<1) { return(X) }
  if(verbose) { cat(sum(nmiss),"missing values being replaced with random alleles .") }
  if(length(which(select))==1) { X[,select] <- miss2(select,cs$P.AA[select],cs$P.AB[select],cs$P.BB[select]); return(X) }
  if(n.cores>1) {
    do.chunk <- function(SEL) {
      ii <- mapply(miss2,SEL,cs$P.AA[SEL],cs$P.AB[SEL],cs$P.BB[SEL])
      cat(".")
      return(ii)
    }
    posz <- which(select)
    Ls <- length(posz)
    split.to <- min(c(round(Ls/100),n.cores),na.rm=T)
    #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
    stepz <- round(seq(from=1,to=Ls+1,length.out=round((split.to+1))))
    if((tail(stepz,1)) != Ls+1) { stepz <- c(stepz,Ls+1) }
    split.to <- length(stepz)-1
    outlist <- sel.range <- vector("list",split.to)
    for (cc in 1:split.to) {
      c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
      sel.range[[cc]] <- posz[c(c1:c2)]
    }
    prv(split.to,stepz,select,Ls,sel.range)
    outlist <- mclapply(sel.range,do.chunk,mc.cores=n.cores)
    for (cc in 1:split.to) {
      X@.Data[,sel.range[[cc]]] <- outlist[[cc]]
    }
  } else {
    ii <- mapply(miss2,which(select),cs$P.AA[select],cs$P.AB[select],cs$P.BB[select])
    X@.Data[,select] <- ii
  }
  cat("done\n")
  return(X)
}



#' Replace missing values in a SnpMatrix object with imputed values
#'
#' This function is a wrapper for the snpStats snp.imputation() and impute.snps() functions.
#' It allows a full imputation with one simple command, and facilitates stratified imputation for
#' subsets of samples, and the parameter 'by' allows subsetting by SNPs, which speeds up the
#' imputation for large sets by running it in smaller batches (as large ones are really slow).
#' The standard use of the snpStats imputation functions will still leave some NA's behind,
#' whereas the option 'random' ensures each missing value is replaced, even if just by an
#' allele chosen at random (using existing frequencies).
# (chris)
#' @param X SnpMatrix object with missing data that you wish to impute
#' @param strata factor, imputation can be done independently for different sample subsets,
#' use a vector here of the same length as the number of samples (nrow(X)) to code this
#' @param by integer, a parameter that is passed to impute.missing() to determine
#' what sized subsets to impute by (smaller subsets makes imputation faster)
#' @param random logical, when imputation is performed using the core function, impute.snps(),
#' usually there are some SNPs that remain missing, due to missingness in the most correlated
#' SNPs used to impute them. Setting this parameter TRUE, will replace this remainder with
#' random values (using the allele frequencies for the non missing data). Setting to FALSE
#' will leave these values as missing in which case you cannot rely on this function returning
#' a fully complete matrix.
#' @param data.frame logical, if TRUE then return the result as a data.frame instead of a
#' SnpMatrix object
#' @param verbose logical, if TRUE, more information about what is being done in the 
#' imputation will be provided. Additionally, if TRUE and the dataset being imputed is 
#' large enough to take a long time, then a progress bar will be displayed using
#' 'loop.tracker()'. If FALSE, no progress bar will be used, regardless of the size of dataset.
#' @param ... further arguments to 'impute.snps()' from snpStats, which is the core function
#' that performs the main imputation.
#' @return Returns a SnpMatrix with no missing values (as any originally present are imputed
#' or randomized), or if random=FALSE, it may contain some missing values. If data.frame=TRUE
#' then the returned object will be a data.frame rather than a SnpMatrix.
#' @export
#' @author Chris Wallace and Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' require(snpStats)
#' test.mat <- rSnpMatrix(nsnp=5,nsamp=10)
#' print(SnpMatrix.to.data.frame(test.mat))
#' out.mat <- impute.missing(test.mat)
#' print(SnpMatrix.to.data.frame(out.mat))
#' test.mat2 <- rSnpMatrix(nsnp=200,nsamp=100,call.rate=.98)
#' head(col.summary(test.mat2))
#' out.mat2 <- impute.missing(test.mat2,by=50)
#' cs <- col.summary(out.mat2)
#' head(cs)
#' # 'random' is true, so there should be no remaining missing values
#' paste(out.of(length(which(cs$Call.rate==1)),nrow(cs)),"SNPs imputed to 0% missing\n")
#' out.mat3 <- impute.missing(test.mat2,random=FALSE)
#' cs <- col.summary(out.mat3)
#' # random was FALSE, so some of these should still have missing values:
#' paste(out.of(length(which(cs$Call.rate==1)),nrow(cs)),"SNPs imputed to 0% missing\n")
#' # for instance, these
#' print(head(cs[which(cs$Call.rate<1),]))
#' # complicated example using both sample and SNP subsets for imputation
#' out.mat4 <- impute.missing(test.mat2, by=100, 
#'                       strata = c(rep(1,25),rep(2,50),rep(3,25)), random=FALSE)
impute.missing <- function (X,  strata = NULL, by=NULL, random=TRUE, data.frame = FALSE,
                            verbose=TRUE,...) {
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","aSnpMatrix")) { stop("X must be a SnpMatrix/aSnpMatrix object") }
  if(Dim(X)[1] < 2 | Dim(X)[2] <2) { stop("X must be at least 2 x 2") }
  N <- as(X, "numeric"); 
  done <- FALSE
  mem.max <- 0.02  #GB
  bp <- 1:ncol(X) # parameter needed for snp.imputation() ?
  if(!is.null(by)) {
    if(!is.numeric(by)) { warning("by must be numeric, ignoring invalid parameter"); by <- NULL }
    by <- round(by) # force integer
    if(by<3) { warning("by is too small to be useful, ignoring parameter"); by <- NULL }
    if((by*1.25)>ncol(X)) { by <- NULL } # ignore 'by' if num SNPs is not that large
  }
  if(any(Dim(X)!=Dim(N))) { stop("'as numeric' caused loss of dimensionality for unknown reason") }
  if (!is.null(strata)) {
    if(length(strata)==nrow(X)) {
      ## Partition by sample subsets, then re-send the subsets to this function for processing
      strata <- as.factor(strata)
      if (length(levels(strata)) > 20) 
        stop("too many levels in strata\n")
      for (i in levels(strata)) {
        if(verbose) { cat("\nImputing sample subset '", i, "' (strata)\n",sep="") }
        wh <- which(strata == i)
        N <- as.data.frame(N)
        N[wh, ] <- impute.missing(X[wh, , drop = FALSE], 
                                  data.frame = TRUE, by=by,...)
      }
      done <- TRUE
    } else {
      #print(length(strata));print(nrow(X))
      warning("'strata' should have length==nrow(X), all samples will be imputed together")
    }
  }
  if (is.null(strata) & is.numeric(by)) {
    ## Partition by SNP subsets, then re-send the subsets to this function for processing
    strata <- as.factor(make.split(ncol(X),by=by))
    # added by nick
    for (i in levels(strata)) {
      wh <- which(strata == i)
      if(verbose) { cat(" imputing SNP subset ", minna(wh), "-",maxna(wh),"\n",sep="") }
      anything <- impute.missing(X[,wh , drop = FALSE], data.frame = TRUE, ...)
      N <- as.data.frame(N)
      N[, wh] <- anything
    }
    done <- TRUE
  } else {
    if(!done) {
      # The core imputation, done using the snpStats functions 'snp.imputation' and 'impute.snps' #
      csumm <- col.summary(X)
      use <- csumm[, "Certain.calls"] == 1
      X2 <- X[, use]
      bp <- bp[use]
      imp <- (csumm[, "Call.rate"] < 1 & !is.na(csumm[, "Call.rate"]))[use]
      #if(sumna(imp)==175) { stst <- strata; byby <- by; prv(stst,byby)}
      if(verbose) { cat(sumna(imp), " SNP",if(sumna(imp)!=1) {"s"} else { "" }," in the set had missing values\n",sep="")  }
      ll <- length(which(imp)); dd <- 1 # mx <- max(which(imp),na.rm=T); 
      if(estimate.memory(X2)<mem.max) { verbose <- FALSE }
      for (i in which(imp)) {
        if(verbose) {
          loop.tracker(dd,ll); dd <- dd + 1 # track this loop as count, these pars not used elsewhere
        }
        supres <- capture.output(rule <- snp.imputation(X2[, -i, drop = FALSE], X2[, 
                                                                                   i, drop = FALSE], bp[-i], bp[i]))
        if (is.null(rule@.Data[[1]])) 
          next
        imp <- impute.snps(rules = rule, snps = X2[, rule@.Data[[1]]$snps, drop = FALSE], ...) 
        wh.na <- which(is.na(N[, i]))
        #gtWWW <- X2
        #prv(gtWWW)
        if(length(wh.na)>0) { 
          N[wh.na, colnames(X2)[i]] <- imp[wh.na]
        }
      }
    }
  }
  if(random) { N <- randomize.missing2(N,verbose=verbose) } # replace any remaining missing values
  if(data.frame) {
    return(as.data.frame(N))
  }
  else {
    #print(Dim(N))
    #print(Dim(X@snps))
    N <- round(N)
    N[is.na(N)] <- -1 # replace missing with -1 [will become zero in next line]
    if(is(N)[1]=="list" | is(N)[1]=="data.frame") { N <- as.matrix(N) } # sometimes it becomes a list for some reason???
    out.m <- (new("SnpMatrix", data = (as.raw(N+1)),
                  nrow = nrow(N), ncol = ncol(N), dimnames = dimnames(N)))
    if(typ=="SnpMatrix") { return(out.m) } else {
      X@.Data <- out.m
      return(X)
    }
  }
}






#' Create a SNP matrix with simulated data
#' 
#' A simple function to simulation random SnpMatrix objects for
#' testing purposes. Does not produce data with an 'LD' structure
#' so is not very realistic. 
#' @param nsnp number of SNPs (columns) to simulate
#' @param nsamp number of samples (rows) to simulate
#' @param call.rate numeric, percentage of SNPs called, e.g, .95 leaves 5% missing data,
#' or 1.0 will produce a matrix with no missing data
#' @param A.freq.fun function, a function to randomly generate frequencies of the 'A'
#' allele, this must produce 'n' results 0 < x < 1 with argument n, the default is
#' for uniform generation.
#' @export
#' @return returns a SnpMatrix object with nsnp columns and nsamp rows, with 1-call.rate% 
#' missing data, with allele frequencies generated by 'A.freq.fun'.
#' @examples
#' require(snpStats)
#' newMat <- rSnpMatrix(call.rate=.92) # specify target of 8% missing data
#' print(as.data.frame(newMat)) # to see actual values
rSnpMatrix <- function(nsamp=10, nsnp=5, call.rate=.95, A.freq.fun=runif) {
  dummy.vec <- rep(nsamp,nsnp)
  call.rate <- force.percentage(call.rate)
  exp.fac <- 10
  ld <- FALSE
  warn.mem <- 0.5 # 0.5GB
  if(!ld) {
    mat <- sapply(dummy.vec,rsnp,A.freq.fun=A.freq.fun,cr=call.rate)
  } else {
    # this didn't really work #
    if(estimate.memory(exp.fac*nsamp*nsnp)>warn.mem) { warning("when ld=TRUE, large simulations can take a long time") }
    mat <- sapply(rep(nsamp,nsnp*exp.fac),rsnp,A.freq.fun=A.freq.fun,cr=call.rate)
    mat <- get.top.n(mat,nsnp)    
  }
  rownames(mat) <- rsampid(nsamp); colnames(mat) <- rsnpid(nsnp)
  #sm <- new("SnpMatrix",as.raw(mat),dimnames=list(rsampid(nsamp),rsnpid(nsnp)))
  sm <- new("SnpMatrix", data = as.raw(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  return(sm)
}




################## endSnpmatrix ##########################


#####################
## Other Functions ##
#####################


#' Function to produce a clean table from logistic regression done via GLM
#' 
#' With input as a glm model result object, returns a clean dataframe with 
#' coefficients, p values, confidence intervals and standard errors. Multiple
#' options for which columns to return and digits to display
#'
#' @param glm.result should be type 'glm', the result of a call to glm() where
#' the family parameter specifies logistic regression, i.e, family=binomial("logit")
#' @param intercept logical, whether to display the intercept as the first row of results
#' @param ci logical, whether to include upper and lower confidence limits in the table of results
#' @param se logical, whether to include standard error (SE) in the table of results
#' @param alpha percentage, the type 1 error rate to apply to calculation of confidence
#' limits, e.g, alpha=.01 results in a 99\% confidence interval.
#' @param o.digits, integer, number of digits to use displaying odds-ratios and standard errors
#' @param p.digits, integer, number of digits to round p-values to, not that a high number ensures
#' that very small p-value are truncated to be displayed as zero.
#' @param label, logical, whether to label output matrices as a list, with name corresponding to 
#' the formula for the model
#' @return a clean dataframe with coefficients, p values, confidence intervals and standard errors. 
#' There are multiple options for which columns to return and number of digits to display.
#' @export
#' @examples
#' # do we need to require(stats)?
#' ph <- c(0,0,0,1,1,1) # phenotype
#' X <- rnorm(6) # independent variable
#' test.glm <- glm(ph ~ X,family=binomial('logit'))
#' logistic.summary(test.glm) # very simple example
#' XX <- sim.cor(10,3)
#' X1 <- XX[,1]; X2 <- XX[,2]; X3 <- XX[,3]
#' ph <-c(rep(0,5),rep(1,5))
#' test.glm <- glm(ph ~ X1 + X2 + X3,family=binomial('logit'))
#' logistic.summary(test.glm,intercept=TRUE) # include intercept
#' logistic.summary(test.glm,ci=TRUE,se=FALSE) # show 95% confidence intervals, hide standard error
#' logistic.summary(test.glm,ci=TRUE,alpha=.01) # get 99% confidence intervals
logistic.summary <- function(glm.result,intercept=FALSE,ci=FALSE,se=TRUE,alpha=.05,o.digits=3,p.digits=299,label=FALSE)
{
  if(!"glm" %in% is(glm.result)) { stop("'glm.result' was not a glm result (of type 'glm'") }
  if(family(glm.result)[[1]]!="binomial" | family(glm.result)[[2]]!="logit") { 
    warning("this function is designed to work for glm logistic regression, i.e; ",
            "glm(family=binomial('logit'))", " invalid output is likely")
  }
  co <- summary(glm.result)$coefficients
  alpha <- force.percentage(alpha,default=.05)
  if(rownames(co)[1]!="(Intercept)") { intercept <- TRUE }
  predz <- rownames(co)
  if(!intercept) { predz <- predz[-1] }
  lab <- paste(summary(glm.result)$call)[2]
  #prv(co)
  if(length(Dim(co))<2) {
    warning("glm result was empty of parameters")
    return(NULL)
  }
  if(!intercept & nrow(co)<2) { 
    warning("intercept was not in use, and there were no other parameters in the models") 
    return(NULL) 
  }
  if(!intercept) { fst <- 2 } else { fst <- 1 }
  p <- co[fst:nrow(co),4]; #print(p)
  o.r <- exp(co[fst:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  
  # outlist <- list(round(o.r,o.digits),p)
  # names(outlist) <- c("OR","p-value")
  SE <- co[fst:nrow(co),2]
  if(ci) {
    #prv(co)
    co1 <- co[fst:nrow(co),1]-(p.to.Z(alpha)*SE)
    co2 <- co[fst:nrow(co),1]+(p.to.Z(alpha)*SE)
    if(sum(co1)<=sum(co2)) { co.l <- co1; co.h <- co2 } else {  co.h <- co1; co.l <- co2 }
    co.l <- exp(co.l); co.h <- exp(co.h)
    out.fr <- cbind(round(o.r,o.digits),round(SE,o.digits),round(co.l,o.digits),round(co.h,o.digits),p)
    colnames(out.fr) <- c("OR","SE","OR-low","OR-hi","p-value")
  } else {
    out.fr <- cbind(round(o.r,o.digits),round(SE,o.digits),p)
    colnames(out.fr) <- c("OR","SE","p-value")
  }
  if(is.null(rownames(out.fr)) & nrow(out.fr)==1) {
    rownames(out.fr) <- predz
  }
  if(!se) {
    if(length(Dim(out.fr))==2) {
      out.fr <- out.fr[,-which(colnames(out.fr) %in% "SE")]
    } else {
      if(length(Dim(out.fr))==1) {
        out.fr <- out.fr[-which(colnames(out.fr) %in% "SE")]
      }
    }
  }
  if(label) { out.fr <- list(out.fr); names(out.fr) <- lab }
  return(out.fr)
}







#' Meta-analysis using odds ratio and standard error from 2 datasets
#' 
#' This function calculates meta analysis odds ratios, standard errors and p-values
#' using results from a table containing odds ratio and standard error data for analyses
#' of 2 different datasets (typically logistic regression, but other analyses can be
#' incorporated if an odds-ratio and SE can be derived, for instance one analysis might
#' be a case control logistic regression GWAS and the other a family TDT analysis).
#' @param X A data.frame with column names which should be entered in the parameters:
#' OR1, OR2, SE1, SE2, and optionally N1, N2. 
#' @param OR1 The column name of X containing odds ratios from the first analysis
#' @param OR2 Same as OR1 above but pertaining to the second analysis
#' @param SE1 The column name of X containing standard errors from the first analysis
#' @param SE2 Same as SE1 above but pertaining to the second analysis
#' @param N1 Only required if method="sample.size". Either the column name in X with the 
#' number of samples in the first analysis, of a vector of the same, or if N's is the
#'  same for all rows, a scalar' value can be entered
#' for each.
#' @param N2 Only required if method="sample.size". Same as N1 above but pertaining to analysis 2
#' @param method character, can be either 'beta', 'z.score' or 'sample.size', and upper/lower
#' case does not matter. 'Beta' is the default and will calculate meta-analysis weights using
#' the inverse variance method (based on standard errors), and will calculate the p-values
#' based on the weighted beta coefficients of the two analyses. 'Z.score' also uses inverse variance
#' but calculates p-values based on the weighted Z scores of the two analyses. 'Sample.size' uses
#' the sqrt of the sample sizes to weight the meta analysis and uses Z scores to calculate p values
#' like 'Z.score' does.#' 
#' @return The object returned should have the same number of rows and rownames as the data.frame
#'  X but columns are the meta analysis stastistics, namely:
#'   OR.meta, beta.meta, se.meta, z.meta, p.meta, which will contain the meta
#' analysis odds-ratio, beta-coefficient, standard error, z-score, and p-values respectively
#' for each row of X.
#' @export
#' @examples
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_Fam=c(1.33,0.95),SE_CC=c(0.02,0.12),SE_Fam=c(0.07,0.5))
#' rownames(X) <- c("rs689","rs23444")
#' X
#' meta.me(X)
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_CC2=c(1.33,0.95),
#'  SE_CC=c(0.02,0.12),SE_CC2=c(0.02,0.05),
#'  n1=c(5988,5844),n2=c(1907,1774))
#' # even with roughly the same number of samples the standard error will determine the influence of
#' # each analysis on the overall odds ratio, note here that the second SE for dataset goes
#' # from 0.5 to 0.05 and as a result the estimate of the odds ratio goes from 1.137 to 0.977,
#' # i.e, from very close to OR1, changing to very close to OR2.
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2") 
#' # sample size and z-score methods give similar (but distinct) results
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="sample.size") 
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="z.score")  # N's will be ignored
meta.me <- function(X,OR1="OR_CC",OR2="OR_Fam",SE1="SE_CC",SE2="SE_Fam",Z1=NA,Z2=NA,
                    N1=NA,N2=NA,method=c("beta","z.score","sample.size")) {
  #N1=18856,N2=7638
  validz <- c("beta","z.score","sample.size")
  method <- tolower(method[1])
  if(!method %in% validz) { method <- validz[1]; warning("invalid method entered, using 'beta' method") }
  if(!is(X)[1]=="data.frame") { stop("X must be a data.frame") }
  cnx <- colnames(X)
  if(is.null(rownames(X))) { rownames(X) <- paste(1:nrow(X)) }
  if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { stop("X must contain column names specified by OR1,OR2,SE1,SE2") } 
  ok <- FALSE
  if(method=="sample.size") {
    if(is.numeric(N1) & is.numeric(N2)) { if(length(N1)!=1 | length(N2)!=1) { 
      stop("N1, N2 should either be scalar integers, or column names containing N's") } else {
        N.coln <- FALSE
      } }
    if(is.character(N1) & is.character(N2)) { 
      if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { 
        stop("N1,N2 must contain either be scalar integers or column names with N's") } else {
          N.coln <- TRUE
        }
    }
  }
  OR.CC <- X[,OR1]
  beta.CC  <- log(X[,OR1])
  se.CC <- X[,SE1]
  OR.family <- X[,OR2]
  beta.family  <- log(X[,OR2])
  se.family <- X[,SE2]
  if(!is.na(Z1) & method!="beta") {
    if(Z1 %in% colnames(X)) {
      z.CC <- X[,Z1]
    } else { warnings("Z1 column not found, ignoring") }
  } else {    
    z.CC <- beta.CC/se.CC
  }
  if(!is.na(Z2) & method!="beta") {
    if(Z2 %in% colnames(X)) {
      z.family <- X[,Z2]
    } else { warnings("Z2 column not found, ignoring") }
  } else {
    z.family <- beta.family/se.family
  }
  inv.CC <- 1 / (se.CC^2)
  inv.family <- 1 / (se.family^2)
  var.meta <- 1 / (inv.CC+inv.family)
  weight.CC <- inv.CC * var.meta
  weight.family <- inv.family * var.meta
  se.meta <- round(sqrt(var.meta), digits=3)
  if(method=="sample.size") {
    if(N.coln) {
      famN <- X[,N2]
      ccN <- X[,N1]
    } else {
      famN <- N2 # 3819*2  #3509*2   #  3819*2   #  10796
      ccN <- N1 # 6683+12173 # including CBR, or for UVA analyses use instead: 9416+6670
    }  
    WeightFam = sqrt(famN)/(sqrt(famN)+sqrt(ccN))
    WeightCC <- 1-WeightFam
    # beta calculated the same way for sample.size method using sample size weights
    beta.meta <- round((WeightCC * beta.CC) + (WeightFam * beta.family),digits=3) # beta based
    z.meta <- round((WeightCC * z.CC) + (WeightFam * z.family),digits=6) # n-based, z-based
  } else {
    # beta calculated the same way for z.score method and beta method using inverse variance
    beta.meta <- round((weight.CC * beta.CC) + (weight.family * beta.family),digits=3) # beta based
    if(method=="z.score") {
      z.meta <- round((weight.CC * z.CC) + (weight.family * z.family),digits=6) # z-based
    } else {
      # default (beta) method
      z.meta <- beta.meta/se.meta
    }
  }  
  OR.meta <- exp(beta.meta)
  p.meta <- 2*pnorm(-abs(z.meta))
  out <- (cbind(OR.meta,beta.meta,se.meta,z.meta,p.meta)) 
  colnames(out) <- c("OR.meta","beta.meta","se.meta","z.meta","p.meta") 
  rownames(out) <- rownames(X)
  return(out)
}


################## end other ##########################



######################
## Simple Functions ##
######################

#' Normalize Lambda inflation factors to specific case-control count
#' 
#' Lambda inflation statistics are influenced by the size of the generating datasets. To facilitate
#' comparison to other studies, this function converts a given lambda from nr cases and mr controls,
#' to n cases and m controls, where m=n=1000 is the most common normalization. All values other than
#' 'Lnm' are forced within this function to be positive integers.
#' @param Lnm numeric, a raw Lambda inflation statistic, generated from nr cases and mr controls
#' @param n integer, desired number of 'virtual' cases to normalise to
#' @param m integer, desired number of 'virtual' controls to normalise to
#' @param nr integer, original number of cases that Lnm was derived from
#' @param mr integer, original number of controls that Lnm was derived from
#' @return A normalized Lambda coefficient
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @seealso \code{\link{lambdas}}
#' @examples
#' require(snpStats) ; data(testdata) # already 'depends'?
#' pheno <- rep(0,nrow(Autosomes))
#' pheno[subject.data$cc=="case"] <- 1
#' raw.L <- lambdas(Autosomes,pheno,output="lambda")
#' L1000 <- lambda_nm(raw.L,1000,1000,200,200)
#' raw.L; L1000
#' lambda_nm(1.56,1000,1000,6500,9300) # lambda1000<lambda when n's>1000,1000
lambda_nm <- function(Lnm,n=1000,m=1000,nr,mr) { 
  if(!is.numeric(Lnm)) { stop("Lnm must be numeric") }
  if(!is.numeric(n)) { stop("n must be numeric") } else { n <- abs(round(n)) }
  if(!is.numeric(m)) { stop("m must be numeric") } else { m <- abs(round(m)) }
  if(!is.numeric(nr)) { stop("nr must be numeric") } else { nr <- abs(round(nr)) }
  if(!is.numeric(mr)) { stop("mr must be numeric") } else { mr <- abs(round(mr)) }
  return(1 + ((Lnm-1)*(((1/nr)+(1/mr))/((1/n)+(1/m)))) )
}


#' Calculate approximate Bayes factors from p values and MAF
#'
#' This is a function to calculate approximate Bayes factors from p
#' values and MAF - for reference see Wakefield, J (2009) Bayes
#' factors for genome-wide association studies: comparison with
#' p-values. 
#' @param p p value
#' @param maf minor allele frequency
#' @param n0 number of controls
#' @param n1 number of cases
#' @param scale0 by default, =n0
#' @param scale1 by default, =n1
#' @references Genetic Epidemiology 33: 79-86.
#' @return the appoximate bayes factor calculated based on
#' the p-value, minor allele frequency and case/control sample
#' sizes
#' @export
#' @author Chris Wallace
#' @examples
#' abf(.0001,.2,9500,6670)
abf <- function(p,maf, n0=1000, n1=1000, scale0=n0, scale1=n1) { 
  # evidence for null - ie low ABF => support for alternative
  z <- qnorm(p/2, lower.tail=FALSE)
  x0 <- 0; x1 <- 1; x2 <- 2 # multiplicative model
  d2 <- (1-maf)^2 * x0^2 + 2*maf*(1-maf)*x1 + maf^2 * x2^2
  d1 <- (1-maf)^2 * x0 + 2*maf*(1-maf)*x1 + maf^2 * x2
  V <- (n0 + n1) / ( n0 * n1 * (d2-d1) )
  ## scale
  scale <- ((n0 + n1)/(scale0 + scale1)) * (scale0/n0) * (scale1/n1)
  V <- V * scale
  W <- ( log(1.5)/qnorm( 0.99, lower.tail=FALSE) )^2
  VW <- V+W
  abf <- 2 * log(sqrt(VW/V) * exp( - z^2 * W / (2 * VW) ))
  return(abf)
}



#' Import a ped file to pedData format used by snpStats
#'
#' PLINK ped files (family information files) are not in the
#' same format as ped files required for use with snpStats, for
#' instance for the tdt.snp() function to conduct a transmission
#' disequilibrium test. This function will import a PLINK style
#' ped file and return a 'pedData' object in the correct form
#' for snpStats and other rpackages. The plink file format is:
#' column 1: family id, column 2: individual id: column 3:
#' father's ID or 0, column 4: mother's ID or 0, column 5: sex of subject,
#' column 7: phenotype of subject.
#' @param file character, the name of a valid PLINK 'ped'/'fam' file
#' @param correct.codes logical, if TRUE, then where coding seems
#' inconsistent, e.g, mother and father both the same, this will
#' try to fix this automatically
#' @param silent logical, when false, output will show what
#' operations and transformations are being done, when TRUE,
#' the function produces no interim text to the console
#' @export
#' @return return a pedData object (a data.frame)
#' @examples
#' \donttest{
#' # fn <- "myfile.ped" # insert name of your own ped file
#' # myPed <- read.pedData(fn); prv(myPed)
#' }
read.pedData <- function(file,correct.codes=TRUE,silent=FALSE) {
  rr <- read.ped.file(file,keepsix = TRUE)
  want <- c("familyid","individual","father","mother","sex","affected")
  have <- colnames(rr)
  if(!silent) { cat(paste("mapping column",have,"==>",want,"\n"),sep="") }
  if(ncol(rr)>6) { warning("expected 6 columns in pedfile, unexpected behaviour could result") }
  if(ncol(rr)<6) { stop("need at least 6 columns in pedfile to map to headings ",paste(want,collapse="")) }
  colnames(rr)[1:length(want)] <- want
  rr <- rr[order(rr$familyid),]
  rr[["member"]] <- unlist(tapply(rep(1,nrow(rr)),factor(rr$familyid),cumsum))
  rr <- rr[,c(2,1,7,3:6)]
  rr <- shift.rownames(rr,T)
  rr[["father"]] <- rr$member[match(rr$father,rownames(rr))]
  rr[["mother"]] <- rr$member[match(rr$mother,rownames(rr))]
  if(correct.codes) {
    badz <- which(rr$father==1 & rr$mother==1)
    rr[badz,"mother"] <- 2
    badz <- which(rr$father==2 & rr$mother==2)
    rr[badz,"father"] <- 1
  }
  return(rr)
}

# prefer to have the following 2 functions in humarray 
# use aSnpMatrix or SnpMatrix plus data.frame of allele.codes
# allele codes must be in the form alt, ref = actually unnecessary unless group is skewed so much MAFs are flipped- update later for this
estimate.HLA <- function(chr6,allele.codes=NULL,snp.names=c("rs3104413", "rs2854275", "rs9273363"),assume.mafs=FALSE) {
   rs3vals <- c("CC","GC","GG")  # rs3104413  C > G (~26% MAF)
   rs2vals <- c("CC","AC","AA") # rs2854275  C > A (~10% MAF)
   rs9vals <- rev(c("AA","AC","CC")) # rs9273363  C > A (~33% MAF) according to immunobase, but paper suggests A is major
   sw <- logical(3) ; sw <- rep(F,3)
   #########################
   chr6 <- chr6[,snp.names]
	if(is(chr6)[1]=="aSnpMatrix") { 
	   allele.codes <- snps(chr6)[,alleles(chr)] 
	} else { 
		if(!is.data.frame(allele.codes)) {
			if(! assume.mafs) { 
				stop("if not using an aSnpMatrix then allele.codes must be a n x 2 data.frame") 
			}
		} else {
			if(is(chr6)[1]!="SnpMatrix") { stop("chr6 must be a SnpMatrix or aSnpMatrix") }
			if(ncol(chr6)!=nrow(allele.codes)) { stop("allele.codes must have rows corresponding to the number of columns (SNPs) in chr6")}
		}
   }
      	   cs <- col.summary(chr6[,snp.names])
   	 #  print(cs)
   if(!assume.mafs) {
   		#if(allele.codes[1,2]=="C") { sw[1] <- T }
   		if(cs$MAF[1]!=cs$RAF[1]) { sw[1] <- T }
   		if(allele.codes[2,2] %in% c("C","G")) { sw[2] <- T }
  		if(allele.codes[3,2] %in% c("C","G")) { sw[3] <- T }
   } else {

   	   for (cc in 1:3) {	if(cs$MAF[cc]!=cs$RAF[cc]) { sw[cc] <- T } }
   }
   if(all(snp.names %in% colnames(chr6))) {
   	   df <- SnpMatrix.to.data.frame(chr6[,snp.names])
   	 #  prv(df)
   	   # make sure the ref is the minor allele
   	   for (cc in 1:3) {
   	   	  if(sw[cc]) { 
   	   	  	 df[,cc] <- 2-df[,cc] 
   	   	  }
   	   }
   	  # prv(df)
   	   rs3104413 <- rs3vals[(df[,1]+1)]; rs3104413[is.na(rs3104413)] <- ""
   	   rs2854275 <- rs2vals[(df[,2]+1)]; rs2854275[is.na(rs2854275)] <- ""
   	   rs9273363 <- rs9vals[(df[,3]+1)]; rs9273363[is.na(rs9273363)] <- ""
   } else {
   	  stop("chr6 must contain the SNPs: ",paste(snp.names,collapse=",")," to calculate HLA types")
   }
   ll <- length(rs3104413)
   results <- data.frame(DR = rep("",ll), DQ = rep("",ll), risk = rep("",ll),stringsAsFactors=FALSE)
   rownames(results) <- rownames(chr6) ; if(all(rownames(chr6)!=rownames(df))) { warning("has row order changed during SnpMatrix conversion?") }
   #for (cc in 1:3) { results[[cc]] <- as.character(results[[cc]])}
   for (dd in 1:ll) {
   	   results[dd,] <- paste(apply.3snp.rules(rs3104413[dd], rs2854275[dd], rs9273363[dd]))
   }
   return(results)
}

apply.3snp.rules <- function(rs3104413 = "GG", rs2854275 = "AC", 
	rs9273363 = "AA") {
	risk <- "normal"
	DQ <- ""; DR <- ""
	if (rs3104413 == "GG") {
		DR <- "DR4/4"
		if (rs9273363 %in% c("AA", "AC")) {
			DQ <- "DR4-DQ8"
			risk <- "high"
		} else {
			DQ <- "x"
			risk = "normal"
		}
	} else {
		if (rs3104413 == "GC") {
			if (rs2854275 == "AC") {
				DR <- "DR3/4"
				if (rs9273363 == "AA") {
					DQ <- "DR3/4-DQ8"
					risk <- "very high"
				} else {
					DQ <- "DR3/4-DQ2" #\DRB1*03:01-DQA1*05:01-DQB1*02:01\"; risk <- \"high\""
					risk <- "normal"
				}
			} else {
				if (rs2854275 == "CC") {
					DR <- "DR4/x"
					if (rs9273363 == "CC") {
						DQ <- "DR4-DQ7"
						risk <- "low"
					} else {
						DQ <- "DR4-DQ8"
						risk <- "normal"
					}
				} else {
					if (rs2854275 == "AA") {
						DR <- "unknown"
						risk <- "unknown"
						DQ <- "unknown"
					} else {
						DR <- DQ <- "missing"
						risk <- "unknown"
					}
				}
			}
		} else {
			if (rs3104413 == "CC") {
				if (rs2854275 == "AA") {
					DR <- "DR3/3"
					risk <- "high"
					DQ <- "DR3/3-DQ2" #\DRB1*03:01-DQA1*05:01-DQB1*02:01\""
				} else {
					DQ <- "x"
					if (rs2854275 == "AC") {
						DR <- "DR3/x"
						risk <- "normal"
					} else {
						if (rs2854275 == "CC") {
							DR <- "DRx/x"
							risk <- "very low"
						} else {
							DR <- "missing"
							risk <- "unknown"
						}
					}
				}
			} else {
				DR <- "missing"
				risk <- "unknown"
			}
		}
	}
	return(c(DR = DR, DQ = DQ, risk = risk))
}



############### end simple functions #######################

##############
## DATASETS ##
##############

#' Example aSnpMatrix (annotSnpStats) object
#'
#' Dataset. aSnpMatrix object are an extension of a SnpMatrix
#' object from the package snpStats where the sample and SNP annotation
#' is built-in and linked to the SnpMatrix object. Some operations are slower
#' as a result, but this a very robust way of keeping the annotation linked
#' to the data, and allows some powerful functions such as strand alignment
#' to be done using the annotSnpStats package. This is an example object
#' with 100 samples and 20 SNPs. A few samples and SNPs have low call rate,
#' one snp is monomorphic, one has HWE deviation, a few samples have low
#' heterozygosity. So there is good scope for testing functions in this
#' package. 
#' @name exAnnotSnp
#' @docType data
#' @format An object of class aSnpMatrix
#' @keywords datasets
#' @examples
#' exAnnotSnp <- reader("~/github/plumbCNV/exAnnotSnp.rda")
#' data(exAnnotSnp)
#' show(exAnnotSnp)
#' prv(exAnnotSnp@.Data)
#' prv(exAnnotSnp@snps)
#' prv(exAnnotSnp@samples)
NULL # need some sort of code to trigger devtools:document to pick up a dataset


#' Example SnpMatrix (snpStats) object
#'
#' Dataset. An example SnpMatrix object with
#' some LD structure
#' @name exAnnotSnp
#' @docType data
#' @format An object of class aSnpMatrix
#' @keywords datasets
#' @examples
#' exAnnotSnp <- reader("~/github/plumbCNV/exSnpMat.rda")
#' data(exSnpMat)
#' show(exSnpMat)
#' prv(exSnpMat)
NULL # need some sort of code to trigger devtools:document to pick up a dataset



