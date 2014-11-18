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
#parse.sexcheck - read a plink sexcheck file and derive the list of samples to exclude
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
###############


# cbind for more than 2 [a]SnpMatrix objects at once
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

# rbind for more than 2 [a]SnpMatrix objects at once
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

# determine whether a snpMatLst is a list of filenames, SnpMatrix objects, or just a [a]SnpMatrix
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


# select all SNPs for a subset of samples in a snpMatLst
sampSel <- function(snpMatLst,samples,dir=NULL) {
  if(length(samples)>1000) { warning("This function is designed for small lists of 'samples', likely to fail for large datasets")}
  sampselfun <- function(snpMat,samples) {
    if(is.character(samples)) { samples <- narm(match(samples,rownames(snpMat))) }
    if(length(samples)>0) { return(snpMat[samples,]) } else { return(NULL) }
  }
  newlist <- fun.snpMatLst(snpMatLst,fun=sampselfun,dir=dir,samples=samples)
  nulz <- sapply(newlist,is.null)
  newlist <- newlist[!nulz] # remove NULLs from this list or else rbind doesn't work properly
  outobj <- do.call("rbind",args=newlist)  
  if(is.character(samples)) { outobj <- outobj[samples,] }
  if(is.character(samples)) { outobj <- outobj[narm(match(samples,rownames(outobj))),] }
  return(outobj)
}

# select all samples for a subset of SNPs in a snpMatLst
snpSel <- function(snpMatLst,snps,dir=NULL) {
  if(length(snps)>50000) { warning("This function is designed for small lists of 'snps', likely to fail for large datasets")}
  if(!is.character(snps)) { warning("This function is designed for character indexing, results may not be what you expect") }
  snpselfun <- function(snpMat,snps) {
    if(is.character(snps)) { 
      snps <- clean.snp.ids(snps)
      snps <- narm(match(snps,clean.snp.ids(colnames(snpMat)))) 
    }
    if(length(snps)>0) { 
      out <- snpMat[,snps]
      #prv(out)
      return(out) 
    } else { return(NULL) }
  }
  newlist <- fun.snpMatLst(snpMatLst,fun=snpselfun,dir=dir,snps=snps)
  nulz <- sapply(newlist,is.null)
  newlist <- newlist[!nulz] # remove NULLs from this list or else cbind doesn't work properly
 # prv(newlist)
  outobj <- do.call("cbind",args=newlist) 
  colnames(outobj) <- clean.snp.ids(colnames(outobj))
  if(is.character(snps)) { outobj <- outobj[,narm(match(snps,colnames(outobj)))] }
  return(outobj)
}


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




if(F) {

WTCCC2.dir <- "/ipswich/data/chrisw/MS-sawcer/WTCCC2"
subdir <- "CCC2_Cases_UKC/"
fulldir <- cat.path(WTCCC2.dir,fn="",pref=subdir)
dat.fn <- list.files(fulldir)
dat.fn <- cat.path(fulldir,dat.fn[grep(".gen.gz",dat.fn)])
samp.fn <- cat.path("/ipswich/data/chrisw/MS-sawcer/WTCCC2","MS_UKC_illumina.sample",pref=subdir)
ms.dr <- "/chiswick/data/ncooper/imputation/MS/"
sml <- read.impute.list(dat.fn=dat.fn, samp.fn=samp.fn, combine=TRUE,
                            snpcol=1, out.pref="TEMP", out.dir=ms.dr, verbose=TRUE)
}


# read impute file(s) into a SnpMatrix[List] format
# multiple files assumed to be same samples, different SNP sets (e.g, chromsomes)
read.impute.list <- function(dat.fn=dat.fn, samp.fn=samp.fn,
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

# apply a function to each element of a snpMatLst
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
if(F) {
chr.list <- vector("list",22)
for (cc in 1:22) {
  (load(paste0("/ipswich/data/gdxbase/ld/ld_4.18/RDATA/ThousandGenomes_05_11_CEU/genotypes_chr",cc,"_ThousandGenomes_CEU_r0511.RData")))
  keepers <- get.ichip.1000(snps$support)
  SNPS <- snps$genotypes[,rownames(keepers)]
  colnames(SNPS) <- keepers[["ichip"]]
  chr.list[[cc]] <- SNPS
}
sml.1000 <- chr.list
}


get.ichip.1000 <- function(support,build=37) {
  sup <- data.frame.to.ranged(support,chr="chromosome",start="position",end="position",build=build)
  keepers <- subsetByOverlaps(sup,get.support(build=build))
  ii <- match(start(keepers),start(get.support()))
  keepers[["ichip"]] <- rownames(get.support())[ii]
  return(keepers)
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

# check whether a specific number of samples and snps will fit in a SnpMatrix object
check.snpmat.size <- function(nsamp,nsnp) {
  maxsize <- ((2^31)-2)
  if(!is.numeric(nsamp) | !is.numeric(nsnp)) { stop("Error: check.data.size takes numeric args")}
  size <- nsamp*nsnp; pc.max <- round(((size/2^31)*100),1)
  valid <- size<maxsize
  if(!valid) { warning(nsamp," x ",nsnp," is ",pc.max,"% of the allowed object size for a SnpMatrix.\n",
          "Suggest splitting the files into smaller chunks (e.g, by sample subset or chromosomes)") }
  return(valid)
}

# convert a snpMatLst to a single SnpMatrix #
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

# get the dimensions of all of the SnpMatrix objects in a snpMatLst, as a matrix, where cols are datasets, first row is nsamp, second nsnp
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


# get a list of the 'colnames' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, combined set of colnames)
colnamesL <- function(snpMatLst,list=FALSE,dir=NULL) {
  alln <- smlapply(snpMatLst,FUN=colnames,dir=dir)
  if(list) { return(alln) }
  allN <- unlist(alln); if(length(snpMatLst)==1) { return(allN) }
  if(anyDuplicated(allN)) { 
    if(all(alln[[1]]==alln[[2]])) { return(unlist(alln[[1]])) } # are all elements the same snps
    warning("duplicate (column) SNP names found across SnpMatrix objects in SnpMatrixList")   
  }
  return(allN)
}

# get a list of the 'rownames' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, combined set of rownames)
rownamesL <- function(snpMatLst,list=FALSE,dir=NULL) {
  alln <- smlapply(snpMatLst,FUN=rownames,dir=dir)
  if(list) { return(alln) }
  allN <- unlist(alln); if(length(snpMatLst)==1) { return(allN) }
  if(anyDuplicated(allN)) { 
    if(all(alln[[1]]==alln[[2]])) { return(unlist(alln[[1]])) } # are all elements the same samples
    warning("duplicate (row) sample names found across SnpMatrix objects in SnpMatrixList") 
  }
  return(allN)
}

# get a list of the 'nrow' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, total nrow)
nrowL <- function(snpMatLst,list=TRUE,dir=NULL) {
  out <- smlapply(snpMatLst,FUN=nrow,c.fun=unlist,dir=dir)
  if(!list)  { 
    oo <- unique(out)
    if(length(oo)==1) { out <- oo } else  { out <- sum(out) }
  }
  return(out)
}

# get a list of the 'ncol' for all SnpMatrix objects in a snpMatLst (or if list=FALSE, total ncol)
ncolL <- function(snpMatLst,list=TRUE,dir=NULL) {
  out <- smlapply(snpMatLst,FUN=ncol,c.fun=unlist,dir=dir)
  if(!list)  { 
    oo <- unique(out)
    if(length(oo)==1) { out <- oo } else  { out <- sum(out) }
  }
  return(out)
}


# 'lapply' function for SnpMatrixList, also works for normal SnpMatrix or aSnpMatrix or XSnpMatrix
#' @param c.fun function to combine the result, noting that the raw form of the result is a list, so
#' if the result should be a vector, then the appropriate combination function would be 'unlist'.
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


# read a plink sexcheck file and derive the list of samples to exclude
parse.sexcheck <- function(fn, excl=0.5, write=TRUE) {
  sc <- read.table(fn,header=T)
  sc <- sc[sc$STATUS!="OK",]
  sc <- sc[!is.na(sc$F),]
  somefails <- ((sc$PEDSEX==1 & sc$SNPSEX==2) | (sc$PEDSEX==2 & sc$SNPSEX==1))
  morefails <- sc$F>excl
  failers <- paste(sc[somefails | morefails,"IID"])
  if(write & length(failers)>0) {
    out.fn <- cat.path(dirname(fn),basename(fn),suf="fail",ext="txt")
    writeLines(failers,con=out.fn)
    cat("wrote samples failing sexcheck to",out.fn,"\n")
  }
  return(failers)
}


# run sample quality control on a SnpMatrix or snpMatLst
doSampQC <- function(snpMatLst=NULL, dir=getwd(), subIDs.plink=NULL, plink=FALSE, call.rate=.95, het.lo=0.1, het.hi=0.4, n.cores=1, proc=1) 
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


# run SNP quality control on a SnpMatrix or snpMatLst
doSnpQC <- function(snpMatLst=NULL, snp.info=NULL, sample.info=NULL, dir=getwd(), plink=FALSE, n.cores=1,
                    call.rate=.95, hwe.p.thr=(10^-7), grp.hwe.z.thr=5.32, grp.cr.thr=.001, group.miss=F,
                    maf.thr=0.005, autosomes.only=T, subIDs.plink=NULL, proc=1) 
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
    } else { stop("if snpMatLst does not contain aSnpMatrix objects, then you must enter argument 'snp.info'")  }
  }
  if(is.null(sample.info)) {
    if(is.aSnpMatrix(snpMatLst,dir=dir)) { 
      sample.info <- sample.info.from.annot(snpMatLst,dir=dir)
    } else {
      sample.info <- sample.info.from.sml(snpMatLst, subIDs.actual=subIDs.actual, dir=dir, plink=plink)
    }
  }
  ## get snp call rates ##
  HD <- switch(snp.mat.list.type(snpMatLst,fail=F),memory=F,disk=T,error=NULL)
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


# apply the thresholds from doSnpQC() to a dataset to get exclusion lists
#' @param maf, numeric, threshold for maf... uses 'less-than', except for zero, which
#' uses 'equal-to'. If you don't want to exclude any SNPs based on maf, use 'NA'.
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




# snpStats::row.summary for a snpMatLst (sample QC summary)
list.rowsummary <- function(snpMatLst,mode="row",dir=getwd(),warn=T,n.cores=1)
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
    cat(" processing elements in",min(n.cores,length(snpMatLst)),"parallel cores ...")
    rowsum.list <- mclapply(snpMatLst,my.summary)
  } else {
    cat(" processing element ")
    for (dd in 1:length(snpMatLst))
    {
      cat(dd,"..",sep="")
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
  cat("done\n")
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

# snpStats::col.summary for a snpMatLst (SNP-qc sumamry)
list.colsummary <- function(snpChrLst,mode="col",dir=getwd(),warn=F,n.cores=1)
{
  # wrapper to make 'list.rowsummary' work for 'col.summary' too
  # warn=F helps to avoid alarming user given known issue with genotypes
  # that are uniformly zero (empty) - chip qc snps, etc
  if(mode!="col" & mode!="row") { mode <- "col" }
  return(list.rowsummary(snpChrLst,mode=mode,dir=dir,warn=F,n.cores=n.cores))
}

# convert snpMatLst which is split by sample subsets to one split by chromosome
convert.smp.to.chr22 <- function(snpMatLst,snp.info,dir="",n.cores=1)
{
  ## convert snp.matrix list separated into different sample sets
  ## with all chromosomes, into 22 chromosome lists with all samples
  ## HD=T saves chromosome list elements to disk (loc) instead of memory
  ## in a situation where data is very large or memory is limited
  ## In HD=T the output will be a list of disk locations rather than an array
  dir <- validate.dir.for(dir,"cr")
  #print(snp.mat.list.type(snpMatLst,fail=T))
  HD <- switch(snp.mat.list.type(snpMatLst,fail=T),memory=F,disk=T,error=NULL)
  chr.set <- chrNums(snp.info); n.chr <- length(chr.set)
  #print(chr.set)
  cat("\nConverting from 'n' group list of all markers to",n.chr,"chromosome list of all samples\n")
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
  cat("\n"); print(list.spec); cat("\n")
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
  one.chr <- function(cc) {
    # do all the procesing for one chromosome
    # snp.info is a RangedData object (IRanges)
    next.chr <- rownames(snp.info[cc]) 
    next.mat <- matrix(raw(),nrow=num.subs,ncol=length(next.chr))
    colnames(next.mat) <- next.chr; rownames(next.mat) <- rnm 
    next.snpmat <- new("SnpMatrix", next.mat)
    rm(next.mat) ; gc()
    for (dd in 1:length(snpMatLst))
    {
      cat(".")
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
      fnm <- cat.path(dir$cr,"chrmat_",suf=cc,ext="RData")
      save(next.snpmat,file=fnm)
      return(fnm)
    } else {
      return(next.snpmat)
    }
  }
  ## run for each chromosome in serial or parallel mode
  if(n.cores>1) {
    cat(" processing chromosomes",paste(range(chr.set),collapse="-"))
    must.use.package("parallel"); snpChrLst <- mclapply(1:n.chr,one.chr,mc.cores=min(n.cores,4))
    cat(" .. done\n")
  } else {
    cat(" processing chromosome "); snpChrLst <- vector("list",n.chr)
    for (j in 1:n.chr) {
      cat(chr.set[j]); snpChrLst[[j]] <- one.chr(j)
    }
    cat("..done\n")
  }
  return(snpChrLst)
}


# convert snpMatLst which is split by chromsome to one split by sample subsets
convert.chr22.to.smp <- function(snpChrLst,snp.info,subIDs.actual,HD=F,group.nums,dir=NULL)
{
  ## convert snp.matrix list separated into 22+ chromosome lists with all samples
  ## into different sample sets with all chromosomes. 
  ## HD=T saves chromosome list elements to disk (loc) instead of memory
  ## in a situation where data is very large or memory is limited
  ## In HD=T the output will be a list of disk locations rather than an array
  HD <- switch(snp.mat.list.type(snpChrLst,fail=T),memory=F,disk=T,error=NULL)
  if (is(snp.info)[1]=="RangedData")
  { must.use.package("genoset",bioC=T)
  } else { stop("snp.info object was not of type 'RangedData'")}
  chr.set <- chrNums(snp.info); n.chr <- length(chr.set)
  cat("\nConverting from",n.chr,"chr list of all samples to 'n' group list of all markers\n")
  num.grps <- length(unique(group.nums))
  snpMatLst <- vector("list", num.grps)
  num.subs <- length(subIDs.actual)
  grp.sizes <- table(group.nums)
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
    cat(paste(" writing group",dd,"\n"))
    next.mat <- matrix(raw(),nrow=grp.sizes[dd],ncol=num.snps)
    colnames(next.mat) <- rownames(snp.info)
    next.grp <- subIDs.actual[group.nums==dd]
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
        cat(" processing chromosome ",chr.set[cc],"...",sep="")
        row.select <- next.grp %in% rownames(snpMatLst[[cc]])
        next.snpmat[row.select,st[cc]:en[cc]] <- snpChrLst[[cc]][next.grp[row.select],]
        cat("..done\n")
      }
    }
    if(HD) {
      fnm <- cat.path(dir$cr,"snpmat_",suf=dd,ext="RData")
      save(next.snpmat,file=fnm)
      snpMatLst[[dd]] <- fnm
    } else {
      snpMatLst[[dd]] <- next.snpmat
    }
    cat(" complete!\n")
  }
  return(snpMatLst)
}

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
      if(any(mis=="aSnpMatrix")) { return(TRUE) } else { return(FALSE) }
    } else {
      if(all(mis=="aSnpMatrix")) { return(TRUE) } else { return(FALSE) }
    }
  } else {
    return(FALSE)
  }
}


# create a full sample.info object from a list of aSnpMatrix objects
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

# create a full snp.info object from a list of aSnpMatrix objects
snp.info.from.annot <- function(aSnpMatLst,dir=getwd()) {
  if(!is.aSnpMatrix(aSnpMatLst,dir=dir)) { stop("not an aSnpMatrix object") }
  if(is.list(aSnpMatLst)) {
    typ <- get.split.type(aSnpMatLst,dir=dir)
    if(typ=="snpSubGroups") {
      outlist <- smlapply(aSnpMatLst,FUN=snp.from.annot,dir=dir)
      out <- do.call("rbind",args=outlist)
      return(out)
    } else {
      return(snp.from.annot(snpMatLst.collate(aSnpMatLst[1],dir=dir)))
    }
  } else {
    return(snp.from.annot(aSnpMatLst))
  }
}

# create a (partial) snp.info object from a single aSnpMatrix object
snp.from.annot <- function(aSnpMat) {
  if(!is(aSnpMat)[1]=="aSnpMatrix") { stop("aSnpMat must be of type 'aSnpMatrix'") }
  df <- aSnpMat@snps
  cn <- colnames(df)
  chr.col <- which(tolower(substr(cn,1,2)) %in% "ch")
  pos.col <- which(tolower(substr(cn,1,2)) %in% "po")
  DF <- data.frame.to.ranged(df,start=cn[pos.col],end=cn[pos.col],chr=cn[chr.col],GRanges=FALSE)
  DF[["QCfail"]] <- rep(0,nrow(DF))
  return(DF)
}

# create a (partial) sample.info object from a single aSnpMatrix object
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


correct.by <- function(X) {
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix")) { stop("X must be a SnpMatrix object") }
  cs <- col.summary(X)
  if(nrow(cs)!=ncol(X)) { stop("col.summary length (",nrow(cs),") did not match length of X (",ncol(X),")")}
  raf <- cs$RAF
  nn <- cs$Calls
  aa <- raf^2
  ab <- raf*(1-raf)*2
  bb <- (raf)^2
  mN <- (aa*1) + (ab*2) + (bb*3)
  ss <- ((bb*nn)*(3-mN)^2) + ((ab*nn)*(2-mN)^2) + ((aa*nn)*(1-mN)^2)
  sD <- sqrt((1/(nn-1))*ss)  
  sD[is.na(sD)] <- 1
  sD[sD<.0001] <- 0.1
  return(list(mean=mN,sd=sD))
}

bigSnpMatrix <- function(X,filename="tempMatrix",tracker=TRUE, limit.ram=FALSE, mn.zero=TRUE, sd.hwe=TRUE, replace.missing=TRUE) {
  typ <- is(X)[1]
  max.gb <- NA
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) {
      stop("X must be a SnpMatrix or aSnpMatrix object")
  }
  if(replace.missing) {
    X <- randomize.missing2(X)
  }
  raw <- X@.Data
  if(mn.zero | sd.hwe) {
    oo <- correct.by(X); mN <- oo$mean; sD <- oo$sd
  } 
  #raw[raw==0] <- NA
  if(!mn.zero) { mN <- 0 }
  if(!sd.hwe) { sD <- 1 }
  nC <- ncol(X); nR <- nrow(X)
  dir <- dirname(filename); filename <- basename(filename)
  if(dir=="") { dir <- getwd() }
  if(tracker) { cat(" creating",nR,"x",nC,"target matrix,",cat.path(dir,filename),"...\n") }
  des <- paste(filename,"dsc",sep=".")
  bck <- paste(filename,"bck",sep=".")
  bigSnp <- big.matrix(nrow=nR,ncol=nC, backingfile=bck, dimnames=list(rownames(X),colnames(X)),
                         backingpath=dir, descriptorfile=des)
  split.to <- 10*round(estimate.memory(bigSnp)) # split into .1GB chunks, save RAM without creating groups too small to process
  #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
  stepz <- round(seq(from=1,to=nC+1,length.out=round((split.to+1))))
  if((tail(stepz,1)) != nC+1) { stepz <- c(stepz,nC+1) }
  split.to <- length(stepz)-1
  for (cc in 1:split.to)
  {
    # within submatrix cols
    c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
    # do the copying
    lilColRange <- c(c1:c2)
    if(tracker) {      loop.tracker(cc,split.to) }
    if(is.finite(sum(lilColRange))) {
      #cat(range(lilColRange)); cat(dim(bigTrans)); cat(dim(bigMat))
      nxt.chunk <- raw[1:nR,lilColRange]
      nxt.chunk <- apply(nxt.chunk,2,as.numeric)
      bigSnp[1:nR,lilColRange] <- (nxt.chunk-mN[lilColRange])/sD[lilColRange]
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
  if(tracker) {
    prv(bigSnp)
  }
  return(bigSnp)
}


# sml <- list.files("/chiswick/data/ncooper/iChipData/",pattern="temp.ichip")
# sml <- as.list(sml)
# dd <- "/chiswick/data/ncooper/iChipData/"

# convert a snpMatLst to a single SnpMatrix #
bigSnpMatrixList <- function(snpMatLst, verbose=TRUE, dir=getwd(),replace.missing=TRUE,
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
        oo <- correct.by(next.chr); mN <- oo$mean; sD <- oo$sd
      } 
      if(!mn.zero) { mN <- 0 }
      if(!sd.hwe) { sD <- 1 }
      #########
      split.to <- 10*round(estimate.memory(next.chr)) # split into .1GB chunks, save RAM without creating groups too small to process
      #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
      nr <- nrow(next.chr); nc <- ncol(next.chr)
      stepz <- round(seq(from=1,to=nc+1,length.out=round((split.to+1))))
      if((tail(stepz,1)) != nc+1) { stepz <- c(stepz,nc+1) }
      split.to <- length(stepz)-1
      for (cc in 1:split.to)
      {
        # within submatrix cols
        c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
        # do the copying
        lilColRange <- c(c1:c2)
        if(is.finite(sum(lilColRange))) {
          #cat(range(lilColRange)); cat(dim(bigTrans)); cat(dim(bigMat))
          nxt.chunk <- raw[1:nr,lilColRange,drop=F]
          nxt.chunk <- apply(nxt.chunk,2,as.numeric)
          if(typ2=="snpSubGroups") {
#            prv(stz[dd],enz[dd])
            #return(list(bigSnp=bigSnp,stz=stz,enz=enz,dd=dd,nr=nr,lilColRange=lilColRange,nxt.chunk=nxt.chunk,mN=mN,sD=sD))
            bigSnp[1:nr,(stz[dd]:enz[dd])[lilColRange]] <- (nxt.chunk-mN[lilColRange])/sD[lilColRange]  
          } else {
            bigSnp[stz[dd]:enz[dd],lilColRange] <- (nxt.chunk-mN[lilColRange])/sD[lilColRange]
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

ld.prune.chr <- function(X,stats="R.squared",thresh=.50) {
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
