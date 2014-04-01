


apply.snp.thresholds <- function(snp.info,callrate.snp.thr=0.95,hwe.thr=10^-5,
                                 grp.cr.thr=.001,grp.hwe.z.thr=4,proc=1) {
  ## now remove samples failing HWE + callrate 95 regardless of source
  # also optionally remove those failing group criteria if columns present in snp.info
  snp.info[["QCfail"]] <- rep(0,times=nrow(snp.info))
  snp.info$call.rate[is.na(snp.info$call.rate)] <- 0
  snp.info$P.hwe[is.na(snp.info$P.hwe)] <- 0
  snp.info$Z.hwe[is.na(snp.info$Z.hwe)] <- 0
  cr.cond <- snp.info$call.rate < callrate.snp.thr
  call.rate.excl.snps <- row.names(snp.info)[cr.cond]
  HWE.cond <- snp.info$P.hwe < hwe.thr
  HWE.exclude.snps <- row.names(snp.info)[HWE.cond]
  to.remove <- unique(c(call.rate.excl.snps,HWE.exclude.snps))
  out.list <- list(snp.info,call.rate.excl.snps,HWE.exclude.snps,cr.cond,HWE.cond)
  names(out.list) <- c("SNP.INFO","CR.EXCL","HWE.EXCL","CR.cond","HWE.cond")
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
                     grp.hwe.exclude.snps,cr.cond,HWE.cond,grp.cr.cond,grp.hwe.cond)
    names(out.list) <- c("SNP.INFO","CR.EXCL","HWE.EXCL","GRP.CR.EXCL","GRP.HWE.EXCL","CR.cond",
                         "HWE.cond","grp.CR.cond","grp.HWE.cond")
  }
  if(length(to.remove)>0)  {
    ## mark those failing on any threshold
    snp.info[["QCfail"]][match(to.remove, rownames(snp.info))] <- proc  }
  snp.info$QCfail[is.na(snp.info$QCfail)] <- 1
  out.list$SNP.INFO <- snp.info
  return(out.list)
}



call.rate.summary <- function(s.info=NULL,snp=(is(s.info)[1]=="RangedData"),cr.vec=NULL,print=T) 
{
  # generate summary of callrate performance for snps/samples
  cat("\nGenerating report")
  if(snp) { cat(" (snp)\n") } else { cat(" (sample)\n") }
  crs.o <- NULL; if(!is.null(cr.vec)) { crs.o <- cr.vec } 
  if(!is.null(s.info)) { if ("call.rate" %in% colnames(s.info)) {
    crs.o <- s.info$call.rate
  } else {
    s.info <- column.salvage(s.info,"call.rate",c("Call.rate","callrate"),ignore.case=T)
    crs.o <- s.info$call.rate
  }}
  if(is.null(crs.o)) { warning("need an object with column 'call.rate' or cr.vec of callrates") ; return(NULL) }
  crs <- narm(crs.o)
  # Report
  thresh.lt <- c(.5,.75,.9,.95)
  thresh.gt <- c(.95,.97,.99,.999)
  counts.lt <- integer(length(thresh.lt))
  counts.gt <- integer(length(thresh.gt))
  # samples or snps (autodetected)
  for (cc in 1:length(thresh.lt))
  { counts.lt[cc] <- length(crs[crs<thresh.lt[cc]]) }
  for (cc in 1:length(thresh.gt))
  { counts.gt[cc] <- length(crs[crs>=thresh.gt[cc]]) }
  denom <- length(crs.o)  
  counts <- c(counts.lt, counts.gt)
  counts.pc <- round(counts/denom*100,2)
  rn.resul <- c(paste("callrate <",thresh.lt),paste("callrate >=",thresh.gt))
  if(snp) {
    result <- data.frame(CallRate=rn.resul,SNPs=counts,Percent=paste(counts.pc,"%",sep=""))
  } else {
    result <- data.frame(CallRate=rn.resul,Samples=counts,Percent=paste(counts.pc,"%",sep=""))
  }
  if(print) { print(result); cat("\n") }
  return(result)
}


colSummary <- function(X,filt=NULL) {
  ## do a snpStats::col.summary() only on selected set of snps
  must.use.package("snpStats",T)
  if(is(X)[1]=="SnpMatrix")  {
    if(is.character(filt) & length(filt)>0) {
      select <- which(colnames(X) %in% filt)
      if(length(select)>0) {
        return(col.summary(X[,select]))
      } else {
        cat("No snp ids from filt in colnames of X\n"); return(col.summary(X))
      }
    } else {
      return(col.summary(X))
    }
  } else {
    warning("not a SnpMatrix, returning NULL");  return(NULL)
  }
}


# detect the format of a snpMatLst (SnpMatrix list)
# returns snpSubGroups, sampleSubGroups, or singleEntry (else error)
get.split.type <- function(snpMatLst,dir=NULL) {
  typz <- sapply(snpMatLst,is)[1,]
  dir <- validate.dir.for(dir,c("cr"))
  if(all(typz=="SnpMatrix"))
  { HD <- F } else {
    if (all(typz=="character")) { HD <- T } else {
      stop("snpMatLst doesn't appear to be a list of SnpMatrix objects",
              "or a list of file locations, row.summary impossible")
    } 
  }
  dimz <- get.snpMat.spec(snpMatLst,dir=dir)
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
    if(!same.num.snps & !same.num.samples) { stop("invalid/incomplete snpMatLst object") }
  }
}



list.rowsummary <- function(snpMatLst,mode="row",dir=getwd(),warn=T,n.cores=1)
{
  # performs 'row.summary' or 'col.summary' snpStats functions
  # on a list of snpMatrix objects
  fail <- F
  typz <- sapply(snpMatLst,is)[1,]
  dir <- validate.dir.for(dir,c("cr"))
  if(all(typz=="SnpMatrix"))
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
  if(mat.type==mat.typez[3]) { return(my.summary(snpMatLst[[1]])) } # a list length 1!
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
        fl.nm <- cat.path(dir$cr,snpMatLst[[dd]])
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
      warning("sample-rownames were corrupted, fixing:\n",
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


list.colsummary <- function(snpChrLst,mode="col",dir=getwd(),warn=F,n.cores=1)
{
  # wrapper to make 'list.rowsummary' work for 'col.summary' too
  # warn=F helps to avoid alarming user given known issue with genotypes
  # that are uniformly zero (empty) - chip qc snps, etc
  if(mode!="col" & mode!="row") { mode <- "col" }
  return(list.rowsummary(snpChrLst,mode=mode,dir=dir,warn=F,n.cores=n.cores))
}


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
  rownames(list.spec) <- c("Samples","SNPs"); colnames(list.spec) <- basename(colnames(list.spec))
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



het.density.plots <- function(data="sampleqc.txt",het.lo=.1,het.hi=0.4,zoom=T,
                              het=NULL,dir=NULL,fn="HZDistribution.pdf",...) {
  # make density plot of Heterozygosity; can enter a file location/data.frame or vecs of heterozygosity
  if(!is.null(dir)) { dir <- validate.dir.for(dir,"cr") }
  if(!is.null(data) & is.null(het)) { 
    val <- F; if(!is.null(dim(data))) { val <- T }
    if(!is.null(dir)) { if(is.file(data,dir$cr,dir)) { 
      data <- find.file(data,dir$cr,dir); val <- T } }
    if(val) {
      tt <- force.frame(data) 
      tt <- column.salvage(tt,"heterozygosity",c("Heterozygosity","Het","HZ","Hzg","Hetz"),T)
      if(all(c("heterozygosity") %in% colnames(tt))) {
        het <- tt$heterozygosity
      } else {
        warning(paste("required column: heterozygosity not found in",data))
      }
    }
  }
  if(is.null(het)) { warning("no heterozygosity data found") ; return(NULL) }
  add.boundary.and.legend <- function(het.lo=NULL,het.hi=NULL) {
    if(!is.null(het.lo) & !is.null(het.hi)) { 
      abline(v=het.lo,lty="dashed",col="blue") 
      abline(v=het.hi,lty="dashed",col="blue")
      legend("topright",legend=c(paste("Cutoffs, ",het.lo,"< Hz <",het.hi,sep="")),
             lty="dashed",col="blue",bg="white",box.lwd=0)
    }
  }
  if(!is.null(dir)) { ofn <- cat.path(dir$cr,fn); pdf(ofn) }
  if(zoom) { xl <- c(0,0.6) }
  par(mfrow=c(1,1))
  plot(density(het),main="Heterozygosity distribution across all samples",
       xlab="Heterozygosity (Hz) score",bty="n",xlim=xl,...)
  add.boundary.and.legend(het.lo,het.hi)
  if(!is.null(dir)) { dev.off() }
}


hwe.density.plots <- function(data="snpqc.txt",hwe.thr=10^-5,zoom=T,
                              Z.hwe=NULL,dir=NULL,fn="HWEDistribution.pdf",...) {
  # make density plot of HWE; can enter a file location/data.frame or vecs of hwe/
  if(!is.null(dir)) { dir <- validate.dir.for(dir,"cr") }
  if(!is.null(data) & is.null(Z.hwe)) { 
    val <- F; if(!is.null(dim(data))) { val <- T } else { if(is.null(dir)) { dir <- getwd() }}
    if(!is.null(dir)) { if(is.file(data,dir$cr,dir)) { 
      data <- find.file(data,dir$cr,dir); val <- T } }
    if(val) {
      tt <- force.frame(data) 
      tt <- column.salvage(tt,"z.hwe",c("z.hwe","hwe.z","zhwe","hwez","hwe"),T)
      if(all(c("z.hwe") %in% colnames(tt))) {
        Z.hwe <- tt$z.hwe
      } else {
        warning(paste("required column: z.hwe not found in",data))
      }
    }
  }
  if(is.null(Z.hwe)) { warning("no Z.hwe data found") ; return(NULL) }
  add.boundary.and.legend <- function(hwe.thr=NULL) {
    if(!is.null(hwe.thr)) { 
      abline(v=c(-1,1)*qnorm(hwe.thr/2),lty="dashed",col="blue") 
      legend("topright",legend=c(paste("HWE cutoff, p<",hwe.thr,sep="")),
             lty="dashed",col="blue",bg="white",box.lwd=0)
    }
  }
  if(!is.null(dir)) { ofn <- cat.path(dir$cr,fn); pdf(ofn) }
  if(zoom) { xl <- c(-10,10) }
  mms <- length(which(is.na(Z.hwe)))
  DD <- density(narm(Z.hwe))
  par(mfrow=c(1,1))
  plot(DD,main="Z-score distribution across all SNPs",
       xlab="HWE Z-score",bty="n",xlim=xl,...)
  if(length(mms>0)) { cat("HWE data for",mms,"monomorphic SNPs was ignored\n")}
  add.boundary.and.legend(hwe.thr=hwe.thr)
  if(!is.null(dir)) { dev.off() }
}



hwe.vs.callrate.plots <- function(data="snpqc.txt",callrate.snp.thr=.95,hwe.thr=10^-5,
                                  Z.hwe=NULL,call.rate=NULL,zoom=T,full=T,dir=NULL,
                                  fn="HWEvsCallrate.pdf",excl=F, incl=F) {
  # make plots of HWE against callrate; can enter a file location/data.frame or vecs of hwe/
  if(!is.null(dir)) { dir <- validate.dir.for(dir,"cr") }
  if(!full & !zoom) { return(NULL) }
  if(!is.null(data) & ((is.null(Z.hwe) | is.null(call.rate)))) { 
    val <- F; if(!is.null(dim(data))) { val <- T }
    if(!is.null(dir)) { if(is.file(data,dir$cr,dir)) { 
      data <- find.file(data,dir$cr,dir); val <- T } }
    if(val) {
      tt <- force.frame(data) 
      tt <- column.salvage(tt,"z.hwe",c("z.hwe","hwe.z","zhwe","hwez","hwe"),T)
      tt <- column.salvage(tt,"call.rate",c("Call.rate","callrate","cr"),T)
      if(all(c("z.hwe","call.rate") %in% colnames(tt))) {
        Z.hwe <- tt$z.hwe; call.rate <- tt$call.rate
      } else {
        warning(paste("required columns: ",paste(c("z.hwe"," call.rate"),collapse=","),", not found in ",data,sep=""))
      }
    }
  }
  if(is.null(Z.hwe) | is.null(call.rate)) { warning("no Z.hwe and/or call.rate data found") ; return(NULL) }
  if(length(Z.hwe)!=length(call.rate)) { cat(" HWE and callrate vectors had unequal length\n"); return(NULL) }
  add.boundary.and.legend <- function(callrate.snp.thr=NULL,hwe.thr=NULL,scl=1) {
    if(!is.null(hwe.thr)) { abline(v=c(-1,1)*qnorm(hwe.thr/2),lty="solid",col="blue") }
    if(!is.null(callrate.snp.thr)) { abline(h=callrate.snp.thr,lty="dashed",col="blue") }
    if(!is.null(hwe.thr) & !is.null(callrate.snp.thr)) {
      legend("topright",legend=c(paste("HWE cutoff, p<",hwe.thr,sep=""),
                                 paste("Call rate cutoff, ",callrate.snp.thr,"%",sep="")),
             lty=c("solid","dashed"),col="blue",bg="white",box.col="white",box.lty="dotted",cex=scl)
    }
  }
  if(!is.null(dir)) { ofn <- cat.path(dir$cr,fn); scl <- 1; pdf(ofn) 
  } else { if(zoom & full) { par(mfrow=c(1,2)); scl <- 0.6 } else { scl <- 1 } }
  if(full) {
    plot(Z.hwe,call.rate,pch=".",xlab="HWE Z-score",ylab="call rate",
         main="full range",ylim=c(1,0),bty="n")
    add.boundary.and.legend(callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,scl=scl)
  }
  if(zoom) {
    plot(Z.hwe,call.rate,pch=".",xlim=c(-10,10),ylim=c(1,callrate.snp.thr-.01),
         xlab="HWE Z-score",ylab="call rate",main="cutoff range",bty="n")
    add.boundary.and.legend(callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,scl=scl)
  }
  if(!is.null(dir)) { dev.off() }
  if(excl | incl) {
    # return excluded sample list
    filt <- which(!(abs(Z.hwe)<hwe.thr & call.rate>callrate.snp.thr))
    if(incl & excl) {
      return(list(incl=!filt,excl=filt))
    } else {
      if(excl) { return(filt) } else { return(!filt) }
    }
  }
}





hz.vs.callrate.plots <- function(data="snpqc.txt",callrate.samp.thr=.95,het.lo=.1,het.hi=0.4,
                                 het=NULL,call.rate=NULL,zoom=T,full=T,dir=NULL,
                                 fn="HZvsCallrate.pdf",snp.info=NULL,excl=F, incl=F) {  
  # make plots of HWE against callrate; can enter a file location/data.frame or vecs of hwe/
  if(!is.null(dir)) { dir <- validate.dir.for(dir,"cr") }
  if(!full & !zoom) { return(NULL) }
  if(!is.null(data) & ((is.null(het) | is.null(call.rate)))) { 
    val <- F; if(!is.null(dim(data))) { val <- T }
    if(!is.null(dir)) { if(is.file(data,dir$cr,dir)) { 
      data <- find.file(data,dir$cr,dir); val <- T } }
    if(val) {
      tt <- force.frame(data) 
      tt <- column.salvage(tt,"call.rate",c("Call.rate","callrate","cr"),T)
      tt <- column.salvage(tt,"heterozygosity",c("Heterozygosity","Het","HZ","Hzg","Hetz"),T)
      if(all(c("heterozygosity","call.rate") %in% colnames(tt))) {
        het <- tt$het; call.rate <- tt$call.rate
      } else {
        warning(paste("required columns: ",paste(c("heterozygosity"," call.rate"),collapse=","),", not found in ",data,sep=""))
      }
    }
  }
  if(is.null(het) | is.null(call.rate)) { warning("no heterozygosity and/or call.rate data found") ; return(NULL) }
  if(length(het)!=length(call.rate)) { cat(" het and callrate vectors had unequal length\n"); return(NULL) }
  add.boundary.and.legend <- function(callrate.samp.thr=NULL,het.lo=NULL,het.hi=NULL,scl=1) {
    if(!is.null(het.lo)) { abline(v=het.lo,lty="solid",col="blue") }
    if(!is.null(het.hi)) { abline(v=het.hi,lty="solid",col="blue") }
    if(!is.null(callrate.samp.thr)) { abline(h=callrate.samp.thr,lty="dashed",col="blue") }
    if(!is.null(het.hi) & !is.null(het.lo) & !is.null(callrate.samp.thr)) {
      legend("topright",legend=c(paste("Hz cutoffs, ",het.lo,"< Hz <",het.hi,sep=""),
                                 paste("Call rate cutoff, ",callrate.samp.thr,"%",sep="")),
             lty=c("solid","dashed"),col="blue",bg="white",box.col="white",box.lty="dotted",cex=scl)
    }
  }
  if(!is.null(dir)) { ofn <- cat.path(dir$cr,fn); scl <- 1; pdf(ofn) 
  } else { if(zoom & full) { par(mfrow=c(1,2)); scl <- 0.6 } else { scl <- 1 } }
  if(full) {
    plot(het,call.rate,pch=".",xlab="Heterozygosity (Hz)",ylab="call rate",
         main="full range",ylim=c(1,0),bty="n")
    add.boundary.and.legend(callrate.samp.thr=callrate.samp.thr,
                            het.lo=het.lo,het.hi=het.hi,scl=scl)
  }
  if(zoom) {
    zoom.range <- mean(het,na.rm=T)+(c(-2,2)*sd(het,na.rm=T))
    plot(het,call.rate,pch=".",xlim=zoom.range,ylim=c(.005+max(call.rate,na.rm=T),callrate.samp.thr-.01),
         xlab="Hz-score",ylab="call rate",main="cutoff range",bty="n")
    add.boundary.and.legend(callrate.samp.thr=callrate.samp.thr,
                            het.lo=het.lo,het.hi=het.hi,scl=scl)
  }
  if(!is.null(dir)) { dev.off() }
  if(excl | incl) {
    # return excluded sample list
    filt <- which(!(het>het.lo & het<het.hi & call.rate>callrate.samp.thr))
    if(incl & excl) {
      return(list(incl=!filt,excl=filt))
    } else {
      if(excl) { return(filt) } else { return(!filt) }
    }
  }
}


