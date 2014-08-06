
get.PCA.subset <- function(dir,pc.to.keep=.13,assoc=F,autosomes=T,big.fn="combinedBigMat.RData",
                           snp.sub.fn="pca.snp.subset.txt",use.current=F,pref="PCAMatrix",
                           descr.fn="pcaSubMat.RData",nprev=0,snp.info=NULL,sample.info=NULL) 
{
  ## extract LRR matrix with subset of SNPs, ready for PCA analysis
  dir <- validate.dir.for(dir,c("ano","big"))
  load.all.libs() # load all main libraries used by plumbCNV
  #if(add.pheno) {
  #  sample.info <- read.sample.info(dir)
  #  phenotype <- read.table(file=dir.fn.cmb(dir$ano,"pheno.lookup.txt"))
  #  sample.info <- add.to.sample.info(sample.info,phenotype,"phenotype")
  #  write.sample.info(sample.info,dir)
  #}
  if(!is.data.frame(sample.info)) { sample.info <- read.sample.info(dir,nprev=nprev) }
  if(is(snp.info)[1]!="RangedData") { snp.info <- read.snp.info(dir,nprev=nprev) }
  #sample.info <- validate.samp.info(sample.info,QC.update=F,verbose=F) #this done done later anyway
  #samp.fn <- "combined.samples.txt"
  if(use.current & is.valid.file(snp.sub.fn,dir$ano,dir)) {
    snps.to.keep <- get.vec.multi.type(snp.sub.fn,dir=dir)
  } else {
    snps.to.keep <- extract.snp.subset(snp.info,sample.info,pc.to.keep=pc.to.keep,assoc=assoc,autosomes=autosomes,
                                       writeResultsToFile=T,big.fn=big.fn,out.fn=snp.sub.fn,dir=dir)
  }
  ###bigMat <- getBigMat(big.fn,dir)
  if(length(snps.to.keep)>100) {
    ##writeLines(colnames(bigMat),paste(dir$ano,samp.fn,sep=""))
    subset.descr <- big.exclude.sort(big.fn,dir=dir,T,tranMode=1,pref=pref,f.snp=snps.to.keep,verbose=F)
  } else {
    stop("Error: list of snps to keep is too small - trying running again with a higher pc.to.keep\n")
  }
  #if(descr.fn!="") {
  #  save(subset.descr,file=dir.fn.cmb(dir$big,descr.fn))
  #} else { warning("submatrix description returned but not saved\n")}
  print(subset.descr)
  return(subset.descr)
}


big.PCA <- function(subDescr,dir,pcs.to.keep=50,SVD=F,LAP=F,save.pcs=T,pcs.fn="PCsEVsFromPCA.RData") 
{
  # run principle components analysis on the SNP subset of the LRR snp x sample matrix
  # various methods to choose from with pro/cons of speed/memory, etc.
  #  must use SNP-subset to avoid LD, destroying main effects, +avoid huge memory requirements
  dir <- validate.dir.for(dir,c("big","pc"))
  pcaMat <- getBigMat(subDescr,dir)
  cat("\nRunning Principle Components Analysis (PCA), using LRR-data subset:\n\n")
  bigMatSummary(pcaMat,name="pcaMat")
  est.mem <- estimate.memory.for.dat(pcaMat)
  cat(" estimated memory required for",nrow(pcaMat),"x",ncol(pcaMat),"matrix:",round(est.mem,2),
      "GB. If this exceeds available,\n  then expect PCA to take a long time or fail!\n")
  subMat <- as.matrix(pcaMat)
  rm(pcaMat)
  # center using row means
  cat(" centering data by row means...")
  subMat <- subMat - rowMeans(subMat)  #matrix(rep(rowMeans(subMat),times=ncol(subMat)),ncol=ncol(subMat))
  cat(" means for first 10 snps:\n")
  print(round(head(rowMeans(subMat)),10)) # show that centering has worked
  subMat[is.na(subMat)] <- 0 # replace missing with the mean
  cat(" replaced missing data with mean (PCA cannot handle missing data)\n")
  #subMat <- t(subMat) # transpose
  dimz <- dim(subMat)
  if(pcs.to.keep > min(dimz)) { 
    # ensure not trying to extract too many pcs
    warning(paste("selected too many PCs to keep [",pcs.to.keep,"], changing to ",min(dimz),"\n",sep="")) 
    pcs.to.keep <- min(dimz)
  } 
  if(!SVD & (dimz[2]>dimz[1])) {
    cat(" PCA using 'princomp' (faster for datasets with more samples than markers)\n")
    print(system.time(result <- princomp(t(subMat))))
    PCs <- result$scores[,1:pcs.to.keep]
    Evalues <- result$sdev
  } else {
    if(!SVD) {
      cat(" PCA by crossproduct and solving ectors\n")
      cat(" obtaining crossproduct of the matrix and transpose XtX...")
      uu <-(system.time(xtx <- crossprod(subMat)))
      cat("took",round(uu[3]/60,1),"minutes\n")
      cat(" obtaining eigen vectors of the crossproduct XtX...")
      uu <-(system.time(result <- eigen((xtx/nrow(subMat)),symmetric=T)))
      cat("took",round(uu[3]/60,1),"minutes\n")
      PCs <- result$vectors[,1:pcs.to.keep]
      Evalues <- result$values
    } else {
      do.fast <- (!LAP & ((require(irlba) & require(bigalgebra))))
      cat(" PCA by singular value decomposition...") # La.svd gives result with reversed dims. (faster?)
      if(!LAP) {
        if(do.fast) {
          uu <-(system.time(result <- irlba(subMat,nv=pcs.to.keep,nu=0,matmul=matmul))) 
        } else {
          cat("[without 'bigalgebra' package, slow for large datasets]\n")
          uu <-(system.time(result <- svd(subMat,nv=pcs.to.keep,nu=0)))
        }
        cat("took",round(uu[3]/60,1),"minutes\n")
        PCs <- result$v[,1:pcs.to.keep]
        Evalues <- result$d
      } else {
        cat("\n [using LAPACK alternative with La.svd]")
        uu <- (system.time(result<- La.svd(subMat,nv=pcs.to.keep,nu=0)))
        cat("took",round(uu[3]/60,1),"minutes\n")
        PCs <- t(result$vt)[,1:pcs.to.keep]  ##?
        Evalues <- result$d
      }
    }
  }
  rownames(PCs) <- colnames(subMat)
  colnames(PCs) <- paste("PC",1:pcs.to.keep,sep="")
  if(save.pcs) {
    ofn <- dir.fn.cmb(dir$pc,pcs.fn)
    cat(paste("~wrote PC data to file:",ofn,"\n"))
    save(PCs,Evalues,file=ofn) }
  out.dat <- list(PCs,Evalues)
  names(out.dat) <- c("PCs","Evalues")
  return(out.dat)
}


LRR.PCA.correct <- function(pca.result,descr.fn,dir,num.pcs=9,n.cores=1,pref="corrected",
                            big.cor.fn=NULL,write=F,sample.info=NULL,correct.sex=F)
{
  ## using results of a PCA analysis, run correction for 'num.pcs' PCs on a dataset
  # uncorrected matrix
  dir <- validate.dir.for(dir,c("big","pc"))
  origMat <- getBigMat(descr.fn,dir)
  cat("\nRunning Principle Components correction (PC-correction), using LRR-dataset:\n")
  bigMatSummary(origMat,name="origMat")
  if(n.cores>1) { multi <- T } else { multi <- F }
  # get filenames now to add to result later
  rN <- rownames(origMat); cN <- colnames(origMat)
  # run pca.correction using ectors (PCs) and alues from LRR.PCA
  if(!is.list(pca.result)) {
    if(is.character(pca.result)) {
      ofn <- dir.fn.cmb(dir$pc,pca.result) 
      if(file.exists(ofn))
      {
        pca.file <- get(load(ofn))
        cat(" loaded PCA values and vectors\n")
        PCs <- pca.file$PCs
      } else {
        stop("Error: file",ofn,"does not exist\n")
      }
    } else {
      if(ncol(origMat) %in% dim(pca.result))
      {
        #given dimensions = number of samples, assume PCs entered as matrix
        PCs <- pca.result
      } else {
        stop("Error: expecting file name or PC matrix: pca.result\n")
      }
    } 
  } else {
    PCs <- pca.result$PCs
  }
  
  # create new matrix same size, ready for corrected values
  nR <- nrow(origMat); nC <- ncol(origMat)
  cat(" creating new file backed big.matrix to store corrected data...")
  pcCorMat <- filebacked.big.matrix(nR,nC, backingfile=paste(pref,"Bck",sep=""),
                                    backingpath=dir$big, descriptorfile=paste(pref,"Descr",sep=""))
  cat("done\n")
  if(!is.filebacked(pcCorMat) | !is.filebacked(origMat)) {
    warning("at least one of the big.matrices is not filebacked, memory problems may be encountered")
  }
  # in result$vectors, PCs are the columns, only need first 10 or so
  # rows are subjects / samples
  col.sel <- 1:ncol(origMat)
  nPCs <- PCs[,1:num.pcs]; sex.txt <- ""
  if(correct.sex) { 
    ## add a column to the PC matrix which will allow to covary for sex too
    coln <- which(tolower(colnames(sample.info)) %in% c("gender","sex")) 
    if(length(coln)>0) { 
      indx <- match(colnames(origMat),rownames(sample.info))
      if(all(!is.na(indx))) {
        sex.vec <- sample.info[indx,coln[1]] 
        sex.vec[is.na(sex.vec)] <- mean(sex.vec,na.rm=T) # replace missing sex with mean 
        nPCs <- cbind(sex.vec,nPCs) ; sex.txt <- " (covarying for sex)"
      } else { warning("could not correct for sex as sample.info was missing some IDs") }
    }
  }
  nPCs <- cbind(rep(1,nrow(nPCs)),nPCs) # this adds intercept term for lm.fit() [remove if using lm() ]
  cat(" correcting by principle components",sex.txt,", taking the LRR lm-residual for each SNP\n",sep="")
  jj <- proc.time()
  nsamples <- ncol(origMat)
  num.snps <- nrow(origMat); sampz <- 1:nsamples
  snps.per.proc <- 400;   flush.freq <- 20
  if(nsamples>10000) { snps.per.proc <- 300 }; if(nsamples>20000) { snps.per.proc <- 150 }
  if(nsamples>40000) { snps.per.proc <- 50 }; if(nsamples>80000) { snps.per.proc <- 10 }
  ## assume if we have lots of cores, we'd also have lots of RAM too
  if(n.cores>5) { snps.per.proc <- snps.per.proc*2; flush.freq <- flush.freq*2 } 
  if(n.cores>15) { snps.per.proc <- snps.per.proc*2 }
  snps.per.proc <- max(snps.per.proc,n.cores) # at least 1 snp per core as a minimum
  stepz <- round(seq(from=1,to=num.snps+1,by=snps.per.proc))
  if((tail(stepz,1)) != num.snps+1) { stepz <- c(stepz,num.snps+1) }
  split.to <- length(stepz)-1
  big.extras <- T # flush memory every 'n' iterations.
  
  # this simple way works (instead of big for-loop) but hogs memory and is no faster
  # [NB: requires transpose of target corrected big matrix dimensions]
  ### pcCorMat <- apply(origMat,1,PC.fn,nPCs=nPCs,col.sel=sampz)
  for (dd in 1:split.to)
  {
    x1 <- stepz[dd]; x2 <- stepz[dd+1]-1 #subset row selection
    # use of this 'sub.big.matrix' structure, stops the memory leak behaviour which spirals
    # the memory relating to 'origMat' out of control. 
    next.rows <- sub.big.matrix(origMat, firstRow=x1, lastRow=x2, backingpath=dir$big )
    # next.rows is now a pointer to a matrix subset, must use 'as.matrix' to coerce to a regular R object 
    if(multi) {
      pcCorMat[x1:x2,] <- PC.fn.mat.multi(as.matrix(next.rows),nPCs,mc.cores=n.cores)
    } else {
      pcCorMat[x1:x2,] <- PC.fn.mat.apply(as.matrix(next.rows),nPCs)
    }
    loop.tracker(dd,split.to)
    ## Every 'flush.freq' iterations clean up the memory, remove the 
    ##  big.matrix object 'pcCorMat' and re-attach it 
    if(dd %% flush.freq == 0) {    
      fl.suc <- flush(pcCorMat) & flush(next.rows)
      if(!fl.suc) { cat("flush failed\n") } 
      gc()  # garbage collection
      if(big.extras) {
        RR <- describe(pcCorMat)
        rm(pcCorMat)
        pcCorMat <- attach.big.matrix(RR,path=dir$big)
      }
    }
    rm(next.rows) # remove the sub-matrix pointer each iteration or this memory builds up 
  }
  
  options(bigmemory.allow.dimnames=TRUE)
  rownames(pcCorMat) <- rN;  colnames(pcCorMat) <- cN 
  ll <- proc.time()
  cat(paste(" LRR PC-Correction took",round((ll-jj)[3]/3600,3),"hours\n"))
  flush(pcCorMat) # should allow names to take  
  cat("\nPC-corrected dataset produced:\n")
  bigMatSummary(pcCorMat,name="pcCorMat")
  
  mat.ref <- describe(pcCorMat)
  if(write) {
    if(is.null(big.cor.fn) | !is.character(big.cor.fn)) {
      big.fn <- paste("describePCcorrect",num.pcs,".RData",sep="")
    } else {
      big.fn <- big.cor.fn[1]
    }
    ofn <- dir.fn.cmb(dir$big,big.fn)
    save(mat.ref,file=ofn)
    cat(paste("~wrote PC-corrected data description file to:\n ",ofn,"\n"))
    return(big.fn)
  } else {
    return(mat.ref)
  }
}


PC.fn.mat <- function(next.rows,nPCs)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time)
  col.sel <- 1:ncol(next.rows)
  for (dd in 1:nrow(next.rows)) {
    # compiled PC.fn should speed up these ops a little
    next.rows[dd,] <- PC.fn(next.rows[dd,],nPCs,col.sel) 
  }  
  return(next.rows)
}


PC.fn.mat.apply <- function(nextrows,nPCs)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time), vectorized version
  # testing shows the for-loop (non-vectorized) to be slightly faster, maybe because of t()
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  col.sel <- 1:ncol(nextrows)
  nextrows <- t(apply(nextrows,1,PC.fn.2,nPCs=nPCs,col.sel=col.sel))
  return(nextrows)
}


PC.fn.mat.multi <- function(nextrows,nPCs,mc.cores=1)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time), vectorized version
  # testing shows the for-loop (non-vectorized) to be slightly faster, maybe because of t()
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  col.sel <- 1:ncol(nextrows)
  nextrows <- lapply(seq_len(nrow(nextrows)), function(i) nextrows[i,]) # multi slows this down
  #nextrows <- multicore::mclapply(nextrows,PC.fn,nPCs=nPCs,col.sel=col.sel,mc.cores=mc.cores)
  nextrows <- multicore::mclapply(nextrows,PC.fn.2,nPCs=nPCs,col.sel=col.sel,mc.cores=mc.cores)
  nextrows <- do.call("rbind",nextrows)
  return(nextrows)
}


PC.fn <- function(next.row,nPCs,col.sel)
{
  # apply PC correction for a single SNP, allowing for missing data.
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals
  return(next.row)
}


PC.fn.2 <- function(next.row,nPCs,col.sel)
{
  # apply PC correction for a single SNP, allowing for missing data.
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  next.row[sel] <- lm.fit(x=nPCs[sel,],y=next.row[sel])$residuals
  return(next.row)
}

