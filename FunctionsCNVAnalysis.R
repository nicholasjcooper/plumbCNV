options(stringsAsFactors=FALSE)


########### DEFINE INTERNAL FUNCTIONS ############

auto.mode <- function(on=T,set=F) {
  if(set) {
    options(plumb.cnv.auto.mode=on)
  }
  auto.gs <- getOption("plumb.cnv.auto.mode"); if(is.null(auto.gs)) { auto.gs <- F }
  return(auto.gs)
}



clean.fn <- function(x,fail=NA,fn=function(x) { x }) {
  # allows an sapply style function to only work on valid values
  if(!is.null(x)) { 
    x <- x[!is.na(x)]; x <- x[(x!="")]; return(fn(x)) 
  } else {  return(fail) } 
}


penn.name.files <- function(dir,write.files=T,ret.fn=F,relative=F) 
{
  ## PennCNV commands are more easily run using files listing the locations of each 
  # sample file. This function generates a set of such files using relative or
  # absolute paths, given a specified number of directories
  ndir <- count.penn.dirs(dir)
  dir.ch <- .Platform$file.sep
  subj.in.each <- lapply(paste(dir$cnv.raw,"p",1:ndir,sep=""),list.files)
  fnms <- cat.path(dir$cnv.raw,paste("p",1:ndir,"names",sep=""),ext="txt")
  my.fnms <- vector("list",length(subj.in.each))
  if(ndir==length(subj.in.each)) {
    for (cc in 1:ndir) {
      rel.fnms <- paste(paste("p",cc,dir.ch,sep=""),subj.in.each[[cc]],sep="")
      abs.fnms <- paste(dir$cnv.raw,rel.fnms,sep="")
      cnv.rel.fnms <- gsub(dir$cnv,paste(".",dir.ch,sep=""),abs.fnms)
      if(relative) { my.fnms[[cc]] <- cnv.rel.fnms } else { my.fnms[[cc]] <- abs.fnms }
      if(write.files) { writeLines(my.fnms[[cc]],con=fnms[cc]) }
    }
    if(write.files) { cat("\n~wrote",ndir,"PennCNV sample selection files to:\n ",dir$cnv.raw,"\n") }
    names(subj.in.each) <- basename(fnms)
  } else {
    stop(paste("Error:",ndir,"PennCNV directories spawned but",length(subj.in.each),"subjects lists"))
  }
  if(ret.fn) {
    return(fnms)
  } else {
    return(my.fnms)
  }
}


import.called.cnvs <- function(dir,out.format=c("Ranges","cnvGSA"),
                               cnv.result.cats=c("allCNV","allDel","allDup","rareDEL","rareDUP"),
                               pheno.file="pheno.lookup.txt",add.phenotype=T)
{
  dir <- validate.dir.for(dir,"cnv.plk")
  num.cats <- length(cnv.result.cats)
  out.format <- out.format[1]
  cnvResults <- vector("list",num.cats)
  for(cc in 1:num.cats) 
  {
    fn <- cat.path(dir$cnv.qc,cnv.result.cats[cc],ext="cnv")
    cnvResults[[cc]] <- switch(out.format,cnvGSA=plink.to.cnvGSA(fn),
                               Ranges=plink.to.Ranges(fn))
  }
  if(add.phenotype) {
    if((out.format=="Ranges") & is.file(pheno.file,dir$ano,dir)) {
      for (dd in 1:num.cats) {
        cnvResults[[dd]] <- add.info.to.ranges(cnv.ranges=cnvResults[[dd]],dir=dir,
                            target="phenotype",info.file=find.file(pheno.file,dir$ano,dir)) 
      }
    } else {
      print(is.file(pheno.file,dir$ano,dir))
    }
  } 
  names(cnvResults) <- cnv.result.cats
  return(cnvResults)
}


get.baf.markers <- function(baf.des.fn, samp.list, dir, low.ram=T) {
  # attach datafile bigmemory object
  dir <- validate.dir.for(dir,c("big"),warn=F)
  bigBAF <- getBigMat(baf.des.fn,dir$big)
  cat("\nRegenerating BAF average data from BAF matrix\n")
  to.keep <- which(colnames(bigBAF) %in% samp.list)
  #mean.sel <- function(X,to.keep) { mean(X[to.keep],na.rm=T) }
  #marker.means <- apply(bigMat2,1,mean.sel,to.keep)
  nC <- length(to.keep); nR <- nrow(bigBAF) ; split.to <- (((nC*nR) %/% (6*10^7))+1)
  marker.means <- numeric(nR); error.times <- NULL
  stepz <- round(seq(from=1,to=nR+1,length.out=round((split.to+1))))
  if((tail(stepz,1)) != nR+1) { stepz <- c(stepz,nR+1) }
  split.to <- length(stepz)-1
  for (cc in 1:split.to)
  {
    x1 <- stepz[cc]; x2 <- ((stepz[cc+1])-1)
    if((x2-x1)>=0) {
      marker.means[x1:x2] <- rowMeans(bigBAF[x1:x2,to.keep],na.rm=T)
    } else {
      error.times <- c(error.times,x1,x2)
    }
    loop.tracker(cc,split.to)
    if(cc %% 10 == 0) {
      # reset memory after every 10 chunks to prevent unnecessary increasing RAM use #
      if(low.ram) {
        RR <- describe(bigBAF); rm(bigBAF); gc(); bigBAF <- attach.big.matrix(RR,path=dir$big)
      }
    }
  }
  odds <- function(x,T) { x[(1:length(x)) %% 2 != 0] }; evens <- function(x,T) { x[(1:length(x)) %% 2 == 0] }
  if(!is.null(error.times)) { warning(paste("found",length(error.times)/2,"illegal row ranges:",
                                      paste(paste(odds(error.times),evens(error.times),sep="-"),collapse=","))) }
  names(marker.means) <- rownames(bigBAF)
  return(marker.means)
}


do.scree.plots <- function(eigenv,dir="",fname="ScreePlotPCA.pdf",elbow=9,n.comp=30,printvar=T,writefig=T,nsamp=NA,...) {
  # plumbCNV wrapper to do SCREE PLOTS AND calculate EIGENVALUE VARIANCE after a PCA
  cat(" generating scree plots for principal components analyses\n")
  dir <- validate.dir.for(dir,"pc")
  if(writefig) {  ofn <- cat.path(dir$pc,fname);   pdf(ofn) }
  out <- pca.scree.plot(eigenv=eigenv,elbow=elbow,n.comp=n.comp,
                         printvar=printvar,nsamp=nsamp,...)
  if(writefig) { dev.off() ;  cat(paste("~wrote file:",ofn,"\n")) }
  return(out)
}


run.PCA.correct <- function(DT=NULL,dir=NULL,pc.to.keep=.13,assoc=F,num.pcs=9,n.store=30,correct.sex=F,
                             comparison=T,big.lrr.fn="combinedBigMat.RData",n.cores=1,
                             snp.sub.fn="pca.snp.subset.txt",restore.mode=F,pref.suf="",
                             sub.pref="PCAMatrix",post.pref="PostPC",cor.pref="corrected",
                             sub.fn="pcaSubMat.RData", SVD=T,LAP=F,big.cor.fn=NULL,add.int=F,
                             pcs.fn="PCsEVsFromPCA.RData", med.fn="postpcmed.RData",
                             comp.gc=F,comp.dlrs=T,comps=c("plate","phenotype"),
                             snp.info=NULL,sample.info=NULL,raw.stat.fn=NULL,ucsc="hg18")
{
  # Run PCA on a subset of the combined dataset
  # then correct the combined dataset according to the first 'n' PC-components 
  # then stats on the resulting distributions can be calculated
  # write results to big.matrix files and pdfs, etc, then
  # returns raw, sample-QC and PC-corrected sets of stats for comparison
  if(n.cores>1) { multi <- T; must.use.package("multicore") } else { multi <- F }
  load.all.libs(big=c("bigmemory","biganalytics","bigalgebra"),more.bio="irlba") # load all main plumbcnv libraries
  if(is.data.tracker(DT)) {
    DTs.required <- c("ngrps","sample.info","snp.info","big.qc","stats")
    is.valid <- req.valid(DT,DTs.required); print(is.valid)
    if(!is.valid) { 
      ## if invalid let it try to work via loading in separate function parameters
      warning("file required for 'run.PCA.correct' invalid in data.tracker")
    } else {
      if(is.null(dir)) {
        dir <- getSlot(DT,"dir")
        #print(dir)
      }
      big.lrr.fn <- getSlot(DT,"big.qc")
      raw.stat.fn2 <- getSlot(DT,"stats")
      if(is.list(raw.stat.fn2) & is.ch(raw.stat.fn2)) 
      { raw.stat.fn2 <- unlist(raw.stat.fn2) } #cat("..unlisted!\n") }
      if(is.character(raw.stat.fn2)) { raw.stat.fn <- raw.stat.fn2 } # otherwise default to main fn arg
    }
  } else { n.grps <- NA }
# print(dir[[1]]) ; 
  dir <- validate.dir.for(dir,c("big","ano","gc"))
  #### if number of PCs to correct for is zero, then PC and correction should be skipped
  # optionally add a suffix to all 'pref's 
  last.two <- T
  if(is.character(pref.suf)) {
    ## allows setting of different suffixes for each output file, e.g, if running parallel versions, etc
    # pref for the PCA subset bigmat, pref for the post pc stats comparison, pref for the pca corrected bigmat,
    # pref for the ector results, pref for gc medians
    if(length(pref.suf)==3) { j <- c(1:3) ; last.two=F } else { 
      if(length(pref.suf)==5) { j <- 1:5 } else { j <- rep(1,5) }
    }
    sub.pref=paste(sub.pref,pref.suf[j[1]],sep="")
    post.pref=paste(post.pref,pref.suf[j[2]],sep="")
    cor.pref=paste(cor.pref,pref.suf[j[3]],sep="")
    if(last.two) {
      pcs.fn=cat.path(fn=pcs.fn,suf=pref.suf[j[4]])
      med.fn=cat.path(fn=med.fn,suf=pref.suf[j[5]])
    }
  }
  if(num.pcs > n.store) { n.store <- num.pcs } # make sure we are extracting enough pcs
  if(n.store>300) { warning(paste("high values of n.store [",n.store,"] cause the PCA to take a long time\n",sep="")) }
  if(pc.to.keep>1) { pc.to.keep <- pc.to.keep/100 } # in case it's in percentage, not fraction
  # load snp and sample.data if not valid from function args
  if(!is.data.frame(sample.info)) { cat("sample.info: "); sample.info <- read.sample.info(dir,nprev=3) }
  if(is(snp.info)[1]!="RangedData") { 
    cat("snp.info: "); snp.info <- read.snp.info(dir,nprev=3) 
    if(!is.data.frame(sample.info)) { stop("couldn't find snp.info object in dir$ano") }
  }
  uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  # create subset datafile
  if(num.pcs<1) {
    corrected.ref <- big.lrr.fn
    pcs.fn <- NULL
  } else {
    subset.descr <- get.PCA.subset(dir=dir,pc.to.keep=pc.to.keep,assoc=assoc,big.fn=big.lrr.fn,
                                   snp.sub.fn=snp.sub.fn,use.current=restore.mode,pref=sub.pref,n.cores=n.cores,
                                   descr.fn=sub.fn,nprev=3,snp.info=snp.info,sample.info=sample.info)
  
    # run the PCA  ### here was 'sub.fn' - what do we do with that now???
    cat("\nRunning Principle Components Analysis (PCA), using LRR-data subset:\n\n")
    pca.result <- big.PCA(bigMat=subset.descr,dir=dir,pcs.to.keep=n.store,
                          SVD=SVD,LAP=LAP,save.pcs=T,pcs.fn=pcs.fn,verbose=T) 
    nsamp <- length(which(sample.info$QCfail==0))
    varz <- do.scree.plots(pca.result$Evalues,dir=dir,elbow=num.pcs,nsamp=nsamp)
    cat("",round(cumsum(varz)[num.pcs]*100,1),"% Est. LRR variance explained by first",num.pcs,"components.\n")
    # correct using the PC-components
    corrected.ref <- LRR.PCA.correct(pca.result,big.lrr.fn,dir=dir,num.pcs=num.pcs,n.cores=n.cores,
                                                pref=cor.pref,big.cor.fn=big.cor.fn,write=T,
                                                sample.info=sample.info,correct.sex=correct.sex, add.int=add.int)
  }
  corrected.bigMat <- getBigMat(corrected.ref,dir) #; rm(corrected.ref) #(in case it was a bigmat)
  # optionally make comparison of pre and post-pc distributions
  if(comparison) {
    med.fn <- cat.path("",med.fn,suf=paste(num.pcs),ext="RData")
    outlist <- post.pc.comparison(corrected.bigMat,dir=dir,snp.info=snp.info,sample.info=sample.info,
                             pref=post.pref,med.chunk.fn=med.fn,ucsc=ucsc,n.cores=n.cores,
                             dlrs=comp.dlrs,gc=comp.gc,comps=comps,raw.stat.fn=raw.stat.fn)
    stat.fn <- cat.path(dir$qc.lrr,outlist$stat.list$tab.file)
  } else { 
    med.fn <- NULL ; stat.fn <- NULL
  }
  #else {
    #return(describe(corrected.bigMat))
  #}
  if(is.data.tracker(DT)) {
    DT <- setSlot(DT,pca=corrected.ref,median=med.fn,eigen=pcs.fn,
                              stat=stat.fn,n.pcs=num.pcs,proc.done=4,warns=F,settings=settings)
    return(DT)
  } else {
    return(corrected.ref)
  }
}



#plumb cnv internal function
load.data.to.bigmat <- function(dat.fn,inputType="TXT",bck,des,dir,Hopt=F,RNopt=T)
{
  # get data stored in a text or binary file and get it into a bigmatrix object format
  dir <- validate.dir.for(dir,c("ano","big"),warn=F)
  if(inputType=="RDATA")
  {
    # attach datafile bigmemory object
    cat("\nLoading RData file...")
    this.mat <- load(dat.fn)
    if (this.mat!="m.by.s.matrix") { 
      m.by.s.matrix <- get(this.mat)
      rm(this.mat) 
    }
    cat("complete\n")
    if(is.null(colnames(m.by.s.matrix)))
    {
      snp.fn <- cat.path(dir$ano,"snpNames.txt")
      snp.list <- readLines(snp.fn)
      if (ncol(m.by.s.matrix)==length(snp.list))
      {
        colnames(m.by.s.matrix) <- snp.list
        cat(paste("*warning: added column names to matrix from",snp.fn,"\n"))
        cat("if replacement with these names is undesired please check/modify the .R source code\n")
      } else {
        cat("Error: matrix has no column names and default file doesn't\n")
        cat("match number of columns, please check/modify the .R source code\n")
        stop()
      }
    }
    cat(" saving datafile as big matrix\n")
    bigMat <- as.big.matrix(m.by.s.matrix, backingfile=bck,
                            backingpath=dir$big, descriptorfile=des)
  } else {
    # assume inputType=="TXT"
    cat("\nLoading TAB file...(this will probably be slow)...")
    read.big.matrix(dat.fn, sep = '\t', header = Hopt,
                    has.row.names=RNopt, ignore.row.names=FALSE,
                    backingfile = bck, backingpath = dir$big,
                    descriptorfile = des, extraCols = NULL)
    cat("complete\n")
  }
}


post.pc.comparison <- function(corrected.bigMat, dir, gc=F, dlrs=T, pref="", snp.info, sample.info=NULL, nSD=3,
             n.cores=1,med.chunk.fn="postpcmed.RData", comps=c("plate"), raw.stat.fn=NULL, ucsc="hg18",restore.mode=F) 
{
  # compare pre-pc and post-pc stats; might want to pass sample.info in for consistency
  cat("\nRunning analyses to compare pre/post PC distributions (can take a long time)\n")
  if(!is(snp.info)[1]=="RangedData") { stop("snp.info was not a valid RangedData object") }
  outlist <- lrr.sample.dists(corrected.bigMat,snp.info=snp.info,dir=dir,n.cores=n.cores,nSD=nSD,
                          gc=gc,dlrs=dlrs,pref=pref,med.chunk.fn=med.chunk.fn,restore.mode=restore.mode,ucsc=ucsc,
                          write.excl=F,extended.report=F,sum.fn="PostPCSampleQCTables")
  cleaned.stat.table <- outlist$stat.table
  # could here compare with the pre-cleaned one
  if(!is.data.frame(sample.info)) { sample.info <- read.sample.info(dir) }
  sample.info <- validate.samp.info(sample.info,dir=dir,QC.update=T,verbose=F,proc=4) #notDEP
  ## optionally creates 3 x 3 plots of raw,sample.qc,pc.corrected data, returns the 3 stat tables
  twc <- three.way.comparison(cleaned.stat.table,sample.info,dir=dir,
                              batch.comps=comps,height=10,width=10,raw.stat.fn=raw.stat.fn) 
  out <- list(twc,outlist); names(out) <- c("comparison.3.ways","stat.list")
  return(out)
}


lrr.sample.dists <- function(bigMat,snp.info,dir,gc=F,dlrs=T,pref="",med.chunk.fn="",restore.mode=F,
                              plotfn="DistributionsBoxPlot.pdf",tabfn="StatsPerSample.tab",
                             nSD=3,lo.cut=c("LB",NA,"LB"),hi.cut=c("UB","UB","UB"),ucsc="hg18",
                             write.excl=T,extended.report=T,sum.fn="SampleQCTables",n.cores=1,...) 
{
  # calculates sample distributions for mean, optionally DLRS/StDev, GC-wave.
  # passes arguments to calc.chr.info via get.chr.info.filt and calculate.gc.for.sample
  # better off passing a complete snp.info object, as will filter out any not in the data automatically
  # lo.cut=c("LB",NA,"LB"),hi.cut=c("UB","UB","UB") are which cutoff to use for the columns of stat.table
  # which are generally Mean, DLRS and GCWave respectively, although if DLRS or GCwave turned off, then this
  # can change
  load.all.libs() # load all main plumbcnv libraries
  dir <- validate.dir.for(dir,c("qc.lrr","big","ano"))

  #print(dim(bigMat)); print(dim(snp.info)); print(length(which(!rownames(bigMat) %in% rownames(snp.info))))
  if(!is(snp.info)[1]=="RangedData") { stop("snp.info was not a valid RangedData object") }
  snp.info <- snp.info[(rownames(snp.info) %in% rownames(bigMat)),]  # snp.info object matching the bigmatrix rows exactly
    
  stat.mat <- calc.main.stats.on.big(bigMat,dir=dir,deleteDLRS=T,snp.info=snp.info,
                              skip.dlrs=!dlrs,apply.method=T,n.cores=n.cores)
  #debug(calculate.gc.for.samples)
  if(gc) {
    cat("\nAdding GC wave to LRR stats, be aware this can be slow\n")
    gc.wave <- calculate.gc.for.samples(bigMat,snp.info=snp.info,dir=dir,med.chunk.fn=med.chunk.fn,
          restore.mode=restore.mode,ucsc=ucsc,n.cores=n.cores)
    ## combine GC and main, then report
    c.nms <- rownames(stat.mat)[rownames(stat.mat) %in% row.names(gc.wave)]
    if(!dlrs) { stat.set <- c("Mean","StDev") } else { stat.set <- c("Mean","DLRS") }
    #c("r_GC","S_WF","S_GCWF"), so col 2=S_WF below = wave factor (not just GC wave)
    stat.table <- as.data.frame(cbind(stat.mat[c.nms,stat.set],gc.wave[c.nms,2]))
    colnames(stat.table)[ncol(stat.table)] <- "GCWave"
  } else {
    cat(" calculating LRR stats without GC wave (faster)\n")
    stat.table <- stat.mat
  }
  if(!is.na(plotfn)) {
    ofn <- cat.path(dir$qc.lrr,plotfn,pref=pref)
    LRR.boxplot(stat.table,ofn)
  }
  if(!is.na(tabfn)) {
    ofn <- cat.path(dir$qc.lrr,tabfn,pref=pref)
    write.table(round(stat.table,5),sep="\t",col.names=T,row.names=T,quote=F,file=ofn)
    cat("~wrote table\n",ofn,"\n")
  } else {
    warning("sample-wise statistics table not written to file, no file name given")
  }
  s.tab <- lrr.stats.tab(stat.table,nSD=nSD)

  #cat("\nCohort distribution statistics for LRR-stats:\n"); print(s.tab,digits=4); cat("\n")
  sum.fn <- file(cat.path(dir$qc.lrr,fn=sum.fn,suf=pref,ext="txt"),open="w")
  dual.cat(file=sum.fn,"\n\nCohort distribution statistics for LRR-stats:\n")
  print.large.to.file(s.tab,con=sum.fn,to.scrn=T,digits=4)  #print((res.tab))
  close(sum.fn)
  cat("\n")
  
  if(extended.report) {
    rez <- get.excl.list.from.stats(s.tab,stat.table,dir=dir,lo.stat.nm=lo.cut,hi.stat.nm=hi.cut,
                                    nSD=nSD,writeFiles=write.excl,append=T) # append!
    # now compile venn list and do failing example plots
    do.venn.QC(rez,dir,pref) # DRAW VENN DIAGRAM
    plot.extreme.samples(rez=rez, stat.table=stat.table, bigMat2=bigMat,
                         snp.info=snp.info, dir=dir, pref=pref, ucsc=ucsc)
    lrr.boundary.scatter(stat.table,s.tab,dir=dir,fn.pre=pref,fn.suf="",col="darkblue")
  }
  outlist <- list(stat.table,s.tab,ofn)
  names(outlist) <- c("stat.table","stat.summary","tab.file")
  return(outlist)
}


get.excl.filter.for <- function(object,dir)
{
  # create sample-wise filter for any object, excluding samples in exclusion files directory
  bad.samps <- paste(get.all.samp.fails(dir))
  if (is.list(object))
  {
    for(cc in 1:length(object)) {  object[[cc]] <- get.excl.filter.for(object[[cc]],dir) }
  } else {
    if(is.character(object))
    {
      object <- !(object %in% bad.samps)
    } else {
      if(!is.null(colnames(object))) {
        object <- !(colnames(object) %in% bad.samps)
      } else {
        if(!is.null(rownames(object))) {
          object <- !(rownames(object) %in% bad.samps)
        } else {
          if(is.factor(object)){
            object <- paste(object)
            object <- !(object %in% bad.samps)
          } else {
            cat(" object",is(object)[1],"seemed unsuitable for exclusion of ids, was left unchanged.\n")
          }
        }
      }
    }
  }
  return(object)
}
  
  
exclude.combine.big.mats <- function(bigList,dir,pref=paste("SampQC_Combined",sep=""),debug=F) {
  ## apply exclusions from sample QC and combine all matrices from different 
  ## datasets into a single big.matrix object
  dir <- validate.dir.for(dir,c("big"),warn=F)
  pulp <- function(vec,ch="-") { paste(vec,collapse=ch) }
  if(length(bigList)==1)
  {
    bigUnified <- getBigMat(big.exclude.sort(bigList[[1]], dir=dir, pref=pref),dir$big,verbose=F)
    descr <- describe(bigUnified); rm(bigUnified)
    return(descr)
  } else {
    # obtain all dimensions, names, etc
    multi.file.list <- load.big.pointer.list(bigList,dir)
    bigMatLst <- multi.file.list$bigMatLst  
    all.samps <- multi.file.list$combined.samples
    samp.list <- multi.file.list$samples
    smp.szs <- multi.file.list$grp.sizes
    n.data.files <- length(bigMatLst)
    tf <- T
    if(n.data.files>1) {
      # check whether each file in turn has the same rownames
      for (jj in 1:(n.data.files-1)) { tf <- (tf & all(rownames(bigMatLst[[jj]])==rownames(bigMatLst[[jj+1]]))) }
    }
    if(tf) {
      snp.labels <- rownames(bigMatLst[[1]])          } else {
      stop("Error: Rownames (snp labels) of each big.matrix file were different") 
    }
    ## need mask filters for each input sub-dataset
    ## T/F masks remove bad samples from each subset
    f.SL <- get.excl.filter.for(samp.list,dir) # separate filters per subgroup
    f.AS <- get.excl.filter.for(all.samps,dir) # all samples filter
    smp.szsF <- sapply(f.SL,function(X) { length(which(X)) })
    cat(" original sample sizes in separate matrices:",paste(smp.szs,collapse=","),"\n")
    cat(" new sample sizes in combined matrix:",paste(smp.szsF,collapse=","),"\n")
    ## need to adjust final matrix size accordingly and all refs
    
    if(all((diff(as.numeric(sapply(bigMatLst,nrow))))==0))
    {
      des <- paste(pref,"descrFile",sep="_")
      bck <- paste(pref,"bckFile",sep="_")
      nR <- length(snp.labels)
      nC <- length(all.samps);  nCf <- length(all.samps[f.AS])
      cat(" combing matrices, expect this to take some time\n")
      bigUnified <- big.matrix(nR,nCf, backingfile=bck,
                               backingpath=dir$big, descriptorfile=des)
      d2 <- d1 <- 0
      for (dd in 1:n.data.files)
      {
        split.to <- max(1,(smp.szs[dd]) %/% 300) # save RAM without creating groups too small to process
        stepz <- round(seq(from=1,to=smp.szs[dd]+1,length.out=round((split.to+1))))
        if((tail(stepz,1)) != smp.szs[dd]+1) { stepz <- c(stepz,smp.szs[dd]+1) }
        split.to <- length(stepz)-1
        if (debug) { print(stepz) ; print(dim(bigUnified)) }
        #^ check this in FN!
        cat(" copying file",dd,"of",n.data.files,"to combined big.matrix object:\n")
        for (cc in 1:split.to)
        {
          # within submatrix cols
          c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
          # within combined matrix cols
          #d1 <- c1+c(0,(cumsum(smp.szs)))[dd]; d2 <- d1 + (c2-c1)
          # do the copying
          lilColRange <- c(c1:c2)[ (f.SL[[dd]])[c(c1:c2)] ] # filter samps in the old sub-matrix
          nvalid <- length(which((f.SL[[dd]])[c(c1:c2)]))
          d1 <- d2+1; d2 <- (d1 + nvalid - 1)
          bigColRange <- c(d1:d2) # filter samps in the new big matrix
          if(debug) {
            cat(" comb'd range:",pulp(range(bigColRange)),"length:",diff(range(bigColRange))+1,"\n")
            cat(" subset range:",pulp(range(lilColRange)),"length:",diff(range(lilColRange))+1,"valid:",nvalid,"\n\n")
          } else {
            loop.tracker(cc,split.to)
          }
          if(is.finite(sum(lilColRange))) {
            bigUnified[1:nR,bigColRange] <- bigMatLst[[dd]][1:nR,lilColRange]
          } else {
            cat(" Warning: empty interval ignored\n")
          }
        }
      }
      cat(" combining complete, converting result to big matrix\n")
    } else {
      stop("Error: number of rows (snps) for each group matrix don't match\n")
    }
    options(bigmemory.allow.dimnames=TRUE)
    colnames(bigUnified) <- all.samps[f.AS]
    cat(" added colnames\n")
    rownames(bigUnified) <- snp.labels  
    cat(" added rownames\n")
    flush(bigUnified) # hopefully this will ensure the row/colnames are added to the file backing
    descr <- describe(bigUnified) 
    return(descr)
  }
}


load.big.pointer.list <- function(descr.list,dir)
{
  ## taking a list of big matrix descriptors (of some sort) return
  # list of pointers to each big.matrix, plus information about the comprised samples
  n.data.files <- length(descr.list)
  bigMatLst <- samps <- vector("list", n.data.files)
  samp.list <- NULL
  smp.szs <- numeric()
  for (cc in 1:n.data.files)
  { 
    bigMatLst[[cc]] <- getBigMat(descr.list[[cc]],dir) 
    samps[[cc]] <- colnames(bigMatLst[[cc]])
    samp.list <- c(samp.list,samps[[cc]])
    smp.szs[cc] <- length(samps[[cc]])
  }
  outlist <- list(bigMatLst,samp.list,smp.szs,samps)
  names(outlist) <- c("bigMatLst","combined.samples","grp.sizes","samples")
  return(outlist)
}


do.quick.LRR.snp.assocs <- function(descr.fn,sample.info=NULL,use.col="phenotype",dir,p.values=T,n.cores=1)
{
  ## use sample.info
  # create list of bigMatrix locations
  # go 'N' snps at time, concatenate into 1 file, run regression assocs
  # for 2 phenotypes/grps gives t values - ordinally equivalent to logistic regression with 2 groups
  if(!all(c(use.col,"QCfail") %in% colnames(sample.info))) { stop(paste("sample.info was invalid for association tests, need",use.col,"and QCfail columns")) }
  cat(" running SNP-wise LRR-intensity tests against",use.col,"to filter most associated\n")
  bigMat <- getBigMat(descr.fn,dir)
  samp.list <- colnames(bigMat)
  tot.samps <- length(samp.list)
  # check that sample.info is valid, if not attempt to fix it
  if(!is.data.frame(sample.info)) { sample.info <- read.sample.info(dir) }
  ## ENSURE ONLY USING SAMPLES IN THE BIGMATRIX ##
  sample.info <- validate.samp.info(sample.info,dir=dir,QC.update=T,verbose=F)
  sample.info <- sample.info[rownames(sample.info) %in% samp.list,]
  if(length(which(sample.info$QCfail!=0))>0) { warning("bigMat still contains QC failing samples") }
  if(ncol(bigMat)!=nrow(sample.info)) { 
    n.mis <- length(which(colnames(bigMat) %in% rownames(sample.info)))
    warning(n.mis,"samples in BigMat were not in sample.info [failure likely]")
  }
  ## determine test to use based on number of phenotypes ##
  n.phenos <- length(table(sample.info[[use.col]],useNA=NULL))
  t.type <- "single"
  if(n.phenos==2) { t.type <- "t.test"}
  if(n.phenos>2) { t.type <- "anova"}
  cat(" found ",n.phenos," ",use.col,"s, ",t.type," will be used to summarise most associated SNPs for LRR vs ",use.col,"\n",sep="")
  three.test <- function(col,pheno) { return(summary(aov(col~pheno))[[1]][["F value"]][1]) }
  two.test <- function(col,pheno) { return((cor.test(col,pheno)$statistic)^2)  }
  ph.test <- switch(t.type,anova=three.test,t.test=two.test,single=NULL)
  if(is.null(ph.test)) { stop("Error: used option for association test by",use.col,"but there is only 1 type in file")}
  
  snp.labels <- rownames(bigMat)
  full.size <- length(snp.labels)
  good.mem.lim <- 10^6
  opt.size <- round(good.mem.lim/tot.samps)
  n.segs <- ceiling(full.size/opt.size)
  seg.starts <- 1+c(0:(n.segs+2))*opt.size
  seg.starts[seg.starts>full.size] <- full.size+1 # +1 ensures final snp included
  seg.starts <- seg.starts[!duplicated(seg.starts)]
  results <- vector("list", n.segs)
  pheno <- sample.info[[use.col]]
  test.seg <- matrix(numeric(),nrow=opt.size,ncol=tot.samps)
  
  # break analysis into chunks that fit in memory
  # NB: avoided parallel computing - it does not seem to add any speed here
  kk <- proc.time()[3]
  if(n.cores>1 & require(multicore)) {
    Fvalues <- bmcapply(bigMat,1,FUN=ph.test,dir=dir$big,by=200,n.cores=n.cores,pheno=pheno)
  } else {  
    for (dd in 1:n.segs)
    {
      row.subset <- (seg.starts[dd]):(seg.starts[dd+1]-1)
      nr <- length(row.subset)
      # test subset of snps for each segment in turn
      test.seg[1:nr,] <- bigMat[row.subset,] 
      results[[dd]] <- apply(test.seg[1:nr,],1,ph.test,pheno=pheno)
      loop.tracker(dd,n.segs)
    }
    Fvalues <- (do.call(c,results))
  }
  #Fvalues <- round(Fvalues,4) # remember t^2 = F
  cat(" took",round((proc.time()[3]-kk)/60,1),"minutes\n")
  ### option to return p values? could also have option for CHI squared.
  ## also should be able to the bigmatrix apply for a speedup
 # p.from.t <- function(t,df) { pt(t, df, lower=FALSE) }
  p.from.f <- function(FF,k,n) {  pf(FF, k, n,lower.tail = FALSE) }
 # pvalues <- sapply(Fvalues^.5,p.from.t,df=tot.samps-1)
  if(is.numeric(Fvalues)) {
    print(summary(Fvalues))
    pvalues <- sapply(Fvalues,p.from.f,k=n.phenos,n=round(tot.samps/n.phenos))
    if(p.values) { out <- pvalues } else { out <- Fvalues }
    if(length(snp.labels)==length(out))
    { names(out) <- snp.labels } else {
      stop("Error: analysis failure as list of SNPs did not match length of test statistics")
    }
  } else { out <- rep(1,length(Fvalues)) }
  return(out)
}


snp.select.least.assoc <- function(descr.fn,pc.to.keep=.05,dir,sample.info=NULL,do.plot=T,n.cores=1)
{
  ##association version of selecting subset of snps for PCA
  # takes the 'pc.to.keep'% least associated with #sample.info$phenotype'
  myPs <- do.quick.LRR.snp.assocs(descr.fn=descr.fn,sample.info=sample.info,dir=dir,p.values=T,n.cores=n.cores)
  myPs <- sort(myPs,decreasing=F)
  cat("\nsummary of p-values for association tests\n"); print(summary(myPs))
  if(do.plot) {
    pdf(cat.path(dir$pc,"distribution_of_pvalues_assoc.pdf"))
    plot(density(myPs),xlab="p-values of associations",main="distribution of association test results")
    abline(v=median(myPs)); dev.off()
  }
  ## HERE CONVERT F's to p values and get summary() and plot(density) to assess bias
  n.to.keep <- round(max(min(1,pc.to.keep),0)*length(myPs))
  kept.snps <- names(myPs)[1:n.to.keep]
  return(kept.snps)
}



extract.snp.subset <- function(snp.info,pc.to.keep=.05,assoc=F,autosomes=T,sample.info=NULL,n.cores=1,
                               writeResultsToFile=T,big.fn="combinedBigMat.RData",out.fn="pca.snp.subset.txt",dir)
{
  # extract small percentage of SNPs as a subset to commit to PCA to correct batch effects
  dir <- validate.dir.for(dir,c("big","ano"),warn=F)
  if(assoc & (is.null(sample.info) | ("phenotype" %in% colnames(sample.info)))) {
    ##association version of selecting subset of snps for PCA
    cat("\nSelecting SNP subset for PCA based on the ",round(pc.to.keep*100,1),"% least associated SNPS\n",sep="")
    keep.snps <- snp.select.least.assoc(sample.info=sample.info,descr.fn=big.fn,pc.to.keep=pc.to.keep,dir=dir,n.cores=n.cores)
  } else {
    cat("\nSelecting SNP subset for PCA based on a ",round(pc.to.keep*100,1),"% random selection of SNPS\n",sep="")
    #subset of SNPs to use (based on even as possible spacing thru genome)
    keep.snps <- random.spacing.snp.select(snp.info=snp.info,pc.to.keep=pc.to.keep,dir=dir,autosomes.only=autosomes)  
  }
  # write to file
  if(writeResultsToFile) {
    ofn <- cat.path(dir$ano,out.fn)
    writeLines(keep.snps,con=ofn)
    cat(paste("~wrote list of SNPs to use in the PCA to file:\n ",ofn,"\n"))
  } else {
    cat(" results not written to file\n")
    print(head(keep.snps)) ; cat("\n\t . . . .\n")
  }
  return(keep.snps)
}


batch.effects.density.plot <- function(plot.stats, grp.dat, grp.lab="Plate", npg=8, 
         ylz=c(25,25), xllz=c(-.2,.0), xlhz=c(.1,.5), dir, extrapref="")
{
  # detailed density plot for each batch (e.g, plate) for mean and DLRS distributions
  dir <- validate.dir.for(dir,c("qc.pl"),warn=F)
  #make the histos for all regions
  grps <- levels(as.factor(grp.dat))
  ng <- ceiling(length(grps)/npg)
  colz <- c("red","green","blue","orange","purple","lightblue","pink","brown")
  refcol <- "black"
  statz <- colnames(plot.stats)
  spotz <- c("topleft","topright")
  
  ofn <- paste(dir$qc.pl,grp.lab,"DistributionsMnDlrs",extrapref,".pdf",sep="")
  pdf(ofn,height=3*ng,width=10)
  
  par(mfcol= c(ng,length(statz)))
  for (tt in 1:2)
  { 
    stat <- paste(statz[tt])
    bg.ref <- density(plot.stats[[stat]])
    yl <- ylz[tt]; xll <- xllz[tt]; xlh <- xlhz[tt]
    spott <- spotz[tt]
    for (pp in 0:(ng-1))
    {
      plot(bg.ref,xlab=paste(stat,"Log-R Ratio"),
           ylab="Frequency amongst samples", col=refcol, lwd=1, ylim=c(0,yl),xlim=c(xll,xlh),
           main=paste("Distribution of",stat,"Log-R Ratio (Call Rate>=.95)"),bty="l")
      XY <- approx(x=bg.ref$x,y=bg.ref$y,n=length(bg.ref$x)*5)
      x <- XY$x; y <- XY$y
      segments(x,0,x,y,col="grey")
      lines(bg.ref, col="white", lwd=1)
      lines(bg.ref, col=refcol, lwd=1,lty="dashed")
      for (cc in 1:npg)
      {
        subrowz <- which(grp.dat==grps[(pp*npg)+cc])
        if(length(subrowz)>1) {
          lines(density(plot.stats[[stat]][subrowz]), col=colz[cc], lwd=1)
        } else {
          # warn if not on the last row of plots
          if (pp!=(ng-1)) {
          cat(" Warning: skipped plotting",stat,"series",cc,"from row",pp,"as the",grp.lab,"had less than 2 samples\n")
          }
        }
      }
      leg.its <- c(paste(grp.lab,narm(grps[(pp*npg)+(1:npg)])),"All (Ref.)")
      legend(spott,legend=narm(leg.its),
             text.col=c(colz[1:(length(leg.its)-1)],refcol),bty="n",ncol=2)
    }
  }
  dev.off()
  cat(paste("~wrote file:",ofn,"\n"))
}


lrr.boundary.scatter <- function(plot.stats,pass.tab,dir,fn.pre="",fn.suf="",...)
{
  # makes scatterplot for each combination of stats versus each other
  # requires 'plot.stats'; also uses 'pass.stats' which is the reference (may be the same as plot.stats)
  colorz <- "darkblue";   extxt <- ""
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  ## PLOT SD/MEAN/DLRS SCATTERS ##
  stat1 <- c("DLRS","GCWave","GCWave")
  stat2 <- c("Mean","Mean","DLRS")
  stat2p <- c("LRR-Mean","LRR-Mean","DLRS")
  stat1p <- c("DLRS","GCWave","GCWave")
  titx <- c("Sample LRR Mean Versus DLRS", "Sample GC Wave Factor Versus LRR Mean",
           "Sample GC Wave Factor Versus DLRS" )
  cat("~wrote file(s):\n")
  for (ss in 1:length(stat1))
  {
    if(all(c(stat1[ss],stat2[ss]) %in% colnames(pass.tab)) & all(c(stat1[ss],stat2[ss]) %in% colnames(plot.stats)))
    { 
      # all req'd data is there for this iteration!
      #SD vs MEAN #SD vs DLRS  # MEAN vs DLRS
      ofn <- cat.path(dir$qc.lrr,paste(stat1[ss],"vs",stat2[ss],sep=""),
                                         suf=fn.suf,pref=fn.pre,ext=".pdf")
      pdf(ofn)
      valz <- c(pass.tab[c("UB","LB"),stat2[ss]],pass.tab[c("UB","LB"),stat1[ss]])
      plot(plot.stats[[stat2[ss]]], plot.stats[[stat1[ss]]], pch=".",
           xlab=stat2p[ss], ylab=stat1p[ss], ...)
      abline(v=valz[1],lty="dotted")
      abline(v=valz[2],lty="dotted")
      abline(h=valz[3],lty="dotted")
      abline(h=valz[4],lty="dotted")
      bnd <- c("Upper","Lower")
      stt <- rep(c(stat2[ss],stat1[ss]),each=2)
      leg.txt <- paste(stt,bnd,"Bound:",round(valz,3))
      legend("topright",legend=c(titx[ss],
                                 leg.txt),bg="white",col=c((colorz)[1],rep("black",4)),
             lwd=c(10,1,1,1,1),lty="dotted" ,cex=.75)
      dev.off()
      cat(" ",basename(ofn),"\n")
    } else {
      warning(paste("skipped plot for ",stat1[ss],"vs",stat2[ss],"as some data was missing"))      
    }
  }
  cat("\n")
}


get.plate.lrr.stats <- function(plt,stat.table,samps=NULL,type="plate")
{
  ## table of plate-wise LRR stats
  if(!is.data.frame(plt)) { plt <- as.data.frame(plt) }
  if(colnames(plt)[1]!="id") { 
    cat(" changed column 1 (",colnames(plt)[1],") to 'id'",sep="")
    colnames(plt)[1] <- "id" 
  }
  which.col <- which(tolower(colnames(plt)) %in% type)[1]
  if(!is.na(which.col)) { 
    colnames(plt)[which.col] <- "myBatch" 
  }
  if(is.null(samps)) { samps <- plt$id }
  plt <- plt[plt$id %in% samps,]
  if(nrow(plt)<2) { stop("Error: no samples were found in plate index table") }
  plate.lrr.stats <- as.data.frame(matrix(ncol=ncol(stat.table)*2,nrow=length(unique(plt$myBatch))))
  colnames(plate.lrr.stats) <- paste(rep(colnames(stat.table),each=2),rep(c("(Av)","(SD)"),2))
  rownames(plate.lrr.stats) <- levels(as.factor(plt$myBatch))
  
  for(j in 1:ncol(stat.table))
  {
    plate.lrr.stats[,1+((j-1)*2)] <- tapply(stat.table[plt$id,j],as.factor(plt$myBatch),mean,na.rm=T)
    plate.lrr.stats[,(j*2)] <- tapply(stat.table[plt$id,j],as.factor(plt$myBatch),sd,na.rm=T)
  }
  return(plate.lrr.stats)
}


dual.cat <- function(file="",...,to.scrn=T,to.file=T) {
  # simultaneously apply the cat function to the screen and an open connection
  if(to.scrn){ cat(...,file="") }
  if(to.file){ if(isOpen(file)) { cat(...,file=file) } }
}


print.large.to.file <- function(tab,con="",digits=3,rlab="",rownums=F,...,to.scrn=T) 
{
  #print a matrix to file nicely
  linez <- print.large(tab,row=nrow(tab),col=ncol(tab),
                       rlab=rlab,digits=digits,rownums=rownums,...,ret=T)
  for (j in 1:(length(linez))) {
    dual.cat(paste(paste(linez[[j]],collapse=" "),"\n",sep=""),file=con,to.scrn=to.scrn)
  }
}


make.QC.summary.table <- function(sample.list,dir,pass.fn="PassQCSamples.txt",sum.fn="QCsummaryTable.txt",to.scrn=T)
{
  # create exclusion comparison table
  # show overlap between samples excluded by various QC metrics
  # use just annotated exclusion lists
  dir <- validate.dir.for(dir,c("ano","qc.lrr"),warn=F)
  sampwise <- sample.bad.count.table(dir,sample.list,type=2)
  keep.samples <- rownames(sampwise)[sampwise$TOTAL==0]
  ofn <- cat.path(dir$ano,pass.fn)
   writeLines(keep.samples,con=ofn)
  cat(paste("~wrote file of sample IDs passing all QC to:\n ",ofn,"\n"))
  
  # prepare index/names for loop of n x n comparisons
  tab.list <- colnames(sampwise[,-which(colnames(sampwise)=="TOTAL")])
  tl <- length(tab.list)
  # unique row at bottom will show samples excluded only on each criteria
  res.tab <- matrix(nrow=(tl+1),ncol=tl)
  # fill in QC table to show # samples excluded by combinations of metrics
  for (ii in 1:(tl+1))
  {
    for (jj in 1:tl)
    {
      if (ii!=(tl+1)) {
        res.tab[ii,jj] <- with ( sampwise, length(which(get(tab.list[ii]) & get(tab.list[jj]))) )
      } else {
        #last row, uniques counts..
        uniqz <- (rowSums(sampwise)==2)
        res.tab[ii,jj] <- with ( sampwise, length(which(uniqz & get(tab.list[jj]))) )
      }
    }
  }
  # tidy up table for display
  colnames(res.tab) <- tab.list
  row.names(res.tab) <- c(tab.list,"No Others")
  tot.uniq <- sum((sampwise$TOTAL))
  
  # make failure counts table
  rss <- rowSums(sampwise)
  failedN <- table(rss)
  failedN.tab <- cbind(failedN,100*round(failedN/sum(failedN),3))
  nmz <- names(failedN); nmz[nmz==0] <- "none"; nmz[nmz==max(4,max(rss,na.rm=T))] <- "All"
#  rownames(failedN.tab) <- c("None","1","2","3","All") #,"4"
  rownames(failedN.tab) <- nmz
  colnames(failedN.tab) <- c("count","%")
  
  if(nchar(paste(sum.fn))>1) {
    ofn <- cat.path(dir$qc.lrr,sum.fn,ext=".txt")
    my.fn <- file(ofn,open="w") # write to file
  } else {
    my.fn <- "" # write to display
  }
  
  ## Print tables to screen / and or file
  dual.cat("\ncounts\n",file=my.fn,to.scrn=to.scrn)
  print.large.to.file(res.tab,con=my.fn,to.scrn=to.scrn)  #print((res.tab))
  dual.cat("\npercent of failers\n",file=my.fn,to.scrn=to.scrn)
  fpc.tab <- ((round((res.tab/tot.uniq),3)))
  print.large.to.file(fpc.tab,con=my.fn,to.scrn=to.scrn) #print(fpc.tab)
  dual.cat("\npercent of callrate passing sample\n",file=my.fn,to.scrn=to.scrn)
  pc.tab <- (round((res.tab/length(sample.list)),3))
  print.large.to.file(pc.tab,con=my.fn,to.scrn=to.scrn) #print(pc.tab)
  # print counts table
  dual.cat("\nFrequency count of number of separate QC indices failed per sample\n",
           file=my.fn,to.scrn=to.scrn)
  print.large.to.file(failedN.tab,con=my.fn,to.scrn=to.scrn)
  
  if(my.fn!="") { close(my.fn); cat("\n~wrote QC summary to:\n ",ofn,"\n") }
  out.list <- list(keep.samples,res.tab,failedN.tab)
  names(out.list) <- c("pass.samples","qcsummary","failcounts")
  return(out.list)
}


LRR.boxplot <- function(stats.table,f.name=NULL,label="Sample")
{
  # make boxplot of LRR for each statistic in 'stats.table'; if file name present to file, 
  # else to screen; add appropriate labels
  if(!is.null(f.name)) {  pdf(f.name) }
    boxplot(stats.table,pch=".",main=paste(label,"Distributions of LRR Statistics"),
          xlab=c("LRR Distribution Metric"), ylab="Value (LRR units)",yaxp=c(-.5,1.25,7),bty="l")
  if(!is.null(f.name)) { dev.off() ;  cat(paste("~wrote boxplot:\n ",f.name,"\n")) }
}


calc.main.stats.on.big <- function(des.fn,dir,dlrs.pref="DLRS",deleteDLRS=T,skip.dlrs=F,
                                   snp.info=NULL,n.cores=1,apply.method=T,rmv.chr.bounds=T) 
{
  ## do the calculation of columnwise Mean, stdev and DLRS(optionally) for a big.matrix
  dir <- validate.dir.for(dir,c("big"),warn=F)
  must.use.package("genoset",T)
  # attach datafile bigmemory object
  bigMat2 <- getBigMat(des.fn,dir)
  
  if(!skip.dlrs & !apply.method) {
    # make DLRS bigmat
    # imply DLRS description and backing file names
    bck.fn.dlrs <- paste(dlrs.pref,"bckfile",sep="")
    des.fn.dlrs <- paste(dlrs.pref,"descrFile",sep="")
    cat("\nCreating DLRS big.matrix file for calculations\n")
    start.tm <- proc.time()
    ## note that it is sometimes quicker to do this seemingly convoluted creation of matrices
    ## than to calculate using 'apply' with an explicit DLRS function ['colsd' is really fast]
    lastr <- nrow(bigMat2)
    post.set <- sub.big.matrix(bigMat2, firstRow=2,lastRow=lastr, backingpath=dir$big)
    pre.set <- sub.big.matrix(bigMat2, firstRow=1,lastRow=(lastr-1), backingpath=dir$big)
    DLRSMAT <- big.matrix.operation(post.set,pre.set,"-",bck.fn.dlrs,des.fn.dlrs,dir=dir,low.ram=T)
    cat(" [DLRS matrix creation took",round(proc.time()[3]-start.tm[3]),"seconds]\n")
  }
  # calculate mean, sd, median of the big matrix
  cat("\nCalculating sample-wise statistics for LRR\n")
  cat(" calculating column means...")
  uu <- system.time(mN <- colmean(bigMat2,na.rm=T)) # 100 sec
  cat("took",round(uu[3],1),"seconds\n")
  cat(" calculating column SDs...")
  uu <- system.time(sD <- colsd(bigMat2,na.rm=T)) # 30 sec
  cat("took",round(uu[3],1),"seconds\n")
  
  if(!skip.dlrs)
  {
    if(apply.method)
    {
      cat(" calculating column DLRS (apply method) ...")
      dlrs.fast <- function(X,y) { sd(diff(X),na.rm=T) }
      if(rmv.chr.bounds) {
        # remove the differences from the first pos of chrN+1 to last of chrN
        dlrs.fast <- function(X,y) { sd(diff(X)[-y],na.rm=T) }
        must.use.package("genoset",T); 
        if(is.null(snp.info)) { snp.info <- read.snp.info(dir) }
        ind <- match(rownames(bigMat2),rownames(snp.info))
        ind.to.rmv <- chrIndices(toGenomeOrder(snp.info[sort(ind),]))[,"offset"]
        ind.to.rmv <- ind.to.rmv[ind.to.rmv!=0]
      } else {
        dlrs.fast <- function(X,y=NULL) { sd(diff(X),na.rm=T) }
        ind.to.rmv <- NULL
      }
      cmp.there <- require(compiler); if (cmp.there) { dlrs.fast <- cmpfun(dlrs.fast) }
      if(n.cores>1) {
        uu <- system.time(dlrs <- bmcapply(bigMat2,2,FUN=dlrs.fast,dir=dir$big,n.cores=n.cores,y=ind.to.rmv)) # ? sec
      } else {
        uu <- system.time(dlrs <- apply(bigMat2,2,dlrs.fast,y=ind.to.rmv)) # ? sec
      }
      dlrs <- dlrs/sqrt(2)
      cat("took",round(uu[3],1),"seconds\n")
    } else {
      cat(" calculating column dlrs (big.matrix method ...\n")
      uu <- system.time(dlrs <- colsd(DLRSMAT,na.rm=T)) # 80 sec
      cat("took",round(uu[3],1),"seconds\n")
      dlrs <- dlrs/sqrt(2)
      if(deleteDLRS)
      {
        cat(" deleting DLRS bigmemory objects to save disk space...\n")
        bfn <- cat.path(dir$big,des.fn.dlrs)
        unlink(bfn); cat("deleted:",bfn,"\n")
        bfn <- cat.path(dir$big,bck.fn.dlrs)  
        unlink(bfn); cat("deleted:",bfn,"\n")
      } else {
        cat(" manually delete DLRS bigmemory objects to save disk space if desired\n")    
      }
    }
    out.mat <- cbind(mN,sD,dlrs)
    colnames(out.mat) <- c("Mean","StDev","DLRS")
  } else {
    out.mat <- cbind(mN,sD)
    colnames(out.mat) <- c("Mean","StDev")
  }
  rownames(out.mat) <- colnames(bigMat2)
  return(out.mat)
}


## hidden
big.matrix.operation <- function(bigM1,bigM2,operation="-",bck,des,dir,split.to=40,low.ram=F)
{
  # apply simple arithmetic functions to a pair of big matrices
  # splits the operation into 'split.to' sets of rows. only need to increase this if
  # the datafile is more than 20 times the available RAM memory.
  dir <- validate.dir.for(dir,c("big"),warn=F)
  if(try(require(bigalgebra))) {
    cat(" using bigalgebra method to reduce RAM requirement for bigmatrix\n")
    options(bigalgebra.mixed_arithmetic_returns_R_matrix=FALSE)
    use.big.al <- T
  } else {
    use.big.al <- F
    cat(" bigalgebra not installed, may increase RAM requirement for DLRS matrix\n")
  }
  if (!operation %in% c("+","-","/","*"))
  { warning("operation",operation,"may not be supported") }
  
  if(all(dim(bigM1)==dim(bigM2)))
  {
    nR <- nrow(bigM1); nC <- ncol(bigM1)
    cat(" applying function to matrices. expect this to take some time\n")
    #bigR <- eval(call("-",bigM1[1:nR,1:nC],bigM2[1:nR,1:nC]))
    bigO <- big.matrix(nR,nC, backingfile=bck,
                       backingpath=dir$big, descriptorfile=des)
    stepz <- round(seq(from=1,to=nR+1,length.out=round((split.to+1))))
    if((tail(stepz,1)) != nR+1) { stepz <- c(stepz,nR+1) }
    split.to <- length(stepz)-1
    for (cc in 1:split.to)
    {
      x1 <- stepz[cc]; x2 <- ((stepz[cc+1])-1)
      if(use.big.al & operation=="-") {
        bigO[x1:x2,1:nC] <- bigM2[x1:x2,1:nC] - bigM1[x1:x2,1:nC]
      } else {
        bigO[x1:x2,1:nC] <- eval(call(operation,bigM1[x1:x2,1:nC],bigM2[x1:x2,1:nC]))
      }
      loop.tracker(cc,split.to)
      if(cc %% 4 == 0) {
        # reset memory after every 4 chunks to prevent skyrocketing RAM use #
        fl.suc <- flush(bigO) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()  
        if(low.ram) {
          RR <- describe(bigO); rm(bigO); bigO <- attach.big.matrix(RR,path=dir$big)
        }
      }
    }
    cat(" function complete, converting result to big matrix\n")
    return(bigO)
  } else {
    cat(" Warning: matrix dimensions don't match\n")
    return(NULL)
  }
}


do.venn.QC <- function(tab,dir,pref="LRR_VennDiagram") 
{
  # for up to 3 columns of sample statistic pass/failures in a table, do venn diagram
  must.use.package("limma",bioC=T)
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  if(ncol(tab)>=2)
  {
    # set to all ids to get in bottom right number NOT excluded:
    vc <- vennCounts(tab[,1:min(ncol(tab),3)])

      ofn <- paste(dir$qc.lrr,pref,paste(colnames(vc)[-ncol(vc)],collapse="+"),".pdf",sep="")
      pdf(ofn)
      vennDiagram(vc,main=c(rep(" ",5),"Degree of overlap between samples","outside bounds for LRR QC Metrics"))
      text(2.2,-2.7,"# not rejected")
      dev.off()
    cat(paste("~wrote file:",ofn,"\n"))
  } else {
    cat("Error: could not generate Venn diagram, input not a table of >=2 columns")
  }
}


trim.plate.duplicates <- function(plate.lookup,by=c("counts","first")[1])
{
  ## remove duplicate entries in a plate.lookup index dataframe that
  # was generated by 'get.plate.info'
  if(!is.list(plate.lookup) | length(plate.lookup)!=3) {
    warning("expecting a plate.lookup object, a list of length 3")
    cat(" with 1) index; 2) duplicates list (may be empty); 3) counts per plate\n")
    cat("returning unchanged lookup object, made no attempt to detect or remove duplicates\n")
    return(plate.lookup)
  } 
  if(!is.null(plate.lookup[[2]])) {
    pre.len <- nrow(plate.lookup[[1]])
    uniqz <- unique(plate.lookup[[2]][,1])
    cntr <- 0
    for (cc in 1:length(uniqz))
    {
      sub.sel <- which(plate.lookup[[2]][,1]==uniqz[cc])
      sub.tab <- plate.lookup[[2]][sub.sel,]
      if(by=="counts") {
        # keep plate assignment with highest number of samples on it
        # (idea is that this increases the probability this was the true source plate)
        if(sub.tab$count[2]>sub.tab$count[1]) { lose <- 2 } else { lose  <- 1 }
      } else {
        lose <- 2 # i.e, do it by order, keep first occuring plate assignment for each id
      }
      one.to.kill <- paste(plate.lookup[[2]][sub.sel[lose],1:2])
      to.kill <- which(plate.lookup[[1]][,1]==one.to.kill[1] & plate.lookup[[1]][,2]==one.to.kill[2])[1]
      if(length(to.kill)==1) {
        plate.lookup[[1]] <- plate.lookup[[1]][-to.kill,]
        cntr <- cntr+1
      } else {
        stop("Error: duplicated plate not found in master list")
      }
    }
    len.chng <- (pre.len-nrow(plate.lookup[[1]]))
    if(cntr==len.chng) {
      cat("To ensure unique id/plate mapping,",cntr,"duplicates were removed from plate table\n")
      cat("[For best accuracy, you should instead manually remove duplicates from the\n")
      cat("source 'plate.lookup.txt' file based on verified annotation]\n")
    } else {
      cat("warning: number of duplicate plates removed",cntr,
          "doesn't match change in table size",len.chng,"\n")
    }
    plate.lookup[[2]] <- plate.lookup[[3]] <- NULL # clear prev info - order important
    plate.lookup[[3]] <- table(plate.lookup[[1]]$plate)
    return(plate.lookup)
  } else {
    #no duplicates
    return(plate.lookup)
  }
}


get.plate.from.sample.info <- function(sample.info,other.cols=NULL) {
  # extract plate.info object from a sample info object
  if(!is.null(sample.info) & is.data.frame(sample.info)){
    pl.col <- which(tolower(colnames(sample.info)) %in% "plate")
    wl.col <- which(tolower(colnames(sample.info)) %in% "well")
    if(length(pl.col)>0){
      plate.lookup <- cbind(rownames(sample.info),sample.info[,pl.col])
      if(length(wl.col)>0) { plate.lookup <- cbind(plate.lookup,sample.info[,wl.col]) }
      colnames(plate.lookup) <- c("id","plate","well")[1:ncol(plate.lookup)]
    } else {
      warning("couldn't find plate column in sample.info")
      return(NULL)
    }
    if(is.character(other.cols)) {
      which.valid <- other.cols %in% colnames(sample.info)
      if(any(which.valid)){
        valid.cols <- unique(other.cols[which.valid])
        for (cc in 1:length(valid.cols)) {
          plate.lookup[[valid.cols[cc]]] <- sample.info[,valid.cols[cc]]
          cat(" added",valid.cols[cc],"to batch effects lookup table\n")
        }
      }
    }
    return(as.data.frame(plate.lookup))
  }  
  warning("sample.info was invalid, return NULL instead of plate.lookup")
  return(NULL)
}


get.plate.info <- function(dir, id.col=1, plate.col=2, well.col=NA,fn="plate.lookup.txt",
                           dup.action=c("","trim","print","return"),prev=F,verbose=T)
{
  ## get a table of plate information from a plate lookup file.
  # fairly flexible but lookup must contain plate ids alongside sample ids as
  # a minimum requirement. may include additional information such as 'well'
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  if(fn %in% list.files(dir$ano))
  {
    plate.lookup <- reader(cat.path(dir$ano,fn))  #read.delim  ,header=T,stringsAsFactors=F)
    if(verbose) { cat("\nRetrieving Plate/Well information from lookup table:\n") }
    if(prev) { print(head(plate.lookup,4)) }
    cN <- colnames(plate.lookup)
    any.ids <- toupper(cN) %in% c("ID","IDS","SAMPLE","SAMPLES","SAMPLEID","SAMPLE-ID","SAMPLE.ID","SAMPLE_ID")
    any.plates <- toupper(cN) %in% c("PLATE","PLATES","PLATE.ID","PLATEID","PLATE_ID","PLATE-ID")
    any.wells <- grep("WELL",toupper(cN),fixed=T)
    if(length(which(any.ids))>0) { id.col <- which(any.ids)[1] }
    if(length(which(any.plates))>0) { plate.col <- which(any.plates)[1] }
    if(length(any.wells)>0) { well.col <- any.wells[1] }
    select.cols <- c(id.col,plate.col)
    if(verbose) { cat(" taking id, plate from columns:",paste(cN[select.cols],collapse=", "),"\n") }
    if (!is.na(well.col)) { select.cols <- c(select.cols,well.col) }
    plate.lookup <- plate.lookup[,select.cols]
    colnames(plate.lookup) <- c("id","plate","well")[!is.na(c(T,T,well.col))]
    counte <- table(plate.lookup[,"plate"])
    plate.lookup[["count"]] <- (counte[match(plate.lookup[,"plate"],names(counte))])
    if(!is.data.frame(plate.lookup)) { plate.lookup <- as.data.frame(plate.lookup) }
    if(anyDuplicated(plate.lookup[,1]))
    {
      warning("duplicate records found in plate info")
      dupz <- plate.lookup[(plate.lookup[,1] %in% plate.lookup[duplicated(plate.lookup[,1]),1]),]
      if (dup.action[1]=="print") { print(dupz) }
      if (dup.action[1]=="return") { return(list(plate.lookup,dupz,counte)) }
      if (dup.action[1]=="trim") {
        plate.obj <- list(plate.lookup,dupz,counte)
        plate.obj <- trim.plate.duplicates(plate.obj,by="counts")
        return(plate.obj)
      }
    } 
    plate.obj <- list(plate.lookup,NULL,counte)
    return(plate.obj)
  } else {
    if(auto.mode()) {
      cat("\nUsing plate QC in this pipeline requires an annotation ")
      cat("file containing each sample id and the plate id. This can be ")
      cat("generated using something like the script in getPlateSupport.sh. ")
      cat("Then ensure the correct column numbers are entered for the function ")
      cat("'get.plate.info' to process this annotation file.\n")
      stop(paste("Error: could not find plate file",fn,"in",dir$ano)) 
    } else {
      warning("plate.lookup.txt file not in dir$ano, could potentially cause problems")
      return(NULL)
    }
  }
}


excl.bad.plates <- function(plate.qc.tab,plate.lookup,dir,badPlateThresh=0.33,
                            writeExclList=T,append=T,batch="plate",incl.cr=T)
{
  ## get list of plates with more than 'badPlateThresh' percent of samples failing on
  # at least one QC measure. return list of plates, and samples in these plates, by
  # default write the sample information to an exclusion file 'BadPlates.txt'
  taken.out.list <- plate.excl <- NULL
  dir <- validate.dir.for(dir,c("ano","excl"))
  if(incl.cr) {
    bad.plts <- which(plate.qc.tab[["PC+CR"]]>badPlateThresh)
    cr.txt <- " (including callrate) "
  } else {
    bad.plts <- which(plate.qc.tab[["PC"]]>badPlateThresh)
    cr.txt <- " (excluding callrate)"
  }
  n.bad <- length(bad.plts)
  if((!batch %in% colnames(plate.lookup)) & ("plate" %in% colnames(plate.lookup))) {
    warning(paste("couldn't find batch=",batch,", changed to 'plate'",sep=""))
    batch <- "plate" # set to plate as default if 'batch' is not in the lookup
  }
  if(!is.data.frame(plate.lookup)) { plate.lookup <- as.data.frame(plate.lookup) }
  cat("\nTesting for ",batch,"s where more than ",round(100*badPlateThresh,1),"% of samples fail QC",cr.txt,"\n",sep="")
  if(n.bad>0)
  {
    plate.excl <- plate.qc.tab$ID[bad.plts]
    cat(" found ",n.bad," bad ",batch,"s\n",sep="")
    for (ii in 1:n.bad)
    {  
      nxt.plt <- plate.excl[ii] 
      in.plt <- (plate.lookup[,"id"])[plate.lookup[,batch] %in% nxt.plt]
      taken.out.list <- c(taken.out.list,in.plt)
      cat(paste(" removed",length(in.plt),"samples from plate",nxt.plt,"\n"))
    }
  }
  # re-write to original file name
  if(writeExclList & dir$ano!="") {
    bad.nm <- paste("Bad",toheader(batch),"s.txt",sep="")
    ofn <- cat.path(dir$excl,bad.nm)
    # write sample exclude list to annotation directory
    if(!append) {
      writeLines(taken.out.list,con=ofn) 
      cat("\n~wrote bad",batch,"excluded subject list to:\n ",ofn,"\n")
    } else {
      if(file.exists(ofn)) {
        exist.lns <- readLines(ofn) ####
      } else {
        exist.lns <- NULL
      }
      to.wrt <- paste(unique(c(exist.lns,taken.out.list)))
      writeLines(to.wrt,con=ofn) 
      cat(paste("~appended bad",batch,"excluded subject list:\n ",ofn,"\n"))     
    }
  } else {
    cat(batch,"excl. file not written\n") 
    print(head(sub.left)) ;cat("\n\t . . . .\n")
  }
  outlist <- list(taken.out.list,plate.excl)
  names(outlist) <- c("samples.from.bad.plates","bad.plates")
  return(outlist)
} 


chr.ab.report <- function(chr.stat,chrWarns,dir,writeExclList=F,makeGraphs=F,pref="",append=F)
{
  # make chromosomal abberation graphs, write sample exclusion file 'ChrAb.txt' to dir$ano
  nC <- ncol(chr.stat$chr.mean)
  dir <- validate.dir.for(dir,c("ano","qc.cab","excl"))
  excluded.samps <- sort(table(unlist(chrWarns)))
  # put stats in summary table for display
  chr.stats.table <- data.frame(LRR.Mean.Av=round(colMeans(chr.stat$chr.mean,na.rm=T),3),
                                LRR.Mean.Sd=round(apply(chr.stat$chr.mean,2,sd,na.rm=T),3),
                                LRR.StDev.Av=round(colMeans(chr.stat$chr.sd,na.rm=T),3),
                                LRR.StDev.Sd=round(apply(chr.stat$chr.sd,2,sd,na.rm=T),3))
  rownames(chr.stats.table) <- paste("Chr",1:nC,sep="")
  cat("\nTable of LRR-Mean and LRR-SD for each chromosome [headers are chr # in file]\n")
  print(chr.stats.table,digits=3)
  
  # print number of outliers found for each chromosome
  cat("\nNumber of outliers per chromosome:\n")
  print((sapply(chrWarns,length)))
  cat("Number of unique samples with warnings:",length(unique(unlist(chrWarns))),"\n")
  cat("Total number of warnings:",length((unlist(chrWarns))),"\n")
  
  cat("\nHistogram of number of chromosome outliers per sample:\n")
  textogram(excluded.samps)
  
  # create box plot of distributions per chromosome
  if(makeGraphs) {
    ofn <- paste(dir$qc.cab,"ChrMeansBoxPlot",pref,".pdf",sep="")
    pdf(ofn)
    boxplot(chr.stat$chr.mean,pch=".",xlab="Chromosome number",ylab="LRR-Mean",
            main="Boxplot of LRR-Mean distributions for each chromosome",bty="l")
    dev.off()
    cat(paste("~wrote plot:",ofn,"\n"))
  }
  if(writeExclList & dir$ano!="") {
    ofn <- cat.path(dir$excl,"ChrAb.txt")
    # write sample exclude list to annotation directory
    if(!append) {
      writeLines(names(excluded.samps),con=ofn) 
      cat("\n~wrote chr.ab excluded subject list to:",ofn,"\n")
    } else {
      if(file.exists(ofn)) {
        exist.lns <- readLines(ofn) ####
      } else {
        exist.lns <- NULL
      }
      to.wrt <- unique(c(exist.lns,names(excluded.samps)))
      writeLines(to.wrt,con=ofn) 
      cat(paste("~appended chr.ab excluded subject list:\n ",ofn,"\n"))     
    }
  }
}


get.chr.stats <- function(bigMat,snp.info,dir="",allow.subset=F)
{
 # get the LRR mean and SD for each chromosome (autosome) separately
 # used for detection of chromosomal abberations
 must.use.package(c("bigmemory","biganalytics"))
 must.use.package("genoset",T)
 dir <- validate.dir.for(dir,c("big"),warn=F)
 if((is(snp.info)[1]!="RangedData") | (is(bigMat)[1]!="big.matrix") )
 { stop("Error: parameter 'snp.info' should be RangedData, and 'bigMat' should be a big.matrix") }
 # Start making histograms of LRR-Mean for each chromosome
 snp.info <- snp.info[(rownames(snp.info) %in% rownames(bigMat)),]
 if (allow.subset) { 
   if(nrow(snp.info)!=nrow(bigMat)) {
     cat(" modifying big.matrix to match snp.info object\n")
     rN <- rownames(bigMat)[rownames(bigMat) %in% rownames(snp.info)]
     big.descr <- big.exclude.sort(bigMat,f.samp=NULL,f.snp=rN,tranMode=1,pref="TEMP_",dir=dir,verbose=F)
     rm(bigMat); bigMat <- getBigMat(big.descr,dir$big)
   }
 } else {
   if(nrow(snp.info)!=nrow(bigMat)) { stop("Error: bigMat contained snps not in snp.info") }
 }
 chr.set <- chr.nums.of.info(snp.info)
 nC <- length(chr.set); range.chr <- 1:nC
 dns <- list(colnames(bigMat),paste("chr",chr.set,sep="")) # get rownames,colnames
 # create matrices for LRR-stats
 chr.mean <- chr.sd <- matrix(nrow=ncol(bigMat),ncol=nC,dimnames=dns) # used to be list()
 indx.first <- chrIndices(snp.info)[,"first"]
 indx.last <- chrIndices(snp.info)[,"last"]
 cat(" processing chr: ")
 for (dd in range.chr) 
 { 
  cat(chr.set[dd],"..",sep="")
  LRR.dat <- sub.big.matrix(bigMat, firstRow=indx.first[dd], lastRow=indx.last[dd], backingpath=dir$big)
  chr.mean[,dd] <- colmean(LRR.dat,na.rm=T) # 50 sec
  chr.sd[,dd] <- colsd(LRR.dat,na.rm=T)
 }
 cat("done\n")
 out.list <- list(chr.mean,chr.sd)
 names(out.list) <- c("chr.mean","chr.sd")
 return(out.list)
}
 
 
make.chr.fail.table <- function(chr.mean,pctile.bound=.01)
{
 # using chromosome (autosome) means from 'get.chr.stats', create list of
 # failures per sample for their chr.mean to within confidence limits
 range.chr <- 1:ncol(chr.mean)
 #chr.set <- chr.nums.of.info(snp.info)
 #n.chr <- length(chr.set); range.chr <- 1:n.chr
 headz <- c("ID","LRR.Mean","LoHi","ChrNum")
 listOfFails <- matrix(ncol=4)
 colnames(listOfFails)<- headz
 for (dd in range.chr) 
 { 
  badcuts <- pctile(chr.mean[,dd],pctile.bound)
  lo.peeps <- rownames(chr.mean)[chr.mean[,dd] <= badcuts[1]]
  hi.peeps <- rownames(chr.mean)[chr.mean[,dd] > badcuts[2]]
  if(length(lo.peeps)<1 | length(hi.peeps)<1) {
    cat("no top percentiles identified for chr col",dd,"as too few datapoints\n")
  } else {
    new.chunk <- rbind(cbind(lo.peeps,chr.mean[lo.peeps,dd],
          rep("low",length(lo.peeps)),rep(dd,length(lo.peeps))),
          cbind(hi.peeps,chr.mean[hi.peeps,dd],
          rep("high",length(hi.peeps)),rep(dd,length(hi.peeps))) )
    colnames(new.chunk) <- headz # see def'n above
    listOfFails <- rbind(listOfFails, new.chunk)
    if(is.na(listOfFails[1,1])) { listOfFails <- listOfFails[-1,]}
    # record list of all samples exceeding 'pctile.bound' upper or lower, per chr.
  }
 }
 return(listOfFails)
}

 
do.chr.ab.contrasts <- function(chr.mean,lob=2,hib=2.5,pctile.bound=.01,nC=22)
{
 ## Set up contrast matrices to test for chromosomal aberrations (autosomes) ##
 ## two matrices, one a simple average comparison, target vs rest,
 ## the other a local comparison, target versus adjacent weighted
 ## 2^n
 # generate local chromosomal contrast matrix
 # comparison of target chromosome to weighted mean of adjacent chromosomes
 if(ncol(chr.mean)!=nC) { stop("Object should contain 1 column for each autosome [n=1..",nC,"]") }
 chromocontr <- matrix(integer(),nrow=nC,ncol=nC)
 # weightings for the 3 chromosomes to the left/right of target chromosome
 i <- c(1,2,4)
 versuz <- c(sum(i),sum(i)+i[3],(2*sum(i))-i[1],rep(12,16),(2*sum(i))-i[1],sum(i)+i[3],sum(i))
 stepo <- list(NULL,-i[3],c(-i[2],-i[3]),c(-i[1],-i[2],-i[3]),c(-i[3],-i[2],-i[1]),c(-i[3],-i[2]),-i[3],NULL)
 for (jj in 1:nC) {
   chromocontr[jj,] <- c( rep(0,max(0,jj-4)),stepo[[max(1,min(jj,4))]],
                    versuz[jj], stepo[[max(5,min(jj-14,8))]],
                     rep(0,max(0,19-jj)) ) /versuz[jj]
 }
 # generate global chromosomal contrast matrix
 # simple comparison of target chromosome to mean of all others
 chromocontr2 <- matrix(-1/(nC-1),nrow=nC,ncol=nC)
 diag(chromocontr2) <- 1
 
 # calculate contrast matrices using matrix multiplication #
 # LOCAL
 chr.conts <- t(chromocontr %*% t(chr.mean))
 colnames(chr.conts) <- paste("chr",1:nC,sep="")
 # GLOBAL
 chr.conts2 <- t(chromocontr2 %*% t(chr.mean))
 colnames(chr.conts2) <- paste("chr",1:nC,sep="")
 
 listOfFails <- make.chr.fail.table(chr.mean)	

 # standardize contrasts to z scores to help determine outliers
 chr.contsZ <- apply(chr.conts,2,standardize)
 chr.contsZ2 <- apply(chr.conts2,2,standardize)

 cat("\nTable of Exclusion processing per chromosome...\n")
 chrWarns <- list()
 for (cc in 1:nC)
 {
  # get T/F vector for all samples if contrast values are HIGH (above thresh)
  hiz <- (chr.contsZ[,cc] > hib & chr.contsZ2[,cc] > lob) | (chr.contsZ2[,cc] > hib & chr.contsZ[,cc] > lob)  
  # get T/F vector for all samples if contrast values are LOW (below thresh)
  loz <- (chr.contsZ[,cc] < -hib & chr.contsZ2[,cc] < -lob) | (chr.contsZ2[,cc] < -hib & chr.contsZ[,cc] < -lob)
 # lookup whether samples also failed on initial raw outlier (percentile) status
  failz <- rownames(chr.contsZ) %in% listOfFails[as.numeric(listOfFails[,"ChrNum"])==cc,"ID"]
  # select list of bad subjects for current chromosome as those failing
  # either the high (hiz) or low (loz) threshold, plus are an outlier in absolute terms (failz)
  cat(paste("chr:",cc,"too high:",length(which(hiz)),"too low:",length(which(loz)),
              "excluded top/bottom",pctile.bound*100,"%:",length(which(failz))),"\n")
  chrWarns[[cc]] <- rownames(chr.contsZ)[which( (hiz | loz) & failz)]
 } 
 # set names of warning list components to allow determination of which
 # chromosome the warning came from once it's been converted to a vector
 names(chrWarns) <- paste("chr",1:nC,"_",sep="")
 return(chrWarns)
}


do.med.chunk <- function(ranges,use.big.list,bigMat,dir,ncb,cont) {
  ## internal function for do.median.for.ranges() to calculate contents of 1 loop
  if(cont) {
    # faster method (when each range is continuous)
    r1 <- ranges[1]; r2 <- ranges[2]
    if(use.big.list) {
      # this provides large speed increase when chunks have >200 snps
      LRR.chunk <- sub.big.matrix(bigMat, firstRow = r1, lastRow = r2, backingpath=dir$big)
    } else {
      LRR.chunk <- bigMat[c(r1:r2),c(1:ncb)]
    }
  } else {
    # method if intervals can be discontinuous [e.g, if dataset not sorted by chr,pos] #
    LRR.chunk <- bigMat[ranges,c(1:ncb)] # this seems wrong?
  }
  return(apply(LRR.chunk,2,median,na.rm=T))
}


do.median.for.ranges <- function(ranges.list,bigMat,dir,cont=T,
  use.big.list=NULL, med.sav.nm="", n.cores=1)
{
  # the diskin method to calculate GC requires sample-wise medians for snps
  # in each 1MB window of the genome
  dir <- validate.dir.for(dir,c("big"),warn=F)
  if((is(ranges.list)[1]!="list") | (is(bigMat)[1]!="big.matrix") )
  { 
    if((is(bigMat)[1]=="matrix")) {
      cat(" for GC calc, bigMat was a matrix, attempting to convert to big.matrix...")
      bigMat <- as.big.matrix(bigMat)
      cat("ok\n")
    } else {
      stop("Error: parameter 'ranges.list' should be a list, and 'bigMat' should be a big.matrix") 
    }
  }
  num.ranges <- length(ranges.list)
  ncb <- ncol(bigMat)
  cat(" doing median calculation\n")
  med.store <- matrix(nrow=num.ranges,ncol=ncb)
  jj <- proc.time()
  if(n.cores>1) { multi <- T } else { multi <- F }
  if(multi) { job.count <- 0 ; cc.collect <- numeric(n.cores)}  # start a counter for parallels (if multi)
  # sort numranges if using multicore (ensures multi processors are working on as similar size chunks as possible)
  lens <- sapply(ranges.list,diff); sort.lens <- order(lens)
  ##unsort.lens <- order(((1:length(lens))[order(lens)])) # not needed
  ####
  runs <- vector("list",n.cores)
  if(multi) {
    for (cc in 1:num.ranges) {
      job.count <- job.count + 1
      cc_s <- sort.lens[cc]
      cc.collect[job.count] <- cc_s
      runs[[job.count]] <- multicore::parallel(do.med.chunk(ranges.list[[cc_s]],(cc_s %in% use.big.list),bigMat,dir,ncb,cont))
      if(job.count>=n.cores | cc>=num.ranges) {
        ## collect previous #'n.cores' runs 
        list.set <- collect(runs[1:job.count])
        for (dd in (which(cc.collect>0))) {
          med.store[cc.collect[dd],] <- list.set[[dd]]
        }
        job.count <- 0; runs <- vector("list",n.cores); cc.collect <- numeric(n.cores)
      }  
      loop.tracker(cc, num.ranges)
    }
  } else {
    for (cc in 1:num.ranges) {
      med.store[cc,] <- do.med.chunk(ranges.list[[cc]],(cc %in% use.big.list),bigMat,dir,ncb,cont)
      # add a tracker, saving backups as we go
      loop.tracker(cc, num.ranges, jj,  sav.obj=med.store, sav.fn=med.sav.nm)
    }
  }
  kk <- proc.time()
  cat(" took",round((kk-jj)[3]/3600,2),"hours\n")
  ## save median matrix to pick up later
  if(!(med.sav.nm=="")) {
    save(med.store,file=med.sav.nm) }
  return(med.store)
}


get.ranges.for.median <- function(snp.info,gc.dat,bigMat,min.size=10)
{
 # get list of position ranges for snps in each 1MB window of the genome,
 # ignoring any window with less than 'min.size' snps
 must.use.package(c("genoset","BiocGenerics"),T); must.use.package("bigmemory")
 if((is(gc.dat)[1]!=is(snp.info)[1]) | (is(gc.dat)[1]!="RangedData") )
 { stop("Error: parameters snp.info and gc.dat should be of type 'RangedData'") }

 olp <- findOverlaps(snp.info,gc.dat)

 # get list of snp.info row positions within each gc.dat range
 mbGroups <- tapply(queryHits(olp),subjectHits(olp),c)
 # get number of snps in each of these
 mbGroupsC <- sapply(mbGroups,length)
 too.small <- which(mbGroupsC<min.size)
 if(any(too.small)) { mbGroups <- mbGroups[-too.small]; mbGroupsC <- mbGroupsC[-too.small] }
 if(length(mbGroups)<2) { stop("Error: insufficient snp density for GC calculations. Check data or turn off GC") }
 
 gc.dat.filt <- subsetByOverlaps(gc.dat,snp.info)
 gc.dat.filt <- gc.dat.filt[-too.small,]
 
 # Prepare list of ranges for median calcs (depends on whether ranges will be continuous)
 if( any(rownames(snp.info)!=rownames(bigMat)) )
 { 
  warning("snp.info object rownames don't match big.matrix rownames.\n",
  "suggest using big.exclude.sort() function to reduce the datafile to only the \n",
  "needed snp-set. Otherwise this routine will be greatly slowed by using name indexing.")
  # identify any ranges that will cause a large enough block of memory allocation
  # that should use big memory subset instead of standard R object for subset chunk
  cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
  max.gb <- 1.5  #maximum number of gigabytes of chunk before might cause probs
  bad.row.max <- round(max.gb*((cells.per.gb/as.double(ncol(bigMat)))))
  name.ind <- T
  # convert position info to snp labels
  mbGroupsN <- lapply(mbGroups,function(X,rn) { rn[X] },rn=rownames(snp.info))
  ranges.list <- lapply(mbGroups,function(X,rn) { narm(match(X,rn)) },rn=rownames(bigMat))
  use.big.chunk <- which(sapply(ranges.list,length)>bad.row.max)
  if(any(use.big.chunk)) 
  { 
   stop("Error: some ranges are too large for name indexing, regenerate datafile as suggested above")
  }
 } else {
  #snp.info rownames match bigmatrix snps exactly
  cat(" 100% match between snp.info annotation and big.matrix file\n")
  name.ind <- F
  ranges.list <- lapply(mbGroups,range)
  row.max <- 200  
  use.big.chunk <- which(sapply(ranges.list,diff)>row.max)
 }
 outlist <- list(ranges.list,use.big.chunk,name.ind,gc.dat.filt)
 names(outlist) <- c("ranges.list","big.list","name.ind","gc.dat") 
 return(outlist)
}


remove.qc.failers <- function(tab,sample.info,dir=NULL,fail.col="QCfail") 
{
  # use column 'QCfail' in sample.info to remove failing ids from any frame
  # sample.info could be any dataframe where rownames are sample ids.
  dat <- force.frame(tab) # allows flexible input, forces dataframe result
  if(!fail.col %in% colnames(sample.info)) 
  { 
    if(!is.null(dir)) {
      warning(paste("sample.info object did not have column '",fail.col,"'. Adding it now..",sep=""))
      sample.info <- validate.samp.info(sample.info,QC.update=T,dir=dir,verbose=F)
    } else {
      stop(paste("Error: sample.info object did not have column '",fail.col,"'",sep=""))
    }
  } 
  pass.ids <- rownames(sample.info)[sample.info[[fail.col]]==0]
  index <- find.id.col(dat,ids=pass.ids)$index
  return(dat[index,])
}


samp.update.qc.fail <- function(sample.info,dir,verbose=F,proc=1) {
  dir <- validate.dir.for(dir,c("lrr.dat","ids","ano","excl"))
  # update QC failure flags based on sample exclusion directory
  if(!( "QCfail" %in% colnames(sample.info)))
  {
    if(verbose) { cat(" adding qc failure flag column to sample.info\n") }
    sample.info[["QCfail"]] <- rep(0,times=nrow(sample.info))
    QC.update <- T
  }
  failed.samps <- get.all.samp.fails(dir)
  if(verbose) { cat(" updating QC failure flags to reflect lists in '/SAMPLE_EXCLUDE/'\n") }
  sel <- rownames(sample.info) %in% failed.samps
  sample.info$QCfail[!sel] <- 0
  if(proc!=1){
    sample.info$QCfail[sel & (sample.info$QCfail==0)] <- proc # allow custom numbers showing where fail occurred
  } else { sample.info$QCfail[sel] <- 1 }
  return(sample.info)
}


snp.update.qc.fail <- function(snp.info,dir,verbose=F,proc=1) {
  dir <- validate.dir.for(dir,c("lrr.dat","cr","ano","excl2"))
  # update QC failure flags based on snp exclusion directory
  if(is(snp.info)[1]!="RangedData") { warning("not RangedData"); return(NULL) }
  if(!("QCfail" %in% colnames(snp.info)))
  {
    if(verbose) { cat(" adding qc failure flag column to snp.info\n") }
    snp.info[["QCfail"]] <- rep(0,times=nrow(snp.info))
    QC.update <- T
  }
  failed.snps <- get.all.snp.fails(dir)
  if(verbose) { cat(" updating QC failure flags to reflect lists in '/SNP_EXCLUDE/'\n") }
  sel <- rownames(snp.info) %in% failed.snps
  snp.info$QCfail[!sel] <- 0
  if(proc!=1){
    snp.info$QCfail[sel & (snp.info$QCfail==0)] <- proc # allow custom numbers showing where fail occurred
  } else { snp.info$QCfail[sel] <- 1 }
  return(snp.info)
}


sync.samp.info.with.file.spec <- function(sample.info,dir) {
  sub.id.list <- get.subIDs(dir,"all")
  ID.lists <- sub.id.list$original
  file.info <- get.file.specs(dir)
  file.info <- column.salvage(file.info,"GRP",c("grp","Grp","Group","group","GROUP"))
  
  if(length(unique(narm(file.info$GRP)))!=length(unique(narm(sample.info$grp))))
  {
    warning(paste("number of 'grp's in sample info clashes with file.spec.txt in:\n ",dir$ids),
    " \n updating sample.info$grp to match 'GRP' coding in file.spec.txt:")
    for (jj in 1:length(ID.lists))
    { 
      fl.lst <- paste(file.info[,1])
      nfn <- (sub.id.list$files)[jj]
      next.grp.num <- match(gsub(".ids","",nfn),fl.lst)
      if(is.na(next.grp.num)) {
        next.grp.num <- match(nfn,fl.lst)
        if(is.na(next.grp.num)) { stop(paste("Error: couldn't find",nfn,"in file.spec.txt")) }
      }
      sample.info$grp[rownames(sample.info) %in% ID.lists[[jj]]] <- file.info$GRP[next.grp.num]
    }
  }
  return(sample.info)
}


validate.samp.info <- function(sample.info,dir,QC.update=T,file.spec=F,verbose=F,proc=1)
{
  ## make sure current sample.info object conforms to the expected structure ##
  # update QC failure flags based on sample exclusion directory
  # update 'grp' membership based on file.spec.txt column 'GRP'
  dir <- validate.dir.for(dir,c("lrr.dat","ids","ano"))
  if(!is.data.frame(sample.info)) {
    cat(" attempting to coerce sample info to data.frame\n")
    sample.info <- as.data.frame(sample.info)
  }
  if(!( "grp" %in% colnames(sample.info)))
  {
    stop("Error: sample.info must contain a column named 'grp'")
  } else {
    if(is.null(rownames(sample.info)) | all(rownames(sample.info)==paste(1:nrow(sample.info))))
    { stop("Error: sample.info must have subject Ids as rownames (unique & !=1:nrow)") }
  }
  if(file.spec & auto.mode()) {
    sample.info <- sync.samp.info.with.file.spec(sample.info,dir)  #DEP ?
  }
  if(QC.update) {
    sample.info <- samp.update.qc.fail(sample.info,dir=dir,verbose=verbose,proc=proc)
  }
  return(sample.info)
}


get.subIDs <- function(dir,ret=c("all","lists","combined","files","groups"),verbose=F,combine.grp=F)
{
  # if separate id lists for each cohort are in 'dir' then create list of these
  dir <- validate.dir.for(dir,c("ids","col","ano"))
  if(!auto.mode()) { warning("a subject id parameter has been set null, which is only valid in auto.mode") }
  if(verbose) { cat("\nSearching for subject IDs (automatic mode)\n") }
  file.info <- get.file.specs(dir)
  if(is.null(file.info)) { stop("Error: tried to automatically find subject IDs but not in auto mode") }
  id.fnz <- paste(file.info[,1],".ids",sep="")
  grpz <- unique(file.info$GRP)
  num.filz <- length(id.fnz); num.grpz <- length(grpz)
  subIDs.list <- vector("list", num.filz)
  names(subIDs.list) <- id.fnz
  path.fn <- cat.path(dir$ids,id.fnz,must.exist=T)  # sample id file name
 ## print(path.fn)
  if(verbose) { cat(paste(" reading IDs from: ",id.fnz,"\n",sep="")) }
  ID.list <- lapply(path.fn,readLines)  #readLines(sample.fn)
  cmb.ID.list <- paste(do.call("c",ID.list))
  if(tolower(ret[1]) %in% c("all","groups","group","group.files","group.file","groupfiles","groupfile")) {
    next.grp.files <- next.grp.ids <- vector("list",length(grpz))
    names(next.grp.files) <- names(next.grp.ids) <- grpz
    dd <- 1
    for (cc in grpz) {
      next.fn <- cat.path(dir$ids,id.fnz[which(file.info$GRP %in% cc)])
      next.idz <- lapply(next.fn,readLines)
      next.grp.files[[dd]] <- next.fn
      if(combine.grp) {
        next.grp.ids[[dd]] <- paste(do.call("c",next.idz)) 
      } else {
        next.grp.ids[[dd]] <- next.idz
      }
      dd <- dd + 1
    }
  }
  out.obj <- NULL
  if (tolower(ret[1]) %in% c("group","groups")) { out.obj <- next.grp.ids }
  if (tolower(ret[1]) %in% c("group.file","group.files","groupfile","groupfiles")) 
    { out.obj <- next.grp.files }
  if (tolower(ret[1]) %in% c("list","lists")) { out.obj <- ID.list }
  if (tolower(ret[1]) %in% c("combined","combine")) { out.obj <- cmb.ID.list }
  if (tolower(ret[1]) %in% c("file","files")) { out.obj <- id.fnz }
  if (tolower(ret[1]) %in% c("","all","every")) {
    out.obj <- list(cmb.ID.list,next.grp.ids,ID.list,id.fnz,next.grp.files)
    names(out.obj) <- c("combined","group","original","files","group.files")
  }
  return(out.obj)
}


get.file.specs <- function(dir,fn="file.spec.txt",quiet=T)
{ 
  # check for and read specifications file 'file.spec.txt' from directory dir$lrr.dat
  # input as 'dir' list or otherwise will try to convert to this format from character()
  dir <- validate.dir.for(dir,c("lrr.dat","ano","raw"),warn=F)
  if(!auto.mode() & !quiet) { warning("automatic search for file specs invoked, which is only valid in auto.mode") }
  if (fn %in% c(list.files(dir$lrr.dat),list.files(dir$ano)))
  {
    if (fn %in% c(list.files(dir$lrr.dat))) {
      file.info <- read.table(cat.path(dir$lrr.dat,fn),header=T) 
    } else {
      file.info <- read.table(cat.path(dir$ano,fn),header=T)
    }
    alph.ord <- order(file.info[,1])
    if (any(alph.ord!=1:nrow(file.info)))
    { cat("rearranging ",fn," into alphabetical order by filename")
      file.info <- file.info[alph.ord,]  } #alphabetical order  
  } else {
    linez <- character()
    linez[1] <- paste("If using a genome studio source file in automatic mode, the location:\n",dir$lrr.dat)
    linez[2] <- paste("should contain a file '",fn,"' which is tab delimited file with columns:\n FILE GRP TYPE SAMP SNP A1 A2 LRR BAF")
    linez[3] <- "FILE=file name (no path), GRP = 1:ngrps, TYPE=txt or gzip, then remainder are the column number "
    linez[4] <- "in each raw file that contain data for SAMP(les) SNP(s) (alleles,A1 A2) and LRR, BAF, respectively."
    linez[5] <- "\nThis file could have been used by 'bash_getAColumn.sh to extract the data initially."
    if ("raw" %in% names(dir)) {
      linez[6] <- paste("File locations should be relative to the path:",dir$raw,"\n") }
    if(!quiet) {
      cat(c(linez))
      warning("file.spec.txt missing [needed for automatic mode]")
    }
    file.info <- NULL
  }
  return(file.info)
}


get.file.lens <- function(dir,fn="file.lengths.txt",recalc=F,write=T)
{
  ## read in file lengths from text file 'file.lengths.txt' (generated during bash import)
  dir <- validate.dir.for(dir,c("lrr.dat","ano","col"),warn=F)
  if (fn %in% c(list.files(dir$lrr.dat),list.files(dir$ano)))
  {
    if(!recalc) {
      # have file lengths file, and no request to recalculate
      if (fn %in% c(list.files(dir$lrr.dat))) {
        lens.fn <- cat.path(dir$lrr.dat,fn)
      } else {
        lens.fn <- cat.path(dir$ano,fn)
      }
      cat(" reading file lengths of raw datafiles from",fn,"\n")
      test.valid <- (length(readLines(lens.fn))>0)
      if(test.valid) {
        len.tab <- read.table(lens.fn)
        len.tab <- as.data.frame(len.tab)
        colnames(len.tab) <- c("length","file.name")
      } else {
        recalc <- T
      }
    }
  } else {
    recalc <- T
  }
  if(recalc) {
    cat(" using bash command to get file lengths of:\n ",dir$col,"\n")
    cmd <- paste("wc -l ",dir$col,"*",sep="")
    linez <- system(cmd,intern=T) 
    dir.linez <- unique(c(grep("Is a directory",linez),grep("total",linez)))
    if(length(dir.linez)>0) { linez <- linez[-dir.linez] }
    if(write) { writeLines(paste(linez),cat.path(dir$lrr.dat,fn)) }
    linezp <- strsplit(rmv.spc(linez)," +")
    len.tab <- data.frame(length=sapply(linezp,"[",1),file.name=basename(sapply(linezp,"[",2)))
  }  
  alph.ord <- order(len.tab[,2])
  if (any(alph.ord!=1:nrow(len.tab)))
  { 
    cat(" rearranging",fn,"into alphabetical order by filename\n")
    len.tab <- len.tab[alph.ord,]  
  }
  return(len.tab)
}


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


get.lrr.files <- function(dir,grps=NA,...) {
  # wrapper for get.lrr.dat.file.names
  if(!auto.mode()) 
    { warning("an option for LRR file name has been set null, which is only valid in auto.mode") }
  all.fn <- get.lrr.dat.file.names(dir,...)
  file.info <- get.file.specs(dir)
  if(all(!is.na(grps)) & length(all.fn)>1) {
    select <- which(paste(file.info$GRP) %in% paste(grps))
    if(length(select)==0) {
      stop(paste("No LRR files found for group(s):",paste(grps,collapse=",")))
    } 
    return(all.fn[select])
  } else {
    return(all.fn)
  }
}


get.lrr.dat.file.names <- function(dir,suffix="LRR.dat",baf=F)
{
  ## tries really hard to find the right datafiles to read in for LRR raw data
  dir <- validate.dir.for(dir,c("col","baf.col")[1+baf])
  if(baf) { dir.col <- dir$baf.col } else { dir.col <- dir$col }
  all.in.dir <- list.files(dir.col)  
  file.info <- get.file.specs(dir) # could change to dir.col but this way more guaranteed to find file.spec.txt
  exp.fns <- paste(file.info[,1],suffix,sep=".")
  cat(" searching for raw datafiles in processed long format (e.g, converted from GenomeStudio)\n")
  if(nrow(file.info)==1) {
    # simplest case - only 1 file
    if(file.exists(cat.path(dir.col,suffix))) {
      cat(" found combined data file:",suffix,"\n")
      names.out <- suffix
    } else {
      if(file.exists(cat.path(dir.col,exp.fns))) {
        names.out <- exp.fns
        cat(" found sole data file:",exp.fns,"\n")
      } else {
        any.dats <- grep(".DAT",toupper(all.in.dir),fixed=T)
        if(length(any.dats)>0) {
          warning(paste("was expecting file named ",exp.fns,"or",suffix))
          cat(", was not in",dir.col,", but there were",length(any.dats),"*.dat files\n")
          cat("Will try to proceed with the first file:",all.in.dir[any.dats][1],"\n")
          cat("If this is not the LRR raw datafile then please CTRL-C and review file.spec.txt\n")
          names.out <- cat.path(dir.col,all.in.dir[any.dats][1])
        } else {
          if(length(all.in.dir)==1) {
            warning(paste("was expecting file named ",exp.fns," or ",suffix,
            ", these were not in\n ",dir.col,", but there was only 1 file present\n",
            "Will try to proceed with file: ",all.in.dir[1],"\n",
            "If this is not the LRR raw datafile then please CTRL-C and review file.spec.txt"))
            names.out <- cat.path(dir.col,all.in.dir[1])
          } else {
            warning(paste("found multiple files in",dir.col,",\n and could not tell which if any were the raw LRR datafile"),
            "\nPlease review file.spec.txt"); return(NULL)
          }
        }
      }
    }
    return(names.out)
  } 
  
  if(length(unique(file.info$GRP))<=1) {
    # only 1 grp, so hopefully a combined datafile already, if not need to manage that
    if(file.exists(cat.path(dir.col,suffix))) {
      cat(" found combined data file:",suffix,"\n")
      names.out <- suffix
      return(names.out)
    } else {
      any.dats <- grep(".DAT",toupper(all.in.dir),fixed=T)
      if(length(any.dats)==1) {
        warning(paste("was expecting file named",suffix,". This was not in\n"),
        paste(dir.col,", but there was 1 dat file\n",
            all.in.dir[any.dats][1],"which will be used. \nIf this is not"),
        paste(" the combined LRR raw datafile then please CTRL-C and review file.spec.txt"))
        names.out <- cat.path(dir.col,all.in.dir[any.dats][1])
        return(names.out)
      } else {
        # move onto next part, hoping now there are multiple files 
        # which are meant to be combined into 1 group
      }
    }
  }
  
  if(nrow(file.info)>1) {
    # if we've made it to here there must be more than one cohort in the files or
    #  there are multiple files that need to be combined into one group   
    ##print(file.exists(cat.path(dir$col,file.info[,1])))
    ##print(cat.path(dir$col,file.info[,1]))    
    exp.fnp <- cat.path(dir.col,exp.fns)
    if (any(file.exists(exp.fnp)))
    {
      if (all(file.exists(exp.fnp)))
      {
        cat(" all",length(exp.fnp),"expected data files found in",dir.col,"\n")
        names.out <- exp.fnp
      } else {
        which.left <- length(which(file.exists(exp.fnp)))
        cat("Error: only",length(which.left),"/",nrow(file.info),
            "of the expected datafiles were found in",dir.col,"\n")
        stop("please review file.spec.txt")
      }
      return(names.out)
    } else {
      warning(paste("none of the expected data files found in",dir.col),
      "was looking for:\n",paste(exp.fnp,collapse="\n"),"\n",
      "Please review file.spec.txt"); return(NULL)
    }
  } else {
    warning(paste("No",c("LRR","BAF")[1+baf],"datafiles found with suffix",suffix,"in \n",dir.col)
    ,"\nPlease review file.spec.txt"); return(NULL)
  }
}


time.est.snps <- function(tab,lines.per.min=12*10^6)
{
  # calculates the expected number of time to import a number of lines into snpMatrix format
  times <- round(tab[,1]/lines.per.min,1)
  times <- c(times,sum(times))
  files <- c(tab[,2],"All files")
  result <- cbind(times,files)   
  colnames(result) <- c("~ mins to process","Raw file name")
  return(result)
}


check.snp.read.fn <- function()
{
  # if read.snps.long is ever removed (function is currently deprecated) then
  # this wrapper helps suck up redundant variables and redirect relevant ones to 'read.long'
  must.use.package("snpStats",bioC=T)
  if(!exists("read.snps.long",mode="function",where="package:snpStats"))
  {
    cat("plumbCNV was designed to work with the fast 'read.snps.long' function which\n")
    cat("even at the time was listed as deprecated. It has been detected as removed from snpStats\n")
    cat("so now redirecting to the newer (slower) function 'read.long' to work in its place.\n")
    read.snps.long <- function(files,sample.id,snp.id,fields,codes=NA,no.call="",threshold=.9,
                               diploid = NULL, lower = TRUE, sep = " ", comment = "#", skip = 0,
                               simplify = c(FALSE,FALSE),
                               verbose = FALSE, in.order=TRUE, every = 1000,...) {
      # suck up redundant variables and redirect relevant ones to 'read.long'
      out <- read.long(file=files, samples=sample.id, snps=snp.id,
                       fields = fields,split = "\t| +", gcodes=codes,
                       no.call = no.call, threshold = threshold)
      return(out)
    }
    was.there <- F
  } else {
    was.there <- T
    read.snps.long <- ""
  }
  return(list(was.there,read.snps.long))
}


import.snp.matrix.list <- function(snp.list,dir,data=NULL,samples=NULL,
                                   snp.fields=NULL,HD=F,multi=F,n.cores=1,verbose=T)                       
{
 # Reads in the SNP data for the whole dataset, with a separate snpMatrix entry per cohort
 # setting HD=T uses a lot less RAM, can be used for very large datasets (but slow)
 # option multi allows use of parallel processing, but the main limited is hard disk
 # speed, so for now this is likely to decrease performance if anything.
 num.markers <- length(snp.list)
 # generate locations of the data files:
 dir <- validate.dir.for(dir,c("raw","col","cr","lrr.dat","ano"),warn=T)
 must.use.package("snpStats",bioC=T)
 redirect.fn <- check.snp.read.fn()
 if(!redirect.fn[[1]]) { read.snps.long <- redirect.fn[[2]] }
 # dir$raw should be the location of the raw data files specified by
 ### locations of files specifying raw data (separate file for zipped/unzipped data)
 if(is.null(data)) {
   #DEP ?
   if(auto.mode()) {
     file.info <- get.file.specs(dir)
     fns <- file.info[,1]
     cat("",length(fns),"data files found\n")
     data <- cat.path(dir$raw,fns,must.exist=T)
   }
   if(is.null(data)) { stop("Error: could not find a file for parameter 'data'") }
 } else {
   data <- find.file(data,dir$raw,dir) # assume genotype data would be here, else abs. loc needed
 }
 if(is.null(samples)) 
 { 
   if(auto.mode()) { subs.list <- get.subIDs(dir,"list") } else { stop("Error: parameter 'samples' was blank") }
 } else {
   if(length(samples)>10) {
     subs.list <- samples
   } else {
     subs.list <- sapply(cat.path(dir$ano,samples,must.exist=T),readLines)
   }
 } 
 num.filz <- length(data)
 # define names and column locations of fields in the data file
 field.list <- list() 
 if(is.null(snp.fields)) {
   # no field.list object passed in so try to read from file.spec.txt
   # might want to make this whole bit DEP #DEP
   file.info <- get.file.specs(dir)
   if(!is.null(file.info)) {
     #defined using columns of 'file.spec.txt'
     for (cc in 1:num.filz) {
       field.list[[cc]] <- c(sample = file.info$SAMP[cc], snp = file.info$SNP[cc], 
                          allele1 = file.info$A1[cc], allele2 = file.info$A2[cc]) }
   } else {
    cat("Error: import requires setup file 'file.spec.txt or else 'field.list': a list for each file,\n")
     cat("with elements: sample,snp,allele1,allele2; giving columns numbers for each in the datafile\n")
     stop()
   }
 } else {
   # allow possibility only 1 row of fields to specify for all
   if(nrow(snp.fields)>=num.filz) { ind <- 1:num.filz } else { ind <- rep(1,num.filz) }
   for (cc in 1:num.filz) {
     ee <- ind[cc]
     field.list[[cc]] <- c(sample = snp.fields$SAMP[ee], snp = snp.fields$SNP[ee], 
                           allele1 = snp.fields$A1[ee], allele2 = snp.fields$A2[ee]) 
   }
 }
 ########################
 len.tab <- get.file.lens(dir,write=T)
 cat("\nTime estimates for importing SNP data:\n")
 print(time.est.snps(len.tab),quote=F)
 file.lens.nh <- (len.tab[,1])  #lengths only 
 # get header lengths (may not be uniform)
 header.lens <- get.hdr.lens(data)
 file.lens.sub <- (file.lens.nh)/num.markers
 if(any(!check.data.size(num.markers,file.lens.sub))) { stop("Proposed file too large") }
 cat("Files contain:",paste(file.lens.sub,collapse=","),"samples, respectively\n")
 # go to /CALLRATE directory
 old.dir <- getwd(); setwd(dir$cr)
 ######
 #############id.file.names <- paste(data,".ids",sep="")
 if(!is.list(subs.list)) { subs.list <- list(subs.list) } # force list
 if(length(subs.list)!=length(data)) { stop("Must be same number of datafiles as lists of IDs") }
 # Read each set of SNPs into a SnpMatrix object stored in a list element
 snpMatLst <- list()
 #if HD=T only store list of locations of saved objects (in case of memory limitation)i
 cat("\n")
 #else list contains the objects themselves
 if(!multi) {
   options(warn = -1)
   for (tt in 1:num.filz) {
     kk <- proc.time()
     if(!verbose) { cat(" importing long format snp data for file",tt,"...") }
     snpMat <- read.snps.long(files = data[tt], sample.id = subs.list[[tt]],    
                              snp.id = snp.list, diploid = NULL, 
                              fields = field.list[[tt]], 
                              codes = "nucleotide", sep = "\t", comment = "#", 
                              skip = header.lens[tt], simplify = c(FALSE,FALSE),
                              verbose = verbose, in.order = TRUE, every = num.markers)
     ## used to save only in HD-mode but now saving either way ##
     if(!verbose) { cat("..done\n") }
     ofn <- paste("snpMat",tt,".RData",sep="")
     save(snpMat,file=ofn)
     if(HD) {       
       snpMatLst[[tt]] <- ofn
     } else {
       snpMatLst[[tt]] <- snpMat ; snpMat <- NULL
     }
     jj <- proc.time()
     cat(paste("file",tt,"took",round((jj[3]-kk[3])/60),"minutes\n\n"))
   }
   options(warn = 0)
 } else {
   n.cores <- min(num.filz,n.cores); rdz <- vector("list",n.cores)
   must.use.package("multicore"); c.u <- 0; snpMatLst <- vector("list",num.filz)
   cat(" reading",num.filz,"long format SNP files using",n.cores,"cores in parallel...\n")
   for (tt in 1:num.filz) {
     c.u <- c.u + 1
     rdz[[tt]] <- multicore::parallel(read.snps.long(files = data[tt], sample.id = subs.list[[tt]],    
                                                snp.id = snp.list, diploid = NULL, 
                                                fields = field.list[[tt]], 
                                                codes = "nucleotide", sep = "\t", comment = "#", 
                                                skip = header.lens[tt], simplify = c(FALSE,FALSE),
                                                verbose = F, in.order = TRUE, every = num.markers))
     if(c.u>=length(rdz)) { 
       snpMatLst[(tt-length(rdz)+1):tt] <- collect(rdz) ; rdz <- vector("list",min(n.cores,num.filz-tt)) ; c.u <- 0 
     }
   }	
   #snpMatLst <- collect(snpMatLst)
 }
 setwd(old.dir) #put current directory back to where it was
 return(snpMatLst)
}


draw.density.plots <- function(fn=NULL,sample.info,snp.info,samp.result=NULL,snp.result=NULL,
                            callrate.samp.thr=.95, callrate.snp.thr=.95, anot=F, par=NULL) 
{
 # draw density plots for SNP/sample-wise callrate evaluation
 cat(" generating density plots\n")
 # set text parameters for samples (1) and snps (2); i,j = both
 par.list <- list(a1=1,b1=140,x1=28,y1=14,a2=27,b2=144,x2=32,y2=14,i=0.17,j=0.3)
 # any parameters entered as a named list in par will override these defaults^
 if(!is.null(par)) { par.list <- update.list.with.list(par.list,par) }
 if(!is.null(sample.info$call.rate) & !is.null(snp.info$call.rate))
 {
   if(is.character(fn)) { pdf(fn) }
   par(mfrow=c(2,1))
   #samples
   plot(density(sample.info$call.rate),xlim=c(callrate.samp.thr,1),main="A. Sample-wise call rate distribution")
   if(anot & length(dim(samp.result))>0) {
     if(all(dim(samp.result)==c(8,3))) {
       attach(par.list)
       tt <- paste(as.matrix(samp.result,byrow=T))
       tt <- c(colnames(samp.result),tt)
       tt <- tt[c(1,4:11,2,12:19,3,20:27)]
       spn <- 1-callrate.samp.thr
       mdpt <- callrate.samp.thr+ spn/2
       if(median(sample.info$call.rate,na.rm=T)<mdpt) {
         ofs <- callrate.samp.thr + spn*.6 
       } else {
         ofs <- callrate.samp.thr
       }
       text(rep(c(0,i,j)*spn,each=9)+ofs,a1*rep(seq(b1,x1,by=-y1),times=3),labels=tt,cex=0.6,pos=4)
       detach(par.list)
     }
   }
   #snps
   plot(density(snp.info$call.rate),xlim=c(callrate.snp.thr,1),main="B. SNP-wise call rate distribution")
   if(anot & length(dim(snp.result))>0) {
     if(all(dim(snp.result)==c(8,3))) {  
       attach(par.list)
       tt<-paste(as.matrix(snp.result,byrow=T))
       tt <- c(colnames(snp.result),tt)
       tt <- tt[c(1,4:11,2,12:19,3,20:27)]
       spn <- 1-callrate.snp.thr
       mdpt <- callrate.snp.thr+ spn/2
       if(median(snp.info$call.rate,na.rm=T)<mdpt) {
         ofs <- callrate.snp.thr + spn*.6 
       } else {
         ofs <- callrate.snp.thr
       }
       text(rep(c(0,i,j)*spn,each=9)+ofs,a2*rep(seq(b2,x2,by=-y2),times=3),labels=tt,cex=0.6,pos=4)
       detach(par.list)
     }
   }
   if(is.character(fn)) { dev.off(); cat(paste("~wrote density figure to file:",fn,"\n")) }
 } else {
   warning("invalid snp.info or sample.info objects, could not produce graph")
 }
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
   s.info <- column.salvage(cs,"call.rate",c("Call.rate","callrate"),ignore.case=T)
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


update.list.with.list <- function(orig.list, add.list)
{
 # idea is start with a full list, and if any elements with
 # the same name are in add.list, and they are the same type,
 # then the original list entry will be updated with the add.list
 # entry. the idea is to facilitate easy control of a large number
 # of parameters (which have defaults) while only having to specify
 # those to change, not all of them
 if(is(add.list)[1]=="list" & is(orig.list)[1]=="list")
 {
   for (cc in 1:length(add.list)) {
     nxt.nm <- names(add.list)[cc]
     if(nxt.nm %in% names(orig.list))
     { 
       if(is(add.list[[cc]])[1]==is(orig.list[[nxt.nm]])[1]) { 
         orig.list[[nxt.nm]] <- add.list[[cc]] 
       } 
     }
   }
 }
 return(orig.list)
}


doSampQC <- function(dir, subIDs.actual, plink=T, callrate.samp.thr=.95, het.lo=0.1,het.hi=0.4,snpMatLst=NULL, n.cores=1, proc=1) 
{
 # Do sample callrate quality control on a snpMatrix object (or just use plink files
 # if the QC has already been performed in plink)
 ## get sample call rates ##
 dir <- validate.dir.for(dir,c("ano","cr.plk"),warn=F)
 if(is(snpMatLst[[1]])[1]=="SnpMatrix") {
   group.nums <- rep(c(1:length(snpMatLst)),sapply(snpMatLst,dim)[1,])
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
     group.nums <- rep(c(1:length(snpMatLst)),times=lsp[1,])
   }
 }
 if(length(group.nums)!=length(subIDs.actual)) { stop(paste("Error: group was length",length(group.nums),
 	"but IDs had length",length(subIDs.actual),"so couldn't construct a data.frame"))}
 sample.info <- data.frame(grp=group.nums,row.names=subIDs.actual,QCfail=rep(0,times=length(group.nums)))
 if (plink | is.null(snpMatLst)){
   if(!plink) { cat(" snpMatLst is NULL but plink=FALSE, trying to see whether plink QC data is available\n") }
   if (!"snpdataout.irem" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink samples-removed-by-QC file: snpdataout.irem, in",dir$cr.plk)) }
   if (!"snpdataout.imiss" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink samples-QC output file: snpdataout.imiss, in",dir$cr.plk)) }
   call.rate.excl.grp <- paste(read.table(paste(dir$cr.plk,"snpdataout.irem",sep=""),header=F)[,2])
   sample.info[call.rate.excl.grp,"QCfail"] <- proc
   cat(paste(" plink sample QC removed",length(call.rate.excl.grp),"samples\n"))
   imisser <- read.table(paste(dir$cr.plk,"snpdataout.imiss",sep=""),header=T)
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
 } else {
   sample.qc <- list.rowsummary(snpMatLst,dir=dir,n.cores=n.cores)
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


get.grp.level.snp.stats <- function(snpMatLst,snp.info,sample.info,n.cores=1,plot.fn=NULL) {
  ### do some snp-qc on each group separately and look for between-grp differences
  if(length(snpMatLst)==length(unique(sample.info$grp))) {
    if(length(snpMatLst)>1) {
      # if more than one group do snp qc separately for each to test for differences between grps
      if(n.cores>1) { 
        # get SNP-qc stats summary by grp
        snp.qc.grpwise <- multicore::mclapply(snpMatLst,colSummary,filt=rownames(snp.info),mc.cores=n.cores) 
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


doSnpQC <- function(dir, plink=T, n.cores=1,
                   callrate.snp.thr=.95, hwe.thr=(10^-7), grp.hwe.z.thr=4, grp.cr.thr=.001, group.miss=F,
                   snpMatLst=NULL, autosomes.only=T, subIDs.actual, snp.info, sample.info=NULL, proc=1) 
{
 # Do SNP callrate/HWE quality control on a snpMatrix object (or just use plink files
 # if the QC has already been performed in plink)
 dir <- validate.dir.for(dir,c("ano","cr.plk","cr"),warn=F)
 ## get snp call rates ##
 if(is(snpMatLst[[1]])[1]=="SnpMatrix") { HD <- F } else { HD <- T }
 if (plink | is.null(snpMatLst)){
   if (!"snpdataout.lmiss" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink snp-QC output file: snpdataout.lmiss, in",dir$cr.plk)) }
   if (!"snpdataout.hwe" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink snp-HWE output file: snpdataout.hwe, in",dir$cr.plk)) }
   lmisser <- read.table(paste(dir$cr.plk,"snpdataout.lmiss",sep=""),header=T)
   snp.info[["call.rate"]] <- 1-(lmisser[match(rownames(snp.info),lmisser$SNP),"F_MISS"])
   hweer <- read.table(paste(dir$cr.plk,"snpdataout.hwe",sep=""),header=T)
   snp.info[["P.hwe"]] <- hweer[match(rownames(snp.info),hweer$SNP),"P"]
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
   exp.snp <- length(chr.nums.of.info(select.autosomes(snp.info))) #,22)
   if(length(snpMatLst)!=exp.snp) {
     snpMatLst <- convert.smp.to.chr22(snpMatLst,snp.info=snp.info,dir=dir,n.cores=n.cores) }
   snp.qc <- list.colsummary(snpMatLst,dir=dir,n.cores=n.cores)
   snp.info[["call.rate"]] <- snp.qc[rownames(snp.info),"Call.rate"]
   snp.info[["maf"]] <- snp.qc[rownames(snp.info),"MAF"]
   snp.info[["het"]] <- snp.qc[rownames(snp.info),"P.AB"]
   snp.info[["BAF"]] <- ((.5*snp.qc[rownames(snp.info),"P.AB"]) + snp.qc[rownames(snp.info),"P.BB"])
   cond <- snp.info$maf<.0005
   cond[is.na(cond)] <- F
   mmz <- rownames(snp.info)[cond]
   ofn <- cat.path(dir$cr,"monomorphic.txt"); cat("~wrote",length(mmz),"suspected monomorphic snps to:\n ",ofn,"\n")
   mafz <- cbind(mmz,round(snp.info[["maf"]][cond],7)); colnames(mafz) <- c("marker.id","minor.allele.frequency")
   write.table(mafz,file=ofn,quote=F,row.names=F,sep="\t")
   zz <- snp.qc[rownames(snp.info),"z.HWE"]
   pp <- (1-pnorm(as.numeric(abs(zz))))*2  #two-tailed p value
   snp.info[["P.hwe"]] <- pp
   snp.info[["Z.hwe"]] <- zz
 }
 out.list <- apply.snp.thresholds(snp.info,callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,
                                  grp.cr.thr=grp.cr.thr,grp.hwe.z.thr=grp.hwe.z.thr,proc=proc)
 return(out.list)
}



get.hdr.lens <- function(fnz,max.feas.hdr=500,firstrowhdr=T,sep.chr="\t") 
{
 # get the lengths of headers in genome studio files
 # assuming the header must be less than 'feas.hdr'/500 lines long.
 # starts from line 500 counts the number of tabs/separators and 
 # looks for when the number first changes. if it doesn't assumes 1st row only is header
 # if firstrowhdr = T, else assume no headers at all
 header.lens <- 1
 num.tabs <- function(str,ch="\t") {  return(-1+sapply(strsplit(str,ch,fixed=T),length)) }

 for (cc in 1:length(fnz))
 {
   tabs.per.line <- num.tabs(readLines(fnz[cc],n=max.feas.hdr),ch=sep.chr)
   full.n.tabs <- tail(tabs.per.line,1)
   fl <- 1; notfnd <- T
   while(notfnd) {
     if(tabs.per.line[fl]!=full.n.tabs) { fl <- fl + 1 } else { notfnd <- F }
   }
   if(firstrowhdr) { header.lens[cc] <- fl } else { header.lens[cc] <- fl-1 }
 }
 return(header.lens)
}


load.all.libs <- function(big=c("bigmemory","biganalytics"),
                          other=c("multicore","lattice","compiler","NCmisc"),...) {
  # loads most of the libraries required by plumbcnv in a fairly discreet manner
  # note that lattice uses multicore function 'parallel' so should be loaded second
  # because plumbCNV always refers to multicore::parallel to prevent confusion
  ih <- F
  if(length(big)>0) { 
    if("bigalgebra" %in% big) {
      try({ if(!require(bigalgebra)) { ih <- T } })
      if(ih) {
        try( { if(!big.algebra.install.help())
        {
          warning("bigalgebra not installed, some functions may be slowed")
        } })
      }
      big <- big[-which(big %in% "bigalgebra")]
    } 
    must.use.package(big,quiet=T)
  }
  must.use.package(other,quiet=T)
  preload.bioC(...)  # load bioC packages without annoying warnings, etc
}


preload.bioC <- function(bio.packs=c("IRanges","BiocGenerics","Biobase","GenomicRanges","genoset"),more.bio=c()) 
{
  ## load bioconductor packages without annoying package loading text ##
  bio.packs <- c(bio.packs,more.bio) # allows adding to, rather than replacing default set
  if(!all(paste("package:",bio.packs,sep="") %in% paste(search()))) {
    cat(" silently loading bioconductor packages:",paste(bio.packs,collapse=", "),"... ")
    suppressWarnings(suppressMessages(must.use.package(bio.packs,T,quietly=T)))
    cat("done\n")
  }
}


snp.mat.list.type <- function(snpMatLst,fail=T)
{
 # plumbCNV snp-qc can be run on a list of snpMatrix objects in memory ('memory') or
 # on a list of RData files on 'disk' to save RAM, this fn detects which type is in use
 typz <- sapply(snpMatLst,is)[1,]
 if(all(typz=="SnpMatrix"))
 { HDt <- "memory" } else {
   if (all(typz=="character")) { HDt <- "disk" } else {
     HDt <- "error"
     if(fail) {  stop("Error: not a valid snpMatLst!") 
     } else { warning("couldn't classify type as was not a valid snpMatLst")  }
   }
 }
 return(HDt)
}


get.snpMat.spec <- function(snpMatLst,fail=T,dir=NULL)
{
 # get dimensions of snpMatLst (SnpMatrix list) regardless
 # of whether it's a set of SnpMatrix objects or list of file locations
 HD <- switch(snp.mat.list.type(snpMatLst,fail),memory=F,disk=T,error=NULL)
 if(is.null(dir)) { dir <- getwd() ; warning("no directory passed to function") }
 if(HD) {
   n.grp <- length(snpMatLst)
   list.spec <- matrix(integer(),ncol=n.grp,nrow=2)
   for (cc in 1:n.grp)
   {
     TF <- is.file(paste(snpMatLst[[cc]]),dir)
     if(TF) { 
       snpMat <- get(paste(load(find.file(paste(snpMatLst[[cc]]),dir)))) 
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



correct.gc.old <- function(next.cols,col.sel,gc,gc_log)
{
  # apply PC correction for a single SNP, allowing for missing data.
  bad1 <- which(is.na(next.cols))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  next.cols[sel] <- lm(next.cols ~ gc + gc_log,na.action="na.exclude")$residuals
  return(next.cols)
}


correct.gc.mat.old <- function(next.cols,gc,gc_log)
{
  # matrix version of correct.gc (used to GC-correct one SNP at a time)
  col.sel <- 1:ncol(next.cols)
  for (dd in 1:nrow(next.cols)) {
    next.cols[dd,] <- correct.gc(next.cols[dd,],col.sel,gc,gc_log) 
  }  
  return(next.cols)
}


remove.samp.from.plink <- function(samps,plink,why="TMcnvs",dir="")
{
	# remove specific sample from a plink file - eg because had too many cnvs
 ext <- c("cnv","fam","cnv.map")
 what <- c("cnvs","samples","")
 for (dd in 1:3) {
   fnn <- cat.path(dir="",plink,ext=ext[dd])
   inp <- readLines(fnn)
   to.rmv <- NULL
   if(length(samps)>0) {
     for (cc in 1:length(samps))
     {
       to.rmv <- c(to.rmv,grep(samps[cc],inp))
     }
   }
   to.rmv <- unique(to.rmv)
   if (length(to.rmv)>0)
   {
     if(any(dir!="")) {
       # if dir is passed then make a sample exclusion file
       dir <- validate.dir.for(dir,"excl")    
       fn <- cat.path(dir$excl,"CNVQC.txt")
       if(file.exists(fn)) { ex <- readLines(fn) } else { ex <- NULL }
       ex <- c(ex,inp[to.rmv])
       writeLines(ex,con=fn)
       cat(" appended CNV-QC failing samples to:",fn,"\n")
     }
     out <- inp[-to.rmv]
     cat(paste("",length(to.rmv),what[dd],"removed from:\n ",fnn,"\n"))
   } else {
     out <- inp
     cat(paste(" no samples removed from:\n ",fnn,"\n"))
   }
   writeLines(out,con=cat.path(dir="",plink,suf=why,ext=ext[dd]))
 }
}


lrr.stats.tab <- function(stats.table,nSD=3)
{
  ## CALCULATE LRR STATS ##
  ## based on a table of sample-wise LRR stats (eg, mean, stdev, dlrs, etc)
  # calculate distribution indices for each LRR-sample-statistic for the whole cohort
  if(!is.numeric(nSD)) { nSD <- 3 }; nSD <- abs(nSD) # ensure positive and numeric
	statz <- colnames(stats.table)
  stats.table <- as.data.frame(stats.table)  #ensure type is dataframe (eg. not matrix)
	TableRowLabels <- boundary.stats(nSD=nSD) # in function user can look at
	nrowz <- length(TableRowLabels)
	ncolz <- length(statz)
	s.tab <- matrix(numeric(),ncol=ncolz,nrow=nrowz)
	
	for (cc in 1:ncolz)
	{
	 nxt.type <- statz[cc]
	 next.stat <- stats.table[[paste(nxt.type)]]
	 sl <- summary(next.stat)
	 IQR <- (sl[5]-sl[2])
	 pcts <- pctile(next.stat)
	 StDev <- sd(next.stat)
	 s.tab[,cc] <- c(sl[4],StDev,sl[c(1:3,5:6)],sl[2]-1.5*IQR,sl[5]+1.5*IQR,sl[4]+c(-2,2,-nSD,nSD)*StDev,pcts)
	}	
	rownames(s.tab) <- TableRowLabels
	colnames(s.tab) <- statz
	return(s.tab)
}


plot.extreme.samples <- function(rez,stat.table,bigMat2,snp.info,dir,pref="",
                                        scl=10^6,ap="",autosomes=F,ucsc="hg18")
{
  # plot the most extreme sample for each combination of QC failure types
  dir <- validate.dir.for(dir,c("ind","big"),warn=F)
  venn.lists <- get.fail.type.list(rez,stat.table[,"Mean"]) 
  ex.id <- get.extreme.examples(stat.table,rez,venn.lists)
  cat("\nPlotting extreme samples for combinations of: Mean vs GC vs DLRS\n")
  ex.id <- ex.id[!is.na(ex.id)] ; ex.id <- ex.id[sapply(ex.id,length)!=0] # remove any categories with missing or length==0
	loopz <- length(ex.id)
  warnz <- NULL
	ofn <- character(loopz)
  if(is(snp.info)[1]!="RangedData") { 
    warning("invalid snp.info object, plot skipped"); return(NULL)
  }
  if(!all(c("gindx") %in% colnames(snp.info))) {
    warning("snp.info object did not contain 'gindx' column, plot skipped")
    return(NULL)
  }
  if(!"color" %in% colnames(snp.info)) {
    snp.info <- add.color.to.snp.info(snp.info,scheme="unique")
  }
  if (autosomes) { snp.info <- select.autosomes(snp.info) }
  uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
  XX <- (snp.info$gindx)/scl
  chrLens <- get.chr.lens(dir,ucsc=ucsc)[chr.set]
  chrStarts <- c(0,cumsum((chrLens)))[1:n.chr] # genome position of chromosome starts
	for (cc in 1:loopz)
	{
		if(!is.na(ex.id[cc])) {
		 loc.col <- match(ex.id[cc],colnames(bigMat2)) 
		 #LRR.dat <- sub.big.matrix(bigMat2, firstCol = loc.col, lastCol = loc.col, backingpath=dir$big)
     LRR.mat <- bigMat2[1:nrow(bigMat2),loc.col]
		 select <- c(1:length(LRR.mat)) ; ccc <- snp.info$color[select]
		 titl <- c(paste("Sample",ex.id[cc]),paste("Extreme LRR",names(ex.id)[cc]))
		 ofn[cc] <- paste(dir$ind,names(ex.id)[cc],ap,"Sample",ex.id[cc],pref,".pdf",sep="")
		 pdf(ofn[cc])
		   par(mar=c(5, 4, 10, 2))
		   plot(XX,LRR.mat,col=ccc,pch=".",cex=1.5,
		        xlab = "Genome Position (Megabases)", ylab="Log R Ratio",
		        main = titl)
		   axis(side=3,at=((chrStarts[1:n.chr]+(chrLens[1:n.chr]/2))/scl), labels=paste(chr.set))
		   mtext ("Chromosome number",side=3,line=2.5,cex=.9)
		 dev.off()
     loop.tracker(cc,loopz)
		} else { 
      #no subjects in this category, skipped..
      warnz <- c(warnz,paste(names(ex.id)[cc]))
      loop.tracker(cc,loopz)
    }
	}
  if(!is.null(warnz)) {
    warning(paste("no subjects found for extreme values of:",paste(warnz,collapse=",")))
  }
	cat(paste("~wrote",length(ofn),"files:\n"))
  cat(dirname(ofn[rev(order(nchar(ofn)))][1]),";\n",sep="")
	cat(paste("  ",rmv.spc(basename(ofn)),"\n",sep=""),sep=""); cat("\n")
  return(NULL)
}	


get.all.samp.fails.old <- function(dir)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl"),warn=F)
  ###src <- dir$excl;   src <- paste(dir$ano,"SAMPLE_EXCLUDE/",sep="")
  fnz <- list.files(dir$excl)
  excl.samps <- character()
  for (cc in 1:length(fnz)) {
    excl.samps <- c(excl.samps,readLines(paste(dir$excl,fnz[cc],sep="")))
  }
  excl.samps <- unique(excl.samps)
  return(excl.samps)
}


get.chr.ab.fail.subset <- function(chrWarns, dir, failerMode="NOTLRR", max.bad.to.plot=22) 
{
 ### SCRIPT TAKES LIST OF SAMPLES WITH CHR ABERRATIONS (using contrasts)
 ##  AND SELECTS A SUBSET TO FLAG/PLOT WHICH CAN BE BASED ON NOT FAILING OTHER CRITERIA
 dir <- validate.dir.for(dir,c("ano","excl"),warn=F)
 n.bad.per.samp <- sort(table(unlist(chrWarns))) # get list of which chrs bad for each sample
 # get previous LRR exclusion lists from files
 fnfn <- cat.path(dir$excl,c("Mean","DLRS","GCWave"),ext="txt")  #"StDev",
 excluded.prev <- character()
 for (cc in 1:length(fnfn)) {
   excluded.prev <- c(excluded.prev,readLines(fnfn[cc]))
 }
 excluded.prev <- unique(excluded.prev)

 ## which samples to use with respect to those excluded already on other stats
 if(failerMode==("NOTLRR")) {
   checklist <- n.bad.per.samp[which(!names(n.bad.per.samp) %in% excluded.prev)]
 }
 if(failerMode==("ONLYLRR")) {
   checklist <- n.bad.per.samp[which(names(n.bad.per.samp) %in% excluded.prev)]
 }
 if(failerMode==("ALL")) { checklist <- n.bad.per.samp }

 # select sample subset with only up to a max number of bad chromosomes
 # eg, might only be interested in looking at those with 1 bad chromosome
 badcheckz <- names(checklist)[checklist<=max.bad.to.plot]

 # decode for each sample which chromosomes are bad and create labels for graphs
 chrLab <- character(length(badcheckz))
 chrN <- list()
 for (mm in 1:length(badcheckz)){
   failstr <- names(unlist(chrWarns))[which(unlist(chrWarns) %in% badcheckz[mm])]
   extr.str <- sapply(strsplit(failstr,"_",fixed=T),"[",1)
   chrN[[mm]] <- gsub("chr","",extr.str)
   if(length(failstr)==1)
   { chrLab[mm] <- paste("Chromosome",chrN[[mm]]) } else {
     chrLab[mm] <- paste("Chromosomes",paste(chrN[[mm]],collapse=","))
   }
 }
 names(chrN) <- badcheckz
 out.list <- list(badcheckz,chrN,chrLab)
 names(out.list) <- c("badcheckz","chrN","chrLab")
 return(out.list)
}


insert.fake.median <- function(pass.stats)
{
	# a hack for now to avoid fixing code so Median can be left out of tables
	cnmz <- colnames(pass.stats)
	if(!"Median" %in% cnmz)
	{
		cnmz <- cnmz[c(1,2,1,3)] ;cnmz[3] <- "Median"
		pass.stats <- pass.stats[,c(1,2,1,3)]
		colnames(pass.stats) <- cnmz	
	}
	return(pass.stats)
}


add.chr.top.axis <- function(select, badChrN, x.co, nC=22, sub=T) 
{
 # on chromosome graphs there is an option in plumbCNV to have the
 # chromosomes on the top axis, this script manages this allowing for
 # different levels of zoom which may cause overlaps. manages it better
 # than the default. also can change a target chr label to red colour.
 badChrz <- ((1:nC)[select] %in% badChrN )
 ## labels get messed as these too close together- so hack to fix
 if(any((c("18","19","21","22")) %in% badChrN ) & nC>20)
 {
   if("19" %in% badChrN )
   { badChrz[20] <- T } else {
     if("21" %in% badChrN )
     { badChrz[c(20,22)] <- NA
     } else { 
       if("22" %in% badChrN )
       { badChrz[21] <- T } else { badChrz[19] <- T }
     }
   }
 }
 axis(side=3,at=x.co[select][!badChrz], labels=paste(1:nC)[select][!badChrz],
      col.axis="black")
 if(length(paste(1:nC)[select][badChrz])>0) {
   axis(side=3,at=x.co[select][badChrz], labels=paste(1:nC)[select][badChrz],
      col.axis="red") }
 if (sub) { mtext ("Chromosome number",side=3,line=2.5,cex=.9) }
 return(list(select,badChrz,badChrN))
}


bw.plot.lrr <- function (samp.chr.means, SX=NULL, snp.info=NULL, samp.chr.sds=NULL, num.labels=T,
                        CI.Z=1.96, grpAv=NA, grpLoCI=NA, grpHiCI=NA, scl=10^6, ucsc="hg18",dir="",...) 
{
 # Black and White Type Plot
 # Alternate graph without the full LRR data
 # This plots the means of each chromosome and their confidence intervals
 # also overlays the mean and std deviation bands of the reference group
 # of all good samples
 # CI.Z = confidence interval z score
 if(is(snp.info)[1]=="RangedData") { 
   uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
   snp.info <- select.autosomes(snp.info)
   chr.set <- chr.nums.of.info(snp.info); nC <- length(chr.set)
 } else { nC <- 22 ; chr.set <- 1:nC }
 chrLens <- get.chr.lens(dir,ucsc=ucsc)[chr.set]
 chrStarts <- c(0,cumsum((chrLens)))[1:nC] # genome position of chromosome starts
 if(is.null(SX)) {  SX <- ((chrStarts[1:nC]+(chrLens[1:nC]/2))/scl) }
 YY <- samp.chr.means
 if(length(YY)!=nC) { warning("number of chromosomes in data doesn't match snp.info") }
 chr.select <- 1:nC
 plot(SX, grpAv, cex=1.5,
      xlab = "Genome Position (Megabases)", ylab="Log R Ratio", 
      ,type="l", lty="dashed", col="black", ...)
 if (!is.na(grpLoCI[1]) & !is.na(grpHiCI[1])) {
   lines(SX,grpLoCI,lty="dotted",col="black")
   lines(SX,grpHiCI,lty="dotted",col="black")
 }
 lines(SX,YY) # subject means plot (main data)
 # confidence intervals based on specified confidence level
 if(!is.null(samp.chr.sds) & !is.null(snp.info)) {
   nECF <- sapply(snp.info,nrow) 
   SS <- CI.Z*(samp.chr.sds/sqrt(nECF))
   arrows(SX,YY,SX,YY+SS,angle=90,length=0.05)
   arrows(SX,YY,SX,YY-SS,angle=90,length=0.05)
 }
 if(num.labels) {
   text(x=SX,y=YY+(.1),labels=paste(round(YY,3)),srt=90,cex=.75) }
 return(SX)
}


plot.chr.ab.samples <- function (dir,bigMat2,chr.stat,chr.ab.samples.obj,snp.info,
                                 colPlot=F,failerMode="NOTLRR",lob=2,hib=2.5,pctile.bound=.01,
                                 max.bad.to.plot=22,medSmooth=F,ratio=10,
                                 rngOn=F,c.rel=T,c.pre.n=2, c.post.n=2, pref=pref) {
  ## PLOT IN ONE OF TWO WAYS (BW or Color)
  ## plot all bad samples in 'checklist'
  dir <- validate.dir.for(dir,c("big","ano","qc.cab"))
  ## from output of function, 'get.chr.ab.fail.subset()'
  if(!is.list(chr.ab.samples.obj) | length(chr.ab.samples.obj)!=3) {
    cat("Invalid chr.ab.samples.obj passed in. Please use get.chr.ab.fail.subset() function to get a valid list\n")
  }
  badcheckz=chr.ab.samples.obj$badcheckz; chrLab=chr.ab.samples.obj$chrLab; chrN=chr.ab.samples.obj$chrN;
  if (colPlot) { cP <- "Col" } else { cP <- "BW" }
  badPlotFileName <- paste("SampleChrAb_",cP,"_",failerMode,"_lo",lob,"_hi",hib,"_pc",
                           pctile.bound,"_n",1,"_",max.bad.to.plot,".pdf",sep="")
  if(!is(snp.info)[1]=="RangedData") { stop("Error: not a valid snp.info object") }
  snp.info <- snp.info[(rownames(snp.info) %in% rownames(bigMat2)),]
  ####scl <- 1000000 # graph scale factor (e.g, megabases=1 million)
  CIz <- hib  # confidence interval z score: 3.29 = 99.9%, 2.575 = 99%, 1.96 = 95%
  chr.mean <- chr.stat$chr.mean; chr.sd <- chr.stat$chr.sd
  # get some stats (mean/CIs) for the group average (minus bad samples)
  RR <- colMeans(chr.mean[-match(badcheckz,rownames(chr.mean)),])
  RRSD <- apply(chr.mean[-match(badcheckz,rownames(chr.mean)),],2,sd,na.rm=T)
  RRL <- RR-(lob*RRSD)
  RRU <- RR+(lob*RRSD)
  # set default to select all chromosomes for plotting
  chr.set <- chr.nums.of.info(snp.info); nC <- length(chr.set)
  if(length(RR)!=nC) { warning("number of chromosomes in data doesn't match snp.info") }
  chr.select <- 1:nC
  if(colPlot) {
    if(!"gindx" %in% colnames(snp.info)) { 
      colPlot <- F
      warning("color plots only work with the extended form of snp.info, with a 'gindx' column")
    }
    if(!"color" %in% colnames(snp.info)) {
      snp.info <- add.color.to.snp.info(snp.info,scheme="unique")
    }
  }
  # if smoothing option in use then adjust y axis limits accordingly
  if (medSmooth) { limz <- c(-.5,.5) } else { limz <- c(-2,2) }
  if (colPlot) { cat("plotting samples in colour (slow)...\n") }
  ofn <- cat.path(dir$qc.cab,badPlotFileName,pref=pref); num.plots <- length(badcheckz)
  pdf(ofn)
  ## loop through each bad sample plotting the desired plots
  if(num.plots>1) { cat("Plotting",num.plots,"samples failing only on chromosomal abberations:\n") }
  for (cc in 1:num.plots)
  {
    titl <- c(paste("Sample",badcheckz[cc]),paste("Aberrant",chrLab[cc]))
    par(mar=c(5, 4, 10, 2))
    # start plotting 
    if (colPlot) {
      # full LRR data for chromosomes (full colour)
      out.list <- col.plot.lrr(ID=badcheckz[cc], bigMat=bigMat2, snp.info=snp.info,
                               centre.chr=chrN[[badcheckz[cc]]], plotAdj=rngOn,
                               samp.chr.means=chr.mean[badcheckz[cc],],ratio=ratio,
                               whole.genome=!crel, c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=medSmooth,
                               ylim=limz, dir=dir,
                               main = titl, xlab = "Genome Position (Megabases)", ylab="Log-R Ratio")
      x.coords <- out.list$x.coords; chr.select <- out.list$chr.select
      legend("topleft",legend=c("LRR-Mean by chromosome",
                                "Raw LRR data [coloured by chromosome]"),cex=.7,
             lty=c("solid",NA),pch=c(NA,"."),col=c("black","grey"),bty="n")
    } else {
      # just the means and CIs for each chromosome, with reference
      x.coords <- bw.plot.lrr(samp.chr.means=chr.mean[badcheckz[cc],], 
                              snp.info=snp.info, samp.chr.sds=chr.sd[badcheckz[cc],], 
                              CI.Z=CIz, grpAv=RR, grpLoCI=RRL, grpHiCI=RRU, ylim=c(-.5,.5), main = titl) 
      legend("topleft",legend=c("LRR-Mean per chromosome [99.9% C.I. bars]",
                              "Average across all good samples [dotted lines: +/-2.5SD]"),cex=.7,
                              lty=c("solid","dashed"),col=c("black","grey"),bty="n")
    }
    # add chromosome labels with bad ones coloured red
    null.result <- add.chr.top.axis(chr.select, badChrN= chrN[[badcheckz[cc]]], x.coords, nC=nC)
    if(num.plots>1) { loop.tracker(cc,num.plots) }
  }
  dev.off()
  cat(paste("~wrote file:",ofn,"\n"))
}


stdev <- function(x,na.rm=T) { sd(x,na.rm=na.rm) }


choose.best.sort.file <- function(search.dir,ref.list)
{
 # hopefully sort directory only contains 1 valid file
 # if there are multiple files, no files, or only invalid files with no matches to ref.list
 # (i.e, to the current big.matrix) this script will try to extract/create the best possible list
 file.list <- paste(search.dir,list.files(search.dir),sep="")
 if(length(grep("SNP_SORT",search.dir,fixed=T))>length(grep("SAMPLE_SORT",search.dir,fixed=T)))
 { guess.mode.is.snp <- T } else { guess.mode.is.snp <- F }
 goodm <- F
 if(length(file.list)>0)
 {
   if(length(file.list)>1) {
     cat(" multiple files found in sort annotation directory - should only be 1\n")
     all.snp <- list(); lens <- integer()
     for (cc in 1:length(file.list))
     {
       all.snp[[cc]] <- (readLines(file.list[cc]))
       lens[cc] <- length(which(ref.list %in% all.snp[[cc]]))
     }
     best.matches <- max(lens,na.rm=T)
     cat("",round(100*(best.matches/length(ref.list)),1),"% of big.matrix indices in sort.list\n")
     if(best.matches>0) {
       out.order <- all.snp[[which(lens==best.matches)[1]]]; goodm <- T
     } 
   } else {
     file1 <- readLines(file.list[1])
     val.len <- length(which(ref.list %in% file1))
     if (val.len>0) {
       cat("",round(100*(val.len/length(ref.list)),1),"% of big.matrix indices in sort.list\n")
       cat(" using",file.list[1],"to sort",c("samples","snps")[1+guess.mode.is.snp],"\n")
       out.order <- file1
       goodm <- T
     } else {
       cat(" Warning: no matches to reference list found in sort file, assuming invalid\n")
     }
   }
 }
 if(!goodm) {
   cat(" no valid files found in sort annotation directory - creating new from 'ref.list'\n")
   if (guess.mode.is.snp) { fnm <- "snpsort.txt" } else { fnm <- "sampsort.txt" }
   ofn <- paste(search.dir,fnm,sep="")
   writeLines(paste(ref.list),con=ofn)
   cat("~wrote file:",ofn,"\n")
   out.order <- ref.list
 }
 return(out.order)
}


check.file.and.subset <- function(fnm, mat, by.row=T, min.size=1000, stop.if.fail=T)
{
 # check that file dir/fnm is a genuine subset of the row/col names (by.row) of matrix 'mat'
 # return subset or if not valid, then return all row/col names of the original 'mat'
 if(by.row)  { ref <- rownames(mat)  } else { ref <- colnames(mat) }
 if(stop.if.fail) { actn <- "aborting script" } else { actn <- "selecting all" }
 if(file.exists(fnm)) {
   out.L <- readLines(fnm)
   #remove headers if required (assumption header unlikely with same number of chars as ids)
   h1 <- out.L[1]
   if(!nchar(h1) %in% nchar(out.L[2:99])) { out.L <- out.L[-1]  ;
                                            cat("\nauto removed header:",h1) }
   if(length(which(out.L %in% ref))<min.size)
   { 
     wrn <- paste(c("column","row")[1+by.row],"subset file",fnm,
                  "did not have at least",min.size,"matches to the reference matrix --> ",actn)
     if(stop.if.fail) { stop(wrn) } else { cat(wrn,"\n"); out.L <- ref }
   }
 } else { 
   wrn <- paste(c("column","row")[1+by.row],"subset file",fnm,"was not found --> ",actn)
   if(stop.if.fail) { stop(wrn) } else { cat(wrn,"\n"); out.L <- ref } 
 }
 return(out.L)
}


sort.exclude.from.annot <- function(bigMat,dir="",dosnp=T,dosamp=T,ordR=T,ordC=T, verb=T)
{
 # sort big matrix according to annotation and exclude those in exclusion files
 dir <- validate.dir.for(dir,c("ano","sort","excl","sort2","excl2"),warn=F)
 cat(" calculating exclusions and order\n")
 samp.names <- colnames(bigMat); snp.names <- rownames(bigMat)
 # order files should only be one file - fix otherwise
 sample.order <- choose.best.sort.file(dir$sort,ref.list=samp.names)
 snp.order   <-  choose.best.sort.file(dir$sort2,ref.list=snp.names)
 ## Get SNP and SAMP exclusion lists ##
 snps.to.cut <- samps.to.cut <- NULL
 if(dosamp) { samps.to.cut <- get.all.samp.fails(dir,verb) }  
 if(dosnp) { snps.to.cut <- get.all.snp.fails(dir,verb) }  
 # use sort/exclusion lists to get reordering vectors
 if(ordR) { to.order.r <- match(snp.order,rownames(bigMat)) } else { to.order.r <- c(1:nrow(bigMat)) }
 if(ordC) { to.order.c <- match(sample.order,colnames(bigMat)) } else { to.order.c <- c(1:ncol(bigMat)) }
 to.order.r <- narm(to.order.r); to.order.c <- narm(to.order.c)
 to.remove.r <- narm(match(snps.to.cut,rownames(bigMat)[to.order.r]))
 to.remove.c <- narm(match(samps.to.cut,colnames(bigMat)[to.order.c]))
 if (length(to.remove.r)>0) { to.order.r <- to.order.r[-to.remove.r] }
 if (length(to.remove.c)>0) { to.order.c <- to.order.c[-to.remove.c] }
 out.list <- list(to.order.r,to.order.c,to.remove.r,to.remove.c,samps.to.cut,snps.to.cut)
 names(out.list) <- c("to.order.r","to.order.c","to.remove.r","to.remove.c","sample.excl","snp.excl")
 return(out.list)
}




big.exclude.sort <- function(des.fn="LRRdescrFile", dir="", deepC=T, tranMode=2, pref="LRRFilt",
                             f.snp="", f.samp="", verbose=T)
{
  # sort and exclude snps/samples from a big.matrix
  must.use.package("bigmemory")
  dir <- validate.dir.for(dir,c("big","ano"),warn=F)
  if (tranMode==1) {  cat("\nFiltering samples/snps with SNP-QC failures:\n")
  } else { cat("\nExcluding samples/snps with SNP-QC failures, sorting SNPs by chr, pos:\n") }
  
  # bigmatrix file names for re-ordered filtered matrix (which is the final output of this script)
  bck.fn.o <- paste(pref,"Sort","bckfile",sep="")
  des.fn.o <- paste(pref,"Sort","descrFile",sep="")
  R.descr <- cat.path(dir$big,des.fn.o,ext=".RData")
  
  bigMat <- getBigMat(des.fn,dir)
  cat(paste(" attached matrix with dims:",paste(dim(bigMat),collapse=","),"\n"))
  # get list of deleting/reordering vectors using annotation files
  if (tranMode==1) {
    #if(f.snp=="") { sif1 <- F } else { sif1 <- T }
    #if(f.samp=="") { sif2 <- F } else { sif2 <- T }
    #snp.L <- check.file.and.subset(fnm=paste(dir$ano,f.snp,sep=""), bigMat, by.row=(rowsAre=="SNP"), stop.if.fail=sif1)
    #samp.L <- check.file.and.subset(fnm=paste(dir$ano,f.samp,sep=""), bigMat, by.row=(rowsAre=="SAMP"), stop.if.fail=sif2)
    trans.list <- select.samp.snp.custom(bigMat,snp=f.snp,samp=f.samp)
    cat(paste(" selected",length(trans.list[[4]]),"listed samples and",length(trans.list[[3]]),"Snps\n"))
  } else {  
    cat(" excluding and ordering listed samples based on annotation directory files; SNP/SAMP\n")
    trans.list <- sort.exclude.from.annot(bigMat,dir=dir,dosnp=T,dosamp=T,ordR=T,ordC=T, verb=verbose)
    ###sort.exclude.from.annot.old(bigMat,snp=T,samp=F,dir=dir$ano,crOff=F)  - original version!
  }
  wrn <- "warning: trans.list was already attached, detaching now..\n"
  while("trans.list" %in% search()) { detach(trans.list); cat(wrn) }
  attach(trans.list)
  if(verbose) {
    cat("\nReordering SNPs and Samples...\n")
    cat("\nINDEXES SUMMARY\n")
    cat(paste(length(to.order.r),"row indexes range is from",min(to.order.r),"to",max(to.order.r),"\n"))
    cat("-->",head(to.order.r),sep=", "); cat ("\n")
    cat(paste(length(to.order.c),"col indexes range is from",min(to.order.c),"to",max(to.order.c),"\n"))
    cat("-->",head(to.order.c),sep=", ")
    cat("\n\n raw big.matrix summary before ordering and exclusion based on SNP-QC:\n\n")
    print.big.matrix(bigMat,"bigMat")
  }
  if(!deepC)
  {
    # this is fast with available RAM (like 20 secs)
    if(verbose) { cat(" running reorder in system memory\n") }
    bigMat1 <- bigMat[to.order.r,to.order.c]
    if(verbose) {
      cat(" adding colnames\n") ; colnames(bigMat1) <- colnames(bigMat)[to.order.c] 
      cat(" adding rownames\n") ; rownames(bigMat1) <- rownames(bigMat)[to.order.r] 
      cat(" converting matrix to big.matrix\n") 
    }
    bigMat2 <- as.big.matrix(bigMat1, backingfile=bck.fn.o,
                             backingpath=dir$big, descriptorfile=des.fn.o)
    if(verbose) { cat(paste(" matrix descr saved as standard description file:",des.fn.o,"\n")) }
    descr <- describe(bigMat2)
  } else {
    #this is slow but creates backing file and will speed up ops later
    cat(" starting deep copy...")
    bigMat2 <- deepcopy(bigMat, cols = to.order.c, rows = to.order.r,
                        backingfile=bck.fn.o,backingpath=dir$big, descriptorfile=des.fn.o )
    cat("done\n")
    if(verbose) { cat("\nAdding names\n") }
    options(bigmemory.allow.dimnames=TRUE)
    colnames(bigMat2) <- colnames(bigMat)[to.order.c]
    if(verbose) { cat(" added colnames\n") }
    rownames(bigMat2) <- rownames(bigMat)[to.order.r]  
    if(verbose) { cat(" added rownames\n") }
    descr <- describe(bigMat2)
    flush(bigMat2) # hopefully this will ensure the row/colnames are added to the file backing
    if(verbose) { cat(paste(" due to use of deep copy option, recommend only to use descr saved as rbinary description file\n")) }
  }
  save(descr,file=R.descr)
  if(verbose) {
    cat(paste(" created big.matrix description file:",des.fn.o,"\n"))
    cat(paste(" created big.matrix backing file:",bck.fn.o,"\n"))
    cat(paste(" created big.matrix binary description file:",basename(R.descr),"\n"))
  } 
  while("trans.list" %in% search()) { detach(trans.list) }
  #return(descr) 
  return(R.descr)
}


thresh.check.valid <- function(nSD,lo.stat.nm=NULL,hi.stat.nm=NULL,mean.thr=NULL,dlrs.thr=NULL,gc.thr=NULL,check.only=F)
{
  #print("nSD");print(nSD)
  ### CHECK VALIDITY OF CHOSEN Mean, DLRS and GC THRESHOLDS ###
  # return validated set, or just T/F if check.only=T
  if(!is.null(lo.stat.nm) & !is.null(hi.stat.nm)) {
    lohi <- T
  } else {
    if(!is.null(mean.thr) & !is.null(dlrs.thr) & !is.null(gc.thr)) {
      lohi <- F
    } else {
      warning("check failed as not enough valid input parameters for test")
      return(F)
    }
  }
  valid <- T # set T which will be changed with any failure
  if(!lohi) {
    lo.stat.nm <- c(mean.thr[1],dlrs.thr[1],gc.thr[1])
    hi.stat.nm <- c(mean.thr[2],dlrs.thr[2],gc.thr[2])
  }
  ss <- !lo.stat.nm %in% boundary.stats(lower=T,nSD=nSD); 
  if(length(which(ss))>0) {
    if(any(!is.na(lo.stat.nm[ss]))) { warning("lower bound stat ",
        paste(narm(lo.stat.nm[ss]),collapse=",")," wasn't in boundary.stats() list: ",
        paste(boundary.stats(lower=T,nSD=nSD),collapse=",")) ; valid <- F }
    lo.stat.nm[ss] <- NA
  }
  ss <- !hi.stat.nm %in% boundary.stats(upper=T,nSD=nSD); 
  if(length(which(ss))>0) {
    if(any(!is.na(hi.stat.nm[ss]))) { warning("upper bound stat ",
        paste(narm(hi.stat.nm[ss]),collapse=",")," wasn't in boundary.stats() list: ",
        paste(boundary.stats(upper=T,nSD=nSD),collapse=",")) ; valid <- F }
    hi.stat.nm[ss] <- NA
  }
  if(!check.only) {
    return(list(lo.stat.nm,hi.stat.nm,valid))
  } else {
    return(valid)
  }
}


get.excl.list.from.stats <- function(stat.dat,stat.frame, dir="",
 lo.stat.nm=c("LB",NA,"LB"),hi.stat.nm=c("UB","UB","UB"),writeFiles=T,append=F,nSD=3)
{
 # calculate which samples in stat.dat should be excluded based on statistical thresholds
 # on the LRR data in stat.frame
 ### CALC EXCL SAMPLES ###
 headers <- colnames(stat.dat)
 excl.list <- list()
 sample.List <- row.names(stat.frame)
 dir <- validate.dir.for(dir,c("ano","excl"),warn=F)
 ch.list <- thresh.check.valid(nSD=nSD,lo.stat.nm=lo.stat.nm,hi.stat.nm=hi.stat.nm,check.only=F)
 lo.stat.nm <- ch.list[[1]]; hi.stat.nm <- ch.list[[2]]
 ## loop through mean, dlrs, gcwave to create exclusion lists ##
 for (dd in 1:length(headers))
 {
   excl.list[[dd]] <- character()

   if (!is.na(lo.stat.nm[dd])) {
     lo.stat <- stat.dat[lo.stat.nm[dd],headers[dd]] 
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[stat.frame[,headers[dd]]<lo.stat]) }

   if (!is.na(hi.stat.nm[dd])) {
     hi.stat <- stat.dat[hi.stat.nm[dd],headers[dd]] 
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[stat.frame[,headers[dd]]>hi.stat]) }

   if(writeFiles) {
     fnfn <- cat.path(dir$excl,headers[dd],ext="txt")
     if(!append) {
       writeLines(excl.list[[dd]],con=fnfn) 
       cat(paste("~wrote sample list:",fnfn,"\n"))
     } else {
       if(file.exists(fnfn)) {
         exist.lns <- readLines(fnfn) ####
       } else {
         exist.lns <- NULL
       }
       to.wrt <- unique(c(exist.lns,excl.list[[dd]]))
       num.added <- length(to.wrt)-length(exist.lns)
       writeLines(to.wrt,con=fnfn) 
       cat(paste("~appended",num.added,"to sample list:",fnfn,"\n"))     
     }
   }
 }
 fail.table <- matrix(logical(),nrow=nrow(stat.frame),ncol=length(lo.stat.nm))
 for (cc in 1:ncol(fail.table)) {
   fail.table[,cc] <- sample.List %in% excl.list[[cc]] }

 rownames(fail.table) <- sample.List; colnames(fail.table) <- headers
 return(fail.table)
}


get.fail.type.list <- function(rez,chk,labels=c("Mean","DLRS","GCWave"))
{
 ## GET LISTS OF SUBJECTS FAILING ON ALL POSSIBLE SETS OF CRITERIA
 ## USE FOR INDIVIDUAL GENOME vs LRR Plot Examples
 venn.lists <- list(); rnr <- rownames(rez)
 venn.lists[[1]] <- rnr[rez[,1] & !rez[,2] & !rez[,3] & chk<mean(chk,na.rm=T)]
 venn.lists[[2]] <- rnr[rez[,1] & !rez[,2] & !rez[,3] & chk>mean(chk,na.rm=T)]
 venn.lists[[3]] <- rnr[!rez[,1] & rez[,2] & !rez[,3]]
 venn.lists[[4]] <- rnr[!rez[,1] & !rez[,2] & rez[,3]]
 venn.lists[[5]] <- rnr[rez[,1] & rez[,2] & !rez[,3]]
 venn.lists[[6]] <- rnr[rez[,1] & !rez[,2] & rez[,3]]
 venn.lists[[7]] <- rnr[!rez[,1] & rez[,2] & rez[,3]]
 venn.lists[[8]] <- rnr[rez[,1] & rez[,2] & rez[,3]] 
 venn.lists[[9]] <- rnr[!rez[,1] & !rez[,2] & !rez[,3]]

 names(venn.lists)[1:2] <- paste(labels[1],c("Lo","Hi"),sep="")
 names(venn.lists)[3:4] <- labels[2:3]
 names(venn.lists)[5:7] <- c(paste(labels[1],labels[2],sep="+"),
                             paste(labels[1],labels[3],sep="+"),
                             paste(labels[2],labels[3],sep="+"))
 names(venn.lists)[8:9] <- c(paste(labels,collapse="+",sep=""),"None")
 return(venn.lists)
}


remove.boundary.snps.window <- function(snp.info,window=0,ucsc="hg18",chrLens=NULL)
{
 ## remove snps which when adding a window exceed the boundaries of their chromosome
 ## e.g, for analyses with windows calculations, e.g 1MB either side of each snp, etc
 uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
 if(is.null(chrLens)) { chrLens <- get.chr.lens(ucsc=ucsc) }
 if(is(snp.info)[1]=="RangedData"){
   chr.set <- chr.nums.of.info(snp.info); chr.n <- length(chr.set)
   to.remove.names <- NULL
   for (cc in 1:chr.n) {  
     too.big <- which(start(snp.info[cc])>(chrLens[cc]-window))
     # [1] 16210 16211 16212 16213 16214 16215 16216 16217 16218 16219 16220
     too.small <- which(start(snp.info[cc])<(window+1))
     to.remove <- c(too.big,too.small)
     if(length(to.remove)>0) {
       to.remove.names <- c(to.remove.names,rownames(snp.info[cc])[to.remove])
     }
   }
   if(length(to.remove.names)>0) {
     to.remove <- match(to.remove.names,rownames(snp.info))
     snp.info <- snp.info[-to.remove,]
   }
 } else {
     warning("not a RangedData object, function failed to remove boundary snps - GC code may crash")
 }
 return(snp.info)
}


sync.snpmat.with.info <- function(snpMatLst,snp.info=NULL,sample.info=NULL,dir=NULL,n.cores=1)
{
 # auto detect whether snpMatLst is a list of SnpMatrix objects or a list of RData
 # file locations and act accordingly. autodetect whether snp and/or sample info inputted.
 #print(length(snpMatLst)); print(is(snpMatLst)); print(dim(snp.info)); print(dim(sample.info)); print(head(dir))
 reorder.fn <- function(snpMatPart,info,snp=T) {
  if(snp) { nms <- colnames } else { nms <- rownames }
  to.keep <- match(rownames(info),nms(snpMatPart))
  to.keep <- to.keep[!is.na(to.keep)]
  if(snp) { return(snpMatPart[,to.keep]) } else { return(snpMatPart[to.keep,]) }
 }
 if(n.cores>1) { 
   require(multicore); l_apply <- function(...) { multicore::mclapply(...,mc.cores=n.cores) } 
 } else { 
   l_apply <- function(...) { lapply(...) }
 }
 if(!is.null(snp.info)) {
   if (is(snp.info)[1]=="RangedData" & is(snpMatLst[[1]])[1]=="SnpMatrix")
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
   if (is(sample.info)[1]=="data.frame" & is(snpMatLst[[1]])[1]=="SnpMatrix")
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


basic.cnv.counts <- function (cnv.hit.list, tot.cnvs=NA, content.txt="reference region", warn=F) {
  # do some basic counting for overlaps between a CNV list and a reference
  # mainly used as a subfunction of count.cnvs()
  n.cnvs <- length(cnv.hit.list)
  hits.by.cnv <- sapply(cnv.hit.list,clean.fn,fail=0,fn=length)
  total.gene.hits <- sum(hits.by.cnv)
  all.hits <- clean.fn(unlist(sapply(cnv.hit.list,unique)),NULL,c) #unique is for cnvs.by.gene
  all.unique.genes <- unique(all.hits)
  total.unique.hits <- length(all.unique.genes)
  cnvs.by.gene <- table(all.hits,dnn=NULL)
  num.cnvs.at.least.one.hit <- length(which(hits.by.cnv!=0))  
  if(is.na(tot.cnvs)) { 
    if(length(which(hits.by.cnv==0))==0) { 
      if(warn) { warning("'tot.cnvs' not provided and 'cnv.ranges' may only contain cnvs with hits, or cnvs passing min.snp criteria")  }
    }
    tot.cnvs <- n.cnvs 
  }
  ## format output results ##
  result1 <- list(hits.by.cnv,all.unique.genes,cnvs.by.gene,
                  total.gene.hits,total.unique.hits,num.cnvs.at.least.one.hit,tot.cnvs)
  nc <- nchar(content.txt)
  if(substring(tolower(content.txt),nc-3) %in% c("xons","ions","enes")) { content.txt <- substr(content.txt,1,nc-1) }  # remove plural s
  descr.names <- c(paste(content.txt," count by CNV",sep=""),
                   paste("list of unique ",content.txt,"s across all CNVs",sep=""),
                   paste("CNV count by ",content.txt,sep=""),
                   paste("number of ",content.txt," hits across all CNVs",sep=""),
                   paste("number of unique ",content.txt,"s across all CNVs",sep=""),
                   paste("number of CNVs with at least one ",content.txt," hit",sep=""),
                   "total number of CNVs")
  names(result1) <- descr.names
  return(result1)
}


sample.cnv.counts <- function(cnv.hit.list, sids, content.txt) {
  # do some basic counting for overlaps between a CNV list and a reference with respect to samples
  # mainly used as a subfunction of count.cnvs()
  s.names <- c("ids with at least 1 CNV", paste("ids with at least 1 ",content.txt," hit",sep=""),
               "number of ids with at least 1 CNV",
               paste("number of ids with at least 1 ",content.txt," hit",sep=""))
  hits.per.cnv <- sapply(cnv.hit.list,clean.fn,fail=0,fn=length)
  all.ids <- sids
  all.unique.ids <- unique(all.ids)
  num.any.cnv <- length(all.unique.ids)
  ids.at.least.one.hit <- unique(sids[hits.per.cnv>0])
  num.ids.with.hits <- length(ids.at.least.one.hit)
  sample.result <- list(all.unique.ids,ids.at.least.one.hit,num.any.cnv,num.ids.with.hits)
  names(sample.result) <- s.names
  return(sample.result)
}


reduce.list.to.scalars.old <- function(ll) {
  if(is.list(ll)) {
    for (cc in length(ll):1) {
      if(is(ll[[cc]])[1]=="list") { 
        ll[[cc]] <- reduce.list.to.scalars.old(ll[[cc]]) 
      } else {
        if(length(ll[[cc]])>1 | !is.numeric(ll[[cc]])) { ll[[cc]] <- NULL }
      }
    }
  }
  return(ll)
}

nmoo <- function(x) { x <- rep(NULL,1) }

reduce.list.to.scalars.2 <- function(ll,max.ln=1000) {
  if(is.list(ll)) {
    if(length(ll)>max.ln) { return(sapply(ll,nmoo)) }
    for (cc in length(ll):1) {
      if(is(ll[[cc]])[1]=="list") { 
        ll[[cc]] <- reduce.list.to.scalars.2(ll[[cc]]) 
      } else {
        if(length(ll[[cc]])>1 | !is.numeric(ll[[cc]])) { ll[[cc]] <- NULL }
      }
    }
  }
  return(ll)
}

rlts <- function(X) {
  if(is(X)[1]=="list") { 
    X <- reduce.list.to.scalars(X) 
  } else {
    if(length(X)>1 | !is.numeric(X)) { X <- NULL }
  }
  return(X)
}

reduce.list.to.scalars <- function(ll) {
  if(is.list(ll)) {
    ll <- rev(lapply(rev(ll),rlts))
  }
  return(ll)
}



get.overlap.stats.chr <- function (next.chr, x, y, xx, yy, name.by.gene=T, rel.query=T, txid=F, prog=F) {
  ## mainly to tidy up function below - does the processing for 1 chromosome
  olp <- findOverlaps(x[[next.chr]],y[[next.chr]])  
  if(length(olp)==0 | length(subjectHits(olp))==0 | length(queryHits(olp))==0) {
    return(list(NULL,NULL,NULL)) } # if no matches return list of 3 NULL elements
  overlap.ranges <- ranges(olp,x[[next.chr]],y[[next.chr]])
  genes.per.cnv.n <- tapply(queryHits(olp),subjectHits(olp),c)
  match.num.to.names <- function(num,ind) { return(ind[num]) }
  if(name.by.gene) {
    if(txid & ("txname" %in% colnames(xx))) {
      # transcript id when using exons
      ind <- xx[next.chr]$txname
      genes.per.cnv <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
    } else {
      ind <- xx[next.chr]$gene
      genes.per.cnv <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
    }
  } else { 
    if(T | length(grep("DGV",toupper(rownames(xx)[1:2])))>0) {
      ind <- rownames(xx[next.chr])
      genes.per.cnv <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
    } else {
      # use row-numbers for matchs unless using DGV ids
      genes.per.cnv <- genes.per.cnv.n 
    }
  }
  gene.overlaps.per.cnv <- tapply(width(overlap.ranges),subjectHits(olp),c)
#  if(T|paste(next.chr)=="21") { print(gene.overlaps.per.cnv) }
  select.overlappers <- as.numeric(names(genes.per.cnv)) # there are more of these than rownames y and they are unique!
#  print(paste("selectoverlappers",length(rownames(yy[next.chr])),length(select.overlappers),length(unique(select.overlappers))))
  if(rel.query) {
    # width relative to query (e.g, Gene %)
 #  gene.lengths <- width(x[[next.chr]])
    ind <- width(xx[next.chr])
    gene.widths <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
    if(length(gene.widths)>=1 & length(gene.overlaps.per.cnv)>=1) {
      pc.per.cnv <- sapply(1:length(gene.widths), function(x) { gene.overlaps.per.cnv[[x]]/gene.widths[[x]] } ,simplify=F)
    } else {
      pc.per.cnv <- gene.widths # i.e, make an equivalent null/na/0 value
    }
  } else {
    # width relative to subject (e.g, CNV %)
    cnv.lengths <- width(y[[next.chr]][select.overlappers])
    if(length(cnv.lengths)>=1 & length(gene.overlaps.per.cnv)>=1) {
      pc.per.cnv <- sapply(1:length(cnv.lengths), function(x) { gene.overlaps.per.cnv[[x]]/cnv.lengths[x] } ,simplify=F)
    } else {
      pc.per.cnv <- cnv.lengths # i.e, make an equivalent null/na/0 value
    }
  }    
 # print(tail(rownames(yy[next.chr])[select.overlappers]))
  names(pc.per.cnv) <- names(genes.per.cnv) <- 
    names(gene.overlaps.per.cnv) <- rownames(yy[next.chr])[select.overlappers] # the dodgy line
  by.chr.list <- list(genes.per.cnv,gene.overlaps.per.cnv,pc.per.cnv)
 # print(head(names(by.chr.list[[1]])))
  if(prog) { cat(".") }
  return(by.chr.list)
}



fill.blank.matches <- function(newlist,missing.val="",full.len=NA)
{
  #ensure a blank entry for cnvs with no overlaps
  # full.len is an optional desired out-length for the list (which if trailing cells are blank might differ from the guess)
##head(paste(names(newlist)))
  to.nums <- as.numeric(paste(names(newlist)))
  if(length(narm(to.nums))<length(newlist)) { warning("some list names were not numbers") }
  if(length(narm(to.nums))<1) { return(newlist) }
  mm <- max(c(to.nums,full.len),na.rm=T);  if(is.infinite(mm)) {  warning("couldn't find end of list (empty?)"); return(newlist) }
  no.gene.list <- which(!paste(c(1:mm)) %in% names(newlist))
  if(length(no.gene.list)>0) {
    none.list <- vector("list",length(no.gene.list))
    names(none.list) <- paste(no.gene.list)
    none.list <- sapply(none.list,function(x) { missing.val } )
    out <- c(newlist,none.list)
    out <- out[order(as.numeric(names(out)))]
  } else {
    out <- newlist
  }
  return(out)
}


add.color.to.snp.info <- function(snp.info,scheme=c("mono","alt","unique","hilite"),
                                  col1="purple",col2="orange",hilite=6) {
  ## add colour scheme for plotting cnv/snp data ##
  # mono = all col1, alternate = alternate between col1 and col2, unique = all different, hilite all col2 except 'hilite' chr(s)
  must.use.package("genoset",T)
  if(!col1 %in% colors()) { col1 <- "purple" }
  if(!col2 %in% colors()) { col2 <- "orange" }
  if(is(snp.info)[1]!="RangedData") { stop("snp.info must be a RangedData object to add colours") }
  snp.info <- toGenomeOrder(snp.info,strict=T)
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
  if(n.chr<1) { stop("there seem to be no chromosomes in the snp.info object") }
  scheme <- scheme[scheme %in% c("mono","alt","unique","hilite")][1]
  scheme[is.na(scheme)] <- "mono"
  if(scheme=="hilite") {
    if(!all(paste(hilite) %in% paste(chr.set))) { hilite <- chr.set[1]; warning("hilite value not valid, defaulting to first chr in set") }
  }
  if(!is(snp.info)[1]=="RangedData") {
    warning("object wasn't a snp.info object, color not added")
    return(snp.info)
  }
  ## define different colour scheme options ##
  monoCol <- function(x,c1,c2) { rep(c1,times=length(x)) }
  altCol <- function(x,c1,c2) { rep(c(c1,c2),times=ceiling(length(x)/2))[1:length(x)] }
  uniqCol <- function(x,c1,c2) { 
    coloz22 <- get.distinct.cols(22) 
    out <- rep("black",times=length(x))
    which.chr <- which(paste(1:22) %in% x)
    if(length(which.chr)>0) { out[1:length(which.chr)] <- coloz22[which.chr] }
    return(out)
  }
  hiCol <- function(x,c1,c2) { out <- monoCol(c2); out[hilite] <- c1; return(out) }
  col.func <- switch(scheme, mono=monoCol, alt=altCol, unique=uniqCol, hilite=hiCol)
  coloz <- col.func(chr.set,col1,col2)
  repz <- table(snp.info$space)
  repz <- as.numeric(repz); repz <- repz[!is.na(repz) & repz>0]
  snp.info[["color"]] <- rep(coloz,times=as.numeric(repz))
  return(snp.info)
}


add.gindx.to.Ranges <- function(snp.info,ucsc="hg18",absolute=T,label="gindx") {
  must.use.package("genoset",T)
  if(!is(snp.info)[1]=="RangedData") {
    warning("object wasn't a RangedData object, genome index not added")
    return(snp.info)
  }
  snp.info <- toGenomeOrder(snp.info,strict=T)
  if(!absolute) {
    snp.info[[paste(label[1])]] <- genoPos(snp.info)
    return(snp.info)
  }
  uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  num.chr <- length(snp.info)
  which.chr <- which(paste(1:22) %in% chrNames(snp.info))
  all.chr.len <- get.chr.lens(dir,ucsc=ucsc)
  if(length(which.chr)>0) {
    chrLens <- all.chr.len
    chr.indx <- chrIndices(snp.info[paste(which.chr)])
    chr.counts <- (chr.indx[,2]-chr.indx[,1]+1)
    chrStarts <- c(0,cumsum((chrLens)))[which.chr]
    chrStarts <- chrStarts[1:length(which.chr)] # trim last position which is an end, not a start
  } else { 
    chrLens <- chrStarts <- NULL
  }
  which.notC <- chrOrder(chrNames(snp.info)[(!chrNames(snp.info) %in% paste(1:22))])
  if(length(which.notC)>0) {
    chr.pos.not <- chrInfo(snp.info[which.notC])
    chr.indx2 <- chrIndices(snp.info[which.notC])
    chr.counts2 <- (chr.indx2[,2]-chr.indx2[,1]+1)
    chrStarts2 <- sum(all.chr.len) + chr.pos.not[,3] # plus offset
  } else {
    chr.counts2 <- chrStarts2 <- NULL
  }
  chr.counts <- c(chr.counts,chr.counts2)
  chrStarts <- c(chrStarts,chrStarts2) # genome position of chromosome starts
  gindx <- start(snp.info[c(paste(which.chr),which.notC)])+rep(chrStarts,times=chr.counts)
  snp.info[[paste(label[1])]] <- gindx
  return(snp.info)
}


filter.snp.info <- function(snp.info,filt.names=NULL,include=T,dir=NULL,verbose=T) {
  # select(include=T) or remove(F) set of 'filt.names' from snp.info
  must.use.package("genoset",T)
  if(!is(snp.info)[1]=="RangedData") {
    warning("object wasn't a RangedData object, genome index not added")
    return(snp.info)
  }
  if(is.null(filt.names) | !is.character(paste(filt.names))) {
    filt.names <- rownames(snp.info)
    if(verbose) { cat("filt.names invalid or blank so selecting all snps and ordering\n") }
  } 
  filt.names <- paste(force.vec(filt.names,dir=dir))
  filt.names.valid <- filt.names[(filt.names %in% rownames(snp.info))]
  if(length(filt.names.valid)>0) {
    pc.ok <- round((length(filt.names.valid)/length(filt.names))*100,1)
    snp.info <- toGenomeOrder(snp.info,strict=T)
    if(include) {
      if(verbose) { cat(pc.ok,"% of 'filt.names' found in snp.info object, selecting these to create new snp.info object\n",sep="") }
      snp.info <- snp.info[rownames(snp.info) %in% filt.names.valid,]
    } else {
      if(verbose) { cat(pc.ok,"% of 'filt.names' found in snp.info object, removing these to create new snp.info object\n",sep="") }
      snp.info <- snp.info[!rownames(snp.info) %in% filt.names.valid,]
    }
  } else {
    warning(paste("no valid names found in filt.names [",paste(filt.names[1:3],collapse=","),",..] for snp.info:\n",sep=""))
    print(head(snp.info))
    cat("returning object unchanged\n")
  }
  return(snp.info)
}


get.chr.info.filt <- function(dir, extended=T, filt.names=NULL,verbose=F,
                             snp.info.fn="snpdata.map",ucsc="hg18",autosomes.only=F,
                             add.gindx=T,absolute=T,snp.fn="snpNames.txt",anot="map3",
                             scheme=c("unique","alt","mono"),col1="black",col2="grey",...)
{
 # calculate chromosome indexes, in particular for the filtered version of snplist rownames(bigMat2)
 # basically does nothing if recalcF == FALSE
 must.use.package(c("genoset","IRanges"),bioC=T)
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 add.dir.if.not <- reader:::add.dir.if.not # get internal function from reader
 snp.info.fn <- add.dir.if.not(snp.info.fn,dir$ano,dir) # would also take entire 'dir'
 #  [ID and position for each SNP within chromosomes]
 snp.info <- import.marker.data(dir,markerinfo.fn=snp.info.fn,snp.fn=snp.fn,
                                                  anot=anot,verbose=verbose,...)
 universe(snp.info) <- paste(ucsc) # annotate UCSC reference
 # potentially filter by a list, add an absolute/relative index, add color information for plotting
 if(!is.null(filt.names)) { snp.info <- filter.snp.info(snp.info,filt.names=filt.names,dir=dir) } 
 if(autosomes.only) { snp.info <- select.autosomes(snp.info) }
 if(extended | add.gindx) { snp.info <- add.gindx.to.Ranges(snp.info,ucsc=ucsc,absolute=absolute) }
 if(all(!is.na(scheme))) { snp.info <- add.color.to.snp.info(snp.info,scheme=scheme,col1=col1,col2=col2) }
 if(verbose) {
   cat(" snp.info object produced with",nrow(snp.info),"snps and",length(snp.info),"chromosomes\n")
 }
 return(snp.info)
}


sample.bad.count.table <- function(dir,sub.list,refresh=T,type=1,extralist=NA,addExtraTo="Mean")
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason for each plate
 dir <- validate.dir.for(dir,c("ano","excl"),warn=F)
 nsamp <- length(sub.list)
 if(type==1) {
   new.fr <- data.frame(CallRate=integer(nsamp),
                        Mean=integer(nsamp),DLRS=integer(nsamp),ChrAb=integer(nsamp),
                        GCWave=integer(nsamp),TOTAL=integer(nsamp))  #StDev=integer(nsamp),
 } else {
   new.fr <- data.frame(Mean=integer(nsamp),
                        DLRS=integer(nsamp),ChrAb=integer(nsamp),
                        GCWave=integer(nsamp),TOTAL=integer(nsamp)) #StDev=integer(nsamp),
 }
 rownames(new.fr) <- sub.list

 if (!is.na(extralist[1]) & addExtraTo %in% colnames(new.fr))
 { extraCol <- match(addExtraTo,colnames(new.fr)) } else { extraCol <- 0 }
 dirfilz <- list.files(dir$excl)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- agrep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   if (length(nxt.col)>0) {
     next.bad[[cc]] <- readLines(cat.path(dir$excl,dirfilz[cc]))
     if (extraCol==nxt.col) { next.bad[[cc]] <- c(next.bad[[cc]],extralist) }
     bad.cnts <- (table(next.bad[[cc]],exclude=NA,useNA="no"))
     new.fr[names(bad.cnts),nxt.col] <- as.vector(bad.cnts)
   }
 }
 TOTAL.bad <- unique(unlist(next.bad))
 bad.cnts <- (table(TOTAL.bad,exclude=NA,useNA="no"))
 new.fr[names(bad.cnts),ncol(new.fr)] <- as.vector(bad.cnts)

 return(new.fr)
}


list.qc.fail.types <- function(dir,filt.list)
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason
 dir <- validate.dir.for(dir,c("ano","excl"),warn=F)
 ns <- 1
 new.fr <- data.frame(CallRate=integer(ns),
                        Mean=integer(ns),DLRS=integer(ns),ChrAb=integer(ns),
                        GCWave=integer(ns),TOTAL=integer(ns))  #StDev=integer(ns),
 dirfilz <- list.files(dir$excl)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- agrep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   next.bad[[cc]] <- readLines(cat.path(dir$excl,dirfilz[cc]))
   if (length(which(!next.bad[[cc]] %in% filt.list))>0) {
     next.bad[[cc]] <- next.bad[[cc]][-which(!next.bad[[cc]] %in% filt.list)]
   } 
   bad.cnts <- length(next.bad[[cc]])
   new.fr[1,nxt.col] <- bad.cnts
 }
 TOTAL.bad <- unique(unlist(next.bad))
 bad.cnts <- length(TOTAL.bad)
 new.fr[1,ncol(new.fr)] <- bad.cnts

 return(new.fr)
}


update.plate.bad.count.table <- function(dir,plate.index=NULL,plate.list="plate.list.txt",
                                         refresh=T,filt.list=NULL,type="plate",incl.cr=T)
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason for each plate
 # plate.index should contain sample ids in column 1, plate ids in column 2
 # would be easy to run for wells by switching cols 2&3
 dir <- validate.dir.for(dir,c("ano","excl"),warn=F)
 if(is.matrix(plate.index)) { plate.index <- as.data.frame(plate.index) }
 plate.col <- which(tolower(colnames(plate.index)) %in% type)[1]
 is.plate <- T
 if(is.na(plate.col)) { 
   plate.col <- 2  # default to col 2 if none named plate
 } else {
   if(tolower(type)!="plate") {
     cat("\nMaking failure count table for",type,"\n")
     is.plate <- F
   }
 }
 if (!is.null(filt.list)) {
   fl <- find.id.col(plate.index,filt.list,ret="index")
   plate.index.sub <- plate.index[fl$index,]
 } else {
   plate.index.sub <- plate.index
 }
 if (is.plate & !refresh & ("plate.stats.txt" %in% list.files(dir$ano))) {
   new.fr <- read.delim(file=cat.path(dir$ano,"plate.stats.txt"))
   plts <- paste(new.fr[,1]) # this should be the column with plate ids
 } else { 
   if(length(plate.list)==1) {
     if((plate.list %in% list.files(dir$ano)) & is.plate)
     {
       plat.fn <- cat.path(dir$ano,plate.list)
       plts  <- readLines(plat.fn)
     } else {
       plts <- paste(unique(plate.index[,plate.col])) # plate ids should be in col 2
     }
   } else {
     plts <- plate.list # assume this is the list of plates/other batch category
   }
   npl <- length(plts)
   new.fr <- data.frame(ID=plts,CallRate=integer(npl),
                        Mean=integer(npl),DLRS=integer(npl),ChrAb=integer(npl),
                         GCWave=integer(npl),TOTAL=integer(npl)) 
   if(!is.null(filt.list) & !is.null(plate.index))
   {
   	per.plate.cnt <- table(plate.index.sub[,plate.col])
   	new.fr[["SIZE"]] <- per.plate.cnt[paste(plts)]
   }
 }
 dirfilz <- list.files(dir$excl)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- grep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   if(length(nxt.col)==1) {
	   next.bad[[cc]] <- readLines(cat.path(dir$excl,dirfilz[cc]))
	   plt.bad <- (plate.index[(plate.index[,1] %in% next.bad[[cc]]),plate.col])
	   plt.bad2 <- (plate.index.sub[(plate.index.sub[,1] %in% next.bad[[cc]]),plate.col])
   	 bad.cnts <- (table(plt.bad))
     if(length(bad.cnts)!=0) {
   	   new.fr[match(names(bad.cnts),new.fr[,1]),nxt.col] <- as.vector(bad.cnts)
     }
   	 cat("",dirfilz[cc]," ",sum(table(plt.bad2)),.Platform$file.sep,
          length(next.bad[[cc]])," in this cohort\n",sep="")
	} else {
		if(length(nxt.col)>1) { 
			stop("Error: more than 1 file matched header. check files") 
		} else {
			warning("file ",dirfilz[cc]," skipped in counts because not in table header")
			next.bad[[cc]] <- ""
		}
	}
   if (!is.null(filt.list) & nchar(next.bad[[cc]][1])>1) {
     excl <- which(!next.bad[[cc]] %in% filt.list)
     if(length(excl)>1) {
       next.bad[[cc]] <- next.bad[[cc]][-excl]  }
   }
 }
 TOTAL.bad <- unique(unlist(next.bad))
 plt.bad <- (plate.index[(plate.index[,1] %in% TOTAL.bad),plate.col])
 plt.bad <- plt.bad[plt.bad!=""]
 bad.cnts <- (table(plt.bad))
 new.fr[match(names(bad.cnts),new.fr[,1]),"TOTAL"] <- as.vector(bad.cnts)
 if(incl.cr) {
   new.fr[["PC+CR"]] <- round((new.fr$TOTAL+new.fr$CallRate)/(new.fr$SIZE+new.fr$CallRate),2)
 } else {
   new.fr[["PC"]] <- round((new.fr$TOTAL)/(new.fr$SIZE),2)
 }
 return(new.fr)
}


do.lens.summary <- function(lenz,dat=NULL,print.summary=T)
{
  # prints a simple summary of an inputted vector of CNV lengths
  my.len <- vector("list",7)
  names(my.len) <- c("length >0 and <=1kb","length >1kb and <=10kb",
  "length >10kb and <=50kb", "length >50kb and <=100kb",
  "length >100kb and <=500kb", "length >500kb and <=5mb","length >5mb ")
  my.len[[1]] <- length(which(lenz>0 & lenz<=1000))
  my.len[[2]] <- length(which(lenz>1000 & lenz<=10000))
  my.len[[3]] <- length(which(lenz>10000 & lenz<=50000))
  my.len[[4]] <- length(which(lenz>50000 & lenz<=100000))
  my.len[[5]] <- length(which(lenz>100000 & lenz<=500000))
  my.len[[6]] <- length(which(lenz>500000 & lenz<=5000000))
  my.len[[7]] <- length(which(lenz>5000000))
  if(!is.null(dat)) {
    dat.want <- dat[,c("IID","CHR","BP1","BP2")]
    cat("\nIDs below 1kb","\n")
    print(cbind(dat.want,lenz)[(which(lenz>0 & lenz<=1000)),])
    cat("IDs >500kb and <=5mb","\n")
    print(cbind(dat.want,lenz)[(which(lenz>500000 & lenz<=5000000)),])
    cat("IDs above 5mb","\n")
    print(cbind(dat.want,lenz)[(which(lenz>5000000)),])
  }
  if(print.summary) { print(my.len) }
  return(my.len)
}


force.pheno.codes <- function(pheno,missing=c(0,-9,-99,99),verbose=F) {
  ## try to cleverly assign phenotype codes
  nph <- as.numeric(paste(pheno)); if(length(narm(nph))>.5*length(pheno)) { pheno <- nph }
  if(!is.numeric(pheno)) {
    ctrl.codes <- c("ctrl","cont","cntr","healthy","norm")
    warning("found phenotype coding as text, assuming any phenotype values are",
          " not cases if they contain following strings: ",paste(ctrl.codes,collapse=","))
    pheno2 <- paste(pheno); pheno <- rep(2,times=length(pheno2))
    for (cc in 1:length(ctrl.codes)) {
      is.ctrl <- grep(ctrl.codes[cc],tolower(pheno2))
      pheno[is.ctrl] <- 1
    }
    pheno[is.na(pheno2)] <- 0
  }
  cur.types <- sort(unique(pheno))
 
  if(length(cur.types)==3 & length(pheno[is.na(pheno)])==0) {
    if(all(cur.types==c(0,1,2))) {
      # coding already seems to be in proper plink format (including missing)
      warning("assuming phenotype coding: 1=controls, 2=cases, 0=missing")
      return(pheno)
    }
  }
  if(length(cur.types)==2 & length(pheno[is.na(pheno)])==0) {
    if(all(cur.types==c(1,2))) {
      # coding already seems to be in proper plink format (excluding missing)
      warning("assuming phenotype coding: 1=controls, 2=cases")
      return(pheno)
    }
  }
  if(any(pheno>2)) {
  pheno[pheno %in% missing] <- NA } # even though we just coded as zero, necessary for 'min' value to work
  lwr <- min(pheno,na.rm=T); 
  if(length(unique(lwr))==1 & lwr==2) { lwr <- 1 }  # if only cases, user may enter just 2s
  # guess at affected and unaffected categories
  warning(" assuming phenotype code ",lwr," are 'unaffected', and >",lwr," are 'affected' samples\n",
    "   if this is wrong, pls ctrl-C and edit 'sampleinfo.tab' to make the minimum value unaffected")
  pheno2 <- pheno; pheno2[pheno==lwr] <- 1
  pheno2[pheno>lwr] <- 2; pheno2[is.na(pheno)] <- 0
  return(pheno2)
}


force.sex.codes <- function(sex,missing=c(0,-9,-99,99),verbose=F) {
  ## try to cleverly assign sex codes
  if(!is.numeric(sex)) {
    if(verbose) { cat(" found sex coding as text, converting to 0,1,2 unknown/male/female") }
    sex2 <- paste(sex); sex <- rep(0,times=length(sex2))
    is.girlz <- grep("F",toupper(sex2))
    is.guyz <- grep("M",toupper(sex2))
    sex[is.guyz] <- 1
    sex[is.girlz] <- 2
  } 
  if(length(unique(sex))<3) {
    if(length(which(sex %in% missing))>0) {
      warning("found only sex codes of ",paste(unique(sex),collapse=","),
              "but missing code(s) are: ",paste(unique(missing),collapse=","),
              "will set those with missing codes to missing. fix sex coding if this is wrong")
    }
  }
  sex[sex %in% missing] <- NA # even though we just coded as zero, necessary for 'min' value to work
  lwr <- min(sex,na.rm=T)
  if(length(unique(lwr))==1 & lwr==2) { lwr <- 1 } # if only female, user may enter just 2s
  # guess at sex coding
  warning("assuming sex code ",lwr," are males, and others are females",
    "\n   if this is wrong, pls ctrl-C and edit 'sampleinfo.tab' to code males as the minimum value")
  sex2 <- rep(2,length(sex))
  sex2[sex==lwr] <- 1; sex2[is.na(sex)] <- 0
  return(sex2)
}


make.fam.file <- function(sample.info,dir,out.fn="plink.fam")
{
  ## create plink family file using sample.info
  #Family ID,Sample ID,Paternal ID,Maternal ID,
  #Sex (1=male; 2=female; other=unknown,"\n")
  #Affection (0=unknown; 1=unaffected; 2=affected,"\n")
  dir <- validate.dir.for(dir,"cnv.plk")
  PED.FAM.TEMP <- sample.info[,rep("grp",6)]
  PED.FAM.TEMP[,6] <- 1
  PED.FAM.TEMP[,c(3,4,5)] <- 0
  PED.FAM.TEMP[,2] <- rownames(sample.info)
  PED.FAM.TEMP[,1] <- 1 #c(1:nrow(PED.FAM.TEMP))
  search.str <- which(tolower(colnames(sample.info)) %in% c("pheno","phenotype","case"))
  if(length(search.str)>0) { 
    pheno <- sample.info[,search.str[1]]
    PED.FAM.TEMP[,6] <- force.pheno.codes(pheno)
  } else {
    warning("no phenotype found in 'sampleinfo.tab'. All samples set to 'unaffected'")
  }
  search.str <- which(tolower(colnames(sample.info)) %in% c("sex","gender"))
  if(length(search.str)>0) { 
    sex <- sample.info[,search.str[1]]
    PED.FAM.TEMP[,5] <- force.sex.codes(sex)
  } else {
    warning("no sex found in 'sampleinfo.tab'. All sexes set to unknown")
  }
  colnames(PED.FAM.TEMP) <- c("FamilyID","SampleID","PaternalID","MaternalID","Sex","Phenotype")
  rownames(PED.FAM.TEMP) <- NULL
  if(out.fn=="") {
    return(PED.FAM.TEMP) 
  } else {
    ofn <- cat.path(dir$cnv.qc,out.fn,ext="fam")
    write.table(PED.FAM.TEMP,file=ofn,quote=F,col.names=F,row.names=F)
    cat(paste("~wrote fam file:",ofn,"\n"))
    return(ofn) 
  }
}

#### to check!!! HERE ###
filter.hifreq.plates <- function(plink.files, sample.info, rem.file=plink.files[1], dir, psd=3,
                                 QC.ON=T, rmv.low.plates=F, suf.for.rmv="pltCNV", 
                                  plot.fn="plateqchist.pdf",rare=F)
{
  # filter samples with too many rare DELs/DUPs (according to specified rare rates) (poisson)
  # or filter samples with too many CNVs (normal)
  dir <- validate.dir.for(dir,c("cnv.plk"),warn=F)
  if(!"QCfail" %in% colnames(sample.info)) { warning("sample.info must have column 'QCfail") }
  nsamps <- length(which(sample.info$QCfail==0))
  rem.file <- paste(rem.file[1])
  pl.path <- cat.path(dir$cnv.qc,plink.files,ext="cnv")
  if(!all(file.exists(pl.path))) { warning("plink file(s) not found") ; return(NULL) }
  if(!is.file(rem.file,dir$cnv.qc)) { warning("'rem.file' (",rem.file,") file not found") ; return(NULL) }
  dat1 <- list()
  for (cc in 1:length(plink.files)) {
    dat1[[cc]] <- read.table(pl.path[cc],header=T)  #; colnames(dat) <- 
  }
  dat <- do.call(rbind,args=dat1)
  plate.list <- get.plate.info(dir,verbose=F); 
  if(is.null(plate.list)) { stop("couldn't do plate QC as 'plate' column not found in sample.info") }
  plate.lookup <- plate.list[[1]];  n.per.plate <- plate.list[[3]]
  dat[["PLATE"]] <- plate.lookup$plate[match(dat$IID,plate.lookup$id)]
  tt <- as.numeric(table(dat$PLATE))[order(as.numeric(table(dat$PLATE)))]
  ss <- names(table(dat$PLATE))[order(as.numeric(table(dat$PLATE)))]
  nn <- (n.per.plate[match(ss,names(n.per.plate))]); #nn <- rep(96,length(ss))
  uu <- tt/nn
  mps <- max(n.per.plate) # get bounds as half and quarter of full plate size
  bounds <- c(round(mps/4),round(mps/2)); sel <- vector("list",3); sdz <- numeric(3)
  sel[[1]] <- nn<bounds[1] ; sel[[2]] <- (nn>=bounds[1] & nn<bounds[2]) ; sel[[3]] <- nn>=bounds[2]
  # if less than 5 samples in subranges of plate size then concatenate into the next range
  if(length(which(sel[[1]]))<5) { sel[[2]] <- sel[[1]] | sel[[2]]; sel[[1]] <- rep(F,length(nn)) }
  if(length(which(sel[[2]]))<5) { sel[[3]] <- sel[[2]] | sel[[3]]; sel[[2]] <- rep(F,length(nn)) }
  for (dd in 1:length(sdz)) { sdz[dd] <- sd(uu[sel[[dd]]],na.rm=T) } # NA when none in range
  mn <- mean(uu,na.rm=T) 
  sdv <- rep(sdz,(sapply(sel,function(x) { length(which(x)) }))) # ordered SDs from ranges 1 to 3 for plotting - NA's will disappear
  ubz <- mn+(psd*unique(sdv)); lbz <- mn-(psd*unique(sdv))
  # depending on whether plate numbers vary a lot, either set one threshold, or multiple
  if(length(uu)>0) {
    pdf(cat.path(dir$cnv.qc,plot.fn))
    if(length(narm(ubz))<2) {
      hist((uu),xlab="samples per plate",main="Number of CNVs per Plate",cex=1.5)
      abline(v=c(narm(lbz),narm(ubz)),lty="dashed")
      pass.pos <- uu<=narm(ubz); pass.neg <- uu>=narm(lbz)
    } else {
      plot(nn[order(nn)],uu[order(nn)], ylim=range(mn-psd*sdz,mn+psd*sdz,na.rm=T)*c(1,1.1),
           ylab="CNVs per sample per plate", xlab="samples per plate",
           ,main="Plates versus number of CNVs",cex=1.5)
      pass.pos <- pass.neg <- rep(T,length(ss))
      for (dd in 1:length(ubz)) { pass.pos[sel[[dd]]] <- uu[sel[[dd]]]<=ubz[dd] }
      for (dd in 1:length(lbz)) { pass.neg[sel[[dd]]] <- uu[sel[[dd]]]>=lbz[dd] }
      # points(nn[!pass.pos],uu[!pass.pos],col="red",pch=".",cex=3)
      # points(nn[!pass.neg],uu[!pass.neg],col="blue",pch=".",cex=3)
      lines(nn[order(nn)],mn+(psd*sdv),lty="dashed")
      lines(nn[order(nn)],mn-(psd*sdv),lty="dotted")
      legend("top",legend=paste(c("over","under"),"sensitivity threshold (",c("+","-"),psd,"SD)",sep=""),
             lty=c("dashed","dotted"),bty="n")
    } 
    dev.off()
  } else { pass.pos <- uu<=narm(ubz); pass.neg <- uu>=narm(lbz) }
  toohi <- cbind(ss[!pass.pos], tt[!pass.pos], nn[!pass.pos], round(uu[!pass.pos],2))
  toolo <- cbind(ss[!pass.neg], tt[!pass.neg], nn[!pass.neg], round(uu[!pass.neg],2))
  colnames(toohi) <- colnames(toolo) <- c("Plate","CNVs","Samples","CNVs.per.Sample")
  toohi <- as.data.frame(toohi); toolo <- as.data.frame(toolo)
  # mean/sd/'3'SD for dels/dups per subj
  cat("\nFiltering samples from plates with a large number of ",{if(rare) "rare " else ""},"CNVs per sample\n",sep="")
  if(rare) { cat("[NB: tagged plates may be same as filtered above for all CNVs]\n") }
  ###dat <- read.table(file.choose(),header=T)
  cat("\nPlates with too many CNVs (oversensitive)\n")
  if(nrow(toohi)>0) {
    print(as.data.frame(toohi),quote=F,row.names=F)
  } else {
    cat("[none exceeded boundary]\n")
  }
  
  cat("\nPlates with too few CNVs (undersensitive)\n")
  if(nrow(toolo)>0) {
    print(as.data.frame(toolo),quote=F,row.names=F)
  } else {
    cat("[none exceeded boundary]\n")
  }
  
  if(QC.ON) {
    to.drop <- toomanyers <- plate.lookup$id[plate.lookup$plate %in% toohi$Plate]
    cat(paste("\n list of",length(toomanyers),"samples from plates with too many CNVs from",paste(plink.files,collapse=",")),":\n")
    print(toomanyers,quote=F); cat("\n")
    if(rmv.low.plates) {
      toofewers <- plate.lookup$id[plate.lookup$plate %in% toolo$Plate]
      cat(paste("\n list of",length(toofewers),"samples from plates with too few CNVs from",paste(plink.files,collapse=",")),":\n")
      print(toofewers,quote=F); cat("\n")
      to.drop <- unique(c(to.drop,toofewers))
    } else { cat(" samples from plates with too few noted but left in the dataset\n\n")}
  } else {
    cat(" NB: reporting thresholds only: Plate CNV QC switched off\n")
    toomanyers <- toofewers <- to.drop <- character(0)
  }
  remove.samp.from.plink(to.drop,cat.path(dir$cnv.qc,rem.file),why=suf.for.rmv,dir=dir)
  cat(" done\n")
  return(to.drop)
}


filter.hifreq.cnv.samples <- function(plink.file, nsamps, rem.file=plink.file, rare="", dir, thr.sd=3,
                                      ndels=NULL, ndups=NULL, del.rate=0.4, dup.rate=0.18, QC.ON=T, pval=0.05)
{
  # filter samples with too many rare DELs/DUPs (according to specified rare rates) (poisson)
  # or filter samples with too many CNVs (normal)
  dir <- validate.dir.for(dir,c("cnv.plk"),warn=F)
  if(is.null(del.rate) & is.null(ndels) | is.null(dup.rate) & is.null(ndups)) {
    stop("Error: must specify at least 1 of dup/del.rate or ndups/ndels")
  }
  rare <- toupper(rare)
  if(!is.null(ndels)) { del.rate <- ndels/nsamps } else { ndels <- round(nsamps*del.rate) }
  if(!is.null(ndups)) { dup.rate <- ndups/nsamps } else { ndups <- round(nsamps*dup.rate) }
  if(rare %in% c("DEL","DUP")) { suf.for.rmv <- paste("_",rare,"rcnt",sep="") } else { suf.for.rmv <- "_CNVcnt" }
  # establish cutoffs for rare deletions/duplications using rates in all samples *2
  dup.cutoff <- which(round(ppois(c(1:10),dup.rate,lower.tail=F)*ndups,4)<pval)[1]-1
  del.cutoff <- which(round(ppois(c(1:10),del.rate,lower.tail=F)*ndels,4)<pval)[1]-1

  dat <- read.table(cat.path(dir$cnv.qc,plink.file,ext="cnv"),header=T)  #; colnames(dat) <- 
  tt <- as.numeric(table(dat$IID))[order(as.numeric(table(dat$IID)))]
  ss <- names(table(dat$IID))[order(as.numeric(table(dat$IID)))]
  if(!rare %in% c("DEL","DUP")) {
    rtxt <- ""
    # mean/sd/'3'SD for dels/dups per subj
    thr <- round(sd(tt)*thr.sd + mean(tt))
    cat("\nFiltering samples with a large number of CNVs\n")
    cat("\n distribution: mean=",round(mean(tt),2),", SD=",round(sd(tt),2),"   ==>   threshold=",
        ceiling(thr),"   [mean + ",thr.sd,"*sd]\n",sep="")
  } else { 
    rtxt <- paste("[",rare,"s]",sep="");
    cat("\nFiltering samples with a large number of rare",rare,"s\n",sep="")
    if(rare=="DUP") { cat(" cutoff for duplicates: ",
                      dup.cutoff,"; p<",pval,", X~poisson(",round(dup.rate,2),")\n",sep="") }
    if(rare=="DEL") { cat(" cutoff for deletions: ",
                      del.cutoff,"; p<",pval,", X~poisson(",round(del.rate,2),")\n",sep="") }
  }
  cat("\nCNV length stats",rtxt,"\n")
  lenz <- (dat$BP2 - dat$BP1)
  print(summary(lenz)) #which(lenz==10048355,"\n")
  cat("\n")
  do.lens.summary(lenz)
  ord.list <- cbind(ss,tt)[order(tt),] # ordered lists of dels/dups by SubID
  colnames(ord.list) <- c("SampleIDs","CNVs")
  cat("\n30 highest CNV samples",rtxt,"\n")
  if(!is.null(dim(ord.list))) {
    if(nrow(ord.list)>30) {
      print(tail(ord.list,30),quote=F)
    } else {
      print(ord.list)
    }
  } else {
    print("none found!\n")
  }
  cat("\nHISTO of CNVs per sample",rtxt,"\n")
  print(textogram(as.numeric(ord.list[,2])))
  if(QC.ON) {
    if (rare=="DEL") { thr <- del.cutoff  }
    if (rare=="DUP") { thr <- dup.cutoff  }
    toomanyers <- ss[(which(tt>thr))]
    cat(paste("\n list of",length(toomanyers),"samples with too many CNVs from",plink.file),":\n")
    print(toomanyers,quote=F); cat("\n")
  } else {
    cat(" NB: reporting thresholds only: CNV QC",rtxt,"switched off\n")
    toomanyers <- character(0)
  }
  remove.samp.from.plink(toomanyers,cat.path(dir$cnv.qc,rem.file),why=suf.for.rmv,dir=dir)
  cat(" done\n")
  return(toomanyers)
}


change.dir.data.tracker <- function(DT,orig.dir,new.dir) {
  # if wanting to change the directory for a session, you must
  # edit the directory in many parts.
  if(is.list(DT)) {
    if(length(DT)>0) {
      for (aa in 1:length(DT)) {
        if(is.list(DT[[aa]])) {
          DT[[aa]] <- change.dir.data.tracker(DT[[aa]],orig.dir,new.dir)
        } else {
          typ <- is(DT[[aa]])[1]
          if(is.null(dim(DT[[aa]]))) {
            DT[[aa]] <- gsub(orig.dir,new.dir,DT[[aa]],fixed=T)
            DT[[aa]] <- as(object=DT[[aa]],typ) # otherwise all get converted to character
          } else {
            # assume matrices won't have directories to change
          }
        }
      }
    }
  }
  return(DT)  
}


blank.data.tracker <- function(dir=NULL,n.grps=1,grp.names=NULL,description=list())
{
  ## setup a blank data.tracker object, optionally with some directory info
  DATATRACKER <- list()
  lgn <- length(grp.names)
  if(lgn>0 & n.grps==1) { n.grps <- lgn }
  blank1 <- list(); blankn <- vector("list",n.grps); blankv <- numeric(n.grps); 
  blank8 <- vector("list",8); names(blank8) <- c("init",paste(0:6))
  if(lgn==n.grps) { names(blankn) <- grp.names }
  bFIELDS <- NULL #data.frame(row.names=names(blankn),sample=2,snp=1,a1=3,a2=4)
  bIDs <- list(samples=character(),snps=character())
  bINFO <- list(snp.info=NULL,sample.info=NULL,gc.genome=NULL,gc.marker=NULL)
  DATATRACKER[[1]] <- list(DESCRIPTION=description,BASE=dir$base,DIR=dir,N.GRPS=n.grps,
                           GRP.NAMES=grp.names,MAP.FILE=NULL,MAP.TYPE=NULL,
                           N.FILES=blankv,PROC.DONE=NULL,SETTINGS=blank8)
  DATATRACKER[[2]] <- list(LRR=blankn,BAF=blank1,SETUP=cat.path(dir$lrr.dat,"file.spec.txt"))
  DATATRACKER[[3]] <- list(GENO=blank1,SNPMAT=blank1,FIELDS=bFIELDS)
  DATATRACKER[[4]] <- list(LRR=blankn,BAF=blank1,IDs=bIDs)
  DATATRACKER[[5]] <- list(LRR=blankn,SAMPLES=blankn,MEDIAN=blankn,STATS=blankn)
  DATATRACKER[[6]] <- list(LRR=blankn)
  DATATRACKER[[7]] <- list(LRR=blankn,EIGEN=blankn)
  DATATRACKER[[8]] <- list(LRR=NULL,BAF=NULL,INFO=bINFO) 
  DATATRACKER[[9]] <- list(LRR=blank1,EIGEN=blank1,MEDIAN=blank1,STATS=blank1) # [min.n.pcs ... max.n.pcs]
  DATATRACKER[[10]] <- list(RESULT=blank1) # [name1,...]
  names(DATATRACKER) <- c("Info","RawDataGrps", "SnpDataGrps", "ShapedDataGrps", "BigRawDataGrps",
                          "BigFiltDataGrps","BigPCADataGrps","BigQCData", "BigCorrectedData",
                          "CnvResults")
  return(DATATRACKER)
}


init.data.tracker <- function(dir,grps=NA,grp.names=NULL,raw.lrr=NULL,raw.baf=NULL,
                              geno.file=NULL,shape.lrr=NULL,shape.baf=NULL,
                              big.lrr=NULL,big.baf=NULL,ucsc="hg18",settings=NULL,
                              snps="snpNames.txt",samples="subIdsALL.txt",sub.samples=NULL,
                              sample.info=NULL,snp.info=NULL,snp.fields=NULL,description=NULL,
                              file.spec="file.spec.txt",map.file="snpdata.map",map.type="map3",
                              plate.fn="plate.lookup.txt",pheno.fn="pheno.lookup.txt",sex.fn="sex.lookup.txt",
                              group.fn=NULL,other.batch=list(),verbose=F,scheme=NA,col1="black",col2="grey")
{
  # Initialise a new DATATRACKER object
  dir <- validate.dir.for(dir,c("base","ano","col","ids","raw"))
  # Info - DESCRIPTION - BASE - DIR
  # description can contain anything you like, e.g, text, or some settings, etc
  varlist <- c("samples","snps","map.file","plate.fn","plate.fn")
  for (cc in 1:length(varlist)) {
    assign(x=varlist[cc],value=add.dir.if.not(get(varlist[cc]),dir$ano,dir)) # add 'dir$ano' if needed
  }
  file.info <- get.file.specs(dir,fn=file.spec)
  if(!is.null(file.info)) {
    if("GRP" %in% toupper(colnames(file.info))) {
      if( (all(is.na(grps))) | (any(grps!=file.info$GRP)) | (length(grps)!=length(file.info$GRP)) ) {
        grps <- (file.info$GRP)
        cat(" changed specified 'grps' in data.tracker to match file.spec.txt\n")
      }
    }
  } else {
    cat(" NB: file.spec.txt not in use, must enter groups as a vector of the same length as number of initial data files\n")
  }
  num.grps <- length(unique(grps))
  ## try to make best guess at number of groups if 'grps' is NA as this will affect datatracker structure
  if(all(is.na(grps))) {
    if(!is.null(sub.samples) | !is.null(raw.lrr) | !is.null(shape.lrr) ) {
      num.grps <- max(length(sub.samples),length(raw.lrr),length(shape.lrr))
    } else {
      num.grps <- 1
      cat("grps is NA, and no datafiles referenced, so defaulting to datatracker expecting 1 file, 1 group\n")
    }
  }
  DATATRACKER <- blank.data.tracker(dir,n.grps=num.grps,grp.names=grp.names,description)
  if(!is.null(samples) & is.null(sample.info) & !is.null(plate.fn)) {
    sample.info <- make.sample.info(dir=dir,samples,phenotype=pheno.fn,plate=plate.fn,
                                    grp=group.fn,sex=sex.fn,other=other.batch)
  } else {
    cat("Warning, sample.info object not created as sample list",samples,", or plate lookup file",plate.fn,"were missing!\n")
  }
  if(!is.null(snps) & is.null(snp.info) & is.character(map.file)) {
    #prv(snps,map.type,map.file,scheme,col1,col2)
    snp.info <- make.snp.info(dir,snpIDs=snps,anot=map.type,snp.info.fn=map.file,make.sort=T,
                              non.autosomes.excl.file=T,verbose=verbose,ucsc=ucsc,scheme=scheme,col1=col1,col2=col2)
  } else {
    cat("Warning, snp.info object not created as the",snps,", or the",map.file,"[chr,pos,id] file were missing!\n")
  }
  sample.info.fn <- write.sample.info(sample.info,dir)
  snp.info.fn <- write.snp.info(snp.info,dir)
  file.info <- get.file.specs(dir,fn=file.spec)
  ## if we have file.spec.txt then many fields can be nicely populated automatically
  if(!is.null(file.info)) {
    if(all(c("FILE","GRP") %in% colnames(file.info))) {
      raw.lrr <- cat.path(dir$raw,file.info$FILE)
      if(all(c("LRR","BAF") %in% colnames(file.info))) {
        if(is.null(raw.baf) & all(!is.na(as.numeric(file.info$BAF)))) { raw.baf <- raw.lrr }
        DATATRACKER$Info$N.FILES <- as.numeric(table(file.info$GRP))
      }
      field.hdrs <- c("SAMP","SNP","A1","A2")
      if(all(field.hdrs %in% colnames(file.info))) {
        snp.fields <- file.info[,field.hdrs]
        if(is.null(geno.file) & all(!is.na(as.numeric((snp.fields[,3]))))) { 
          geno.file <- raw.lrr # assume genotype data is in genome studio file (usually is)
        }
      } 
      if(is.null(sub.samples)) {
        sub.samples <- cat.path(dir$ids,get.subIDs(dir,"files"))
      }
      if(is.null(shape.lrr)) {
        shape.lrr <- get.lrr.files(dir) # will return null if not yet present 
      }
      if(is.null(shape.baf)) {
        shape.baf <- get.lrr.files(dir,suffix="BAF.dat",baf=T)
      }
    }
  }
  ## without file.spec.txt do our best to fill in the essential fields of the data.tracker
  if(is.null(file.info)) {
    DATATRACKER$RawDataGrps$SETUP <- NULL # remove reference to file.spec.txt
    if(is.null(sub.samples)) { 
      cat(" 'sub.samples' is NULL and file.spec.txt not in use. \n default is contents of",dir$ids,"\n")
      sub.samples <- cat.path(dir$ids,list.files(dir$ids))
    }
    n.dats <- length(sub.samples); n.grps <- length(grps)
    if(all(is.na(grps))) {
      # if 'grps' not specified, try to make the best guess #
      if(n.dats>0) { 
        grps <- 1:n.dats 
      } else {
        shape.n <- (length(list.files(dir$col)))
        raw.n <- (length(list.files(dir$raw)))
        if( is.file(samples,dir$ano,dir) & 
          ( (shape.n==1) | (raw.n==1) ) ) {
          grps <- 1; sub.samples <- samples
        } else {
          stop(paste("'grps' was NA, sub.samples was empty, number of files in data",
            "directories was: ",raw.n,shape.n,", and file.spec.txt not in use, so cannot proceed"))
        }
      }
    }
    if(all(n.dats==grps) & n.dats>1 & n.grps==1) {
      # e.g, mistakenly entered '3' instead of '1:3'
      cat(" grps should be a vector same length as sub.samples [",n.dats,
          "], was ",grps," so setting grps to:",paste(1:grps,collapse=","),"\n",sep="")
      grps <- 1:grps
    }
    if((n.grps==1) & (all(grps==1)) & (n.dats>1) ) {
      cat(" multiple datafiles detected, grp=1 specified, assuming all to be treated as a single cohort\n")
      grps <- rep(1,n.dats)
    }
    if((n.grps>1) & (n.grps==n.dats)) {
      ## seemingly successful!
      cat(" found",n.dats,"sample id files and ",n.grps," values for the 'grps' parameter\n")
    }
    # set geno.file if genotypes seem to be in raw.lrr file
    if(is.null(geno.file) & all(!is.na(as.numeric((snp.fields[,3]))))) { geno.file <- raw.lrr }
    # warn about missing snp.fields if it seems like they want it
    if(is.null(snp.fields) & (!is.null(raw.lrr) | !is.null(geno.file))) { 
      warning("snp.fields is NULL but there seems to be data in either ",
          "geno.file or raw.lrr, implying SNP-QC should be run, ",
          "so recommend a restart specifying fields for SNP,SAMP,A1,A2 in geno.file")
    }
    ## whether or not grps was specified, run with what we have and guess n.files values #
    cat(" group lengths:",paste(as.numeric(table(grps)),collapse=","),"\n")
    DATATRACKER$Info$N.FILES <- as.numeric(table(grps))
  }
  if(is.null(geno.file)) { cat(" 'geno.file' is NULL, assuming SNP-QC will be skipped\n") }
  DATATRACKER <- setSlot(DT=DATATRACKER,dir=dir,grps=grps,raw.lrr=raw.lrr,raw.baf=raw.baf,
                                     shape.lrr=shape.lrr,shape.baf=shape.baf,big.lrr=big.lrr,settings=settings,
                                     sub.samples=sub.samples,snps=snps,samples=samples,geno.file=geno.file,
                                     snp.fields=snp.fields,sample.info.fn=sample.info.fn,snp.info.fn=snp.info.fn,
                                     big.baf=big.baf,description=description,map.file=map.file,map.type=map.type)
  # RawDataGrps - LRR - BAF
  # SnpDataGrps - LRR - FIELDS [samp,snp,a1,a2]
  # ShapedDataGrps - LRR - BAF - IDs [samples,snps]
  # BigRawDataGrps - LRR - IDs [samples]
  # BigFiltDataGrps - LRR
  # BigPCADataGrps - LRR
  # BigQCCombData - LRR - BAF - INFO [sample.info, snp.info]
  # BigCorrectedData - LRR [min.n.pcs ... max.n.pcs]
  # CnvResultData - RESULT [name1...]
  return(DATATRACKER)
}


setSlot <- function(DT,dir=NULL,grps=NA,n.pcs=NA,raw.lrr=NULL,raw.baf=NULL,
                                shape.lrr=NULL,shape.baf=NULL,big.lrr=NULL,big.baf=NULL,big.filt=NULL,
                                pcas=NULL,eigens=NULL,sub.samples=NULL,snps=NULL,samples=NULL,
                                geno.file=NULL,snp.fields=NULL,sample.info.fn=NULL,snp.info.fn=NULL,
                                big.qc=NULL,snp.stats=NULL,
                                pca=NULL,eigen=NULL,cnv.result=NULL,result.name=NULL,
                                description=NULL,map.file=NULL,map.type=NULL,
                                medians=NULL,median=NULL,gc.genome=NULL,gc.marker=NULL,
                                stats=NULL,stat=NULL,proc.done=NULL,warns=F,settings=NULL)
{
  # update some field(s) in a datatracker object
  # put filenames etc into groups according to parameter 'grps'
  # if 'grps' not inputted, use other information to infer if shape.lrr is inputted
  # only works if input is sorted by grp e.g, 1,1,2,3,3,3,3,4, etc
  if(!is.data.tracker(DT)) {
    stop("Error: could not update because 'DT' was not a data.tracker object")
  }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,c("base","raw","ano","big","col","ids"))
    if(dir$base!=DT$Info$BASE) {
      warning("directory object did not match existing. No change to data tracker object")
      return(DT)
    }
  }
  if(all(is.na(grps))) {
    ## if 'grps' not inputted, use other information to infer if shape.lrr is inputted
    # only works if input is sorted by grp e.g, 1,1,2,3,3,3,3,4, etc
    if(!is.null(shape.lrr) | !is.null(sub.samples)) {
      exgrps <- DT$Info$N.GRPS
      if(is.numeric(exgrps)) {
        ncheck <- NA; 
        if(is.numeric(DT$Info$N.FILES)) { if(all(DT$Info$N.FILES>0)){ 
          ncheck <- DT$Info$N.FILES 
          grps <- rep(1:exgrps,times=ncheck)
          if(warns) { warning(paste("'grps' was NA, inferred group order of input:",paste(grps,collapse=","))) }
        } }
      }
    }
    if(all(is.na(grps)) & warns) { warning("without specifying 'grps', filenames might not be placed in tracker") }
  }
  if(all(!is.na(grps))) {
    #check grp based input
    exgrps <- DT$Info$N.GRPS
    nmgrps <- DT$Info$GRP.NAMES
    ugrps <- unique(grps)
    nugrps <- length(ugrps)
    if(all(ugrps %in% nmgrps)) { 
      grps <- match(grps,nmgrps); ugrps <- unique(grps)
    }
    if(max(ugrps)>exgrps) {
      stop("Error: illegal group number specified in 'grps', exceeds number of grps")
    }
    if(nugrps<exgrps) {
      if(warns) { warning("fewer grps specified than 'DT$Info$N.GRPS' so assuming update applies to subset") }
    } else {
      if(nugrps!=exgrps) {
        stop(paste("Error: number of grps specified (",nugrps,
                   ") exceeds n.grps in 'DT' (",exgrps,")"),sep="")
      }
    }
    for (cc in ugrps) {
      ## raw files (length 'grps')
      if(is.ch(raw.lrr) & length(raw.lrr)==length(grps))
      { DT$RawDataGrps$LRR[[cc]] <- raw.lrr[which(grps==cc)] }
      if(is.ch(shape.lrr) & length(shape.lrr)==length(grps)) 
      { DT$ShapedDataGrps$LRR[[cc]] <- shape.lrr[which(grps==cc)] }
      if(is.ch(geno.file) & length(geno.file)==length(grps))
      { DT$SnpDataGrps$GENO[[cc]] <- geno.file[which(grps==cc)] }
      if(is.ch(sub.samples) & length(sub.samples)==length(grps))
      { DT$BigRawDataGrps$SAMPLES[[cc]] <- sub.samples[which(grps==cc)] }
      if(is.ch(raw.baf) & length(raw.baf)==length(grps))
      { DT$RawDataGrps$BAF[[cc]] <- raw.baf[which(grps==cc)]  } 
      if(is.ch(shape.baf) & length(shape.baf)==length(grps))
      { DT$ShapedDataGrps$BAF[[cc]] <- shape.baf[which(grps==cc)]  }
      # processed files (length 'ugrps')
      if(is.ch(big.lrr) & length(big.lrr)==length(ugrps))
      { DT$BigRawDataGrps$LRR[[cc]] <- unlist(big.lrr)[which(ugrps==cc)] }
      if(is.ch(big.filt) & length(big.filt)==length(ugrps))
      { DT$BigFiltDataGrps$LRR[[cc]] <- unlist(big.filt)[which(ugrps==cc)] } 
      if(is.ch(pcas) & length(pcas)==length(ugrps))
      { DT$BigPCADataGrps$LRR[[cc]] <- unlist(pcas)[which(ugrps==cc)] }
      if(is.ch(eigens) & length(eigens)==length(ugrps))
      { DT$BigPCADataGrps$EIGEN[[cc]] <- unlist(eigens)[which(ugrps==cc)] }
      if(is.ch(medians) & length(medians)==length(ugrps))
      { DT$BigRawDataGrps$MEDIAN[[cc]] <- unlist(medians)[which(ugrps==cc)] }
      if(is.ch(stats) & length(stats)==length(ugrps))
      { DT$BigRawDataGrps$STATS[[cc]] <- unlist(stats)[which(ugrps==cc)] }
    }
  }
  # update baf grp things
  # BAF raw/shape can be grouped or ungrouped, if not already found to be group, try singular
  if(all(is.na(grps))) {
    ## not 100% sure whether it should be with the [[1]]'s or not below:...
    if(is.ch(raw.baf) & length(raw.baf)==1) 
    { DT$RawDataGrps$BAF[[1]] <- raw.baf }
    if(is.ch(shape.baf) & length(shape.baf)==1) 
    { DT$ShapedDataGrps$BAF[[1]] <- shape.baf } 
    if(is.ch(geno.file) & length(geno.file)==1) 
    { DT$SnpDataGrps$GENO[[1]] <- geno.file }
  }
  # update combined (unique) inputs
  if(is.ch(map.file)) { DT$Info$MAP.FILE <- map.file }
  if(is.ch(map.type)) { DT$Info$MAP.TYPE <- map.type }
  if(is.ch(snp.stats)) { DT$SnpDataGrps$SNPMAT <- snp.stats }
  if(is.ch(description)) { DT$Info$DESCRIPTION <- description }
  if(is.ch(big.qc)) { DT$BigQCData$LRR <- big.qc }
  if(is.ch(big.baf)) { DT$BigQCData$BAF <- big.baf }
  if(is.ch(gc.genome)) { DT$BigQCData$INFO$gc.genome <- gc.genome }
  if(is.ch(gc.marker)) { DT$BigQCData$INFO$gc.marker <- gc.marker }
  if(is.ch(sample.info.fn)) { DT$BigQCData$INFO$sample.info <- sample.info.fn }
  if(is.ch(snp.info.fn)) { DT$BigQCData$INFO$snp.info <- snp.info.fn }
  if(is.ch(snps)) { DT$ShapedDataGrps$IDs$snps <- snps }
  if(is.ch(samples)) { DT$ShapedDataGrps$IDs$samples <- samples }
  # add settings to correspond to process just done - or 'init' if no process specified
  if(is.list(settings)) { 
    if(is.null(proc.done)) {
      DT$Info$SETTINGS[[1]] <- settings 
    } else { 
      pD <- as.numeric(paste(proc.done))
      if(pD %in% 0:6) {
        DT$Info$SETTINGS[[pD+2]] <- settings
      }
    }
  }
  # update number of grps
  get.rw.lrr <- getSlot(DT,"raw.lrr")
  get.sh.lrr <- getSlot(DT,"shape.lrr")
  bst <- get.rw.lrr; if(length(get.sh.lrr)>length(bst)) { bst <- get.sh.lrr }
  if(length(bst)>0) { 
    DT$Info$N.FILES <- sapply(bst,length) 
  }
  if(length(proc.done)>0) {
    if(all(proc.done %in% c(0:6))) {
      DT$Info$PROC.DONE <- as.integer(unique(c(DT$Info$PROC.DONE,proc.done)))
    } else {
      warning("invalid value entered for 'proc.done', this process won't be registered")
    }
  }
  # update snp fields matrix; needs one entry per file, all groups together, must be inputted all grps together
  if(!is.null(dim(snp.fields))) { if(dim(snp.fields)[1]==max(DT$Info$N.FILES,length(grps)))
  { DT$SnpDataGrps$FIELDS <- snp.fields } }
  # update named inputs
  if(all(!is.na(n.pcs))) {
    for(dd in 1:length(n.pcs)) {
      if(is.ch(pca)) { DT$BigCorrectedData$LRR[[paste("pc",n.pcs[dd],sep="")]] <- pca[[dd]] }
      if(is.ch(eigen)) { DT$BigCorrectedData$EIGEN[[paste("pc",n.pcs[dd],sep="")]] <- eigen[[dd]] }
      if(is.ch(median)) { DT$BigCorrectedData$MEDIAN[[paste("pc",n.pcs[dd],sep="")]] <- median[[dd]] }
      if(is.ch(stat)) { DT$BigCorrectedData$STATS[[paste("pc",n.pcs[dd],sep="")]] <- stat[[dd]] }
      if(is.ch(cnv.result)) { DT$CnvResults$RESULT[[paste("pc",n.pcs[dd],sep="")]] <- cnv.result[[dd]] }
    }
  }
  #if(is.character(result.name) & is.ch(cnv.result)) {
   # for(dd in 1:length(result.name)) {
     # DT$CnvResultData$RESULT[[result.name]] <- cnv.result[[dd]]
   # }
  #}
  if(is.data.tracker(DT)) {
    return(DT)
  } else {
    stop("Error: data tracker object became corrupted during update")
  }
}


methods.data.tracker <- function(all=T) {
  dats <- c("map.file","sample.info.fn","snp.info.fn","snp.stats",
            "snp.fields","sample.info","gc.genome","gc.marker","stat","stats")
  ranged <- c("cnv.results","cnv.result","cnv.calls","snp.info")
  bigs <- c("big.lrr","big.filt","big.qc","big.baf","pcas","big.pcc")
  massives <- c("geno","genofile","raw.lrr","raw.baf", "shape.lrr","shape.baf")
  scalars <- c("n.grps","num.grps","n.grp","map.type")
  vec.not.files <- c("group.names","description","n.files","num.files","n.file","proc.done")
  vec.files <- c("eigens","sub.samples","snps","samples","eigen")
  vecs <- c(vec.not.files,vec.files)
  dats.not.files <- c("snp.fields","settings")
  numz <- c("n.grps","num.grps","n.grp","n.files","num.files","n.file")
  strz <- c("group.names","description","map.type")
  fnames <- c(dats,bigs,massives,vec.files,ranged)
  valid <- c(dats,bigs,massives,scalars,vecs,ranged,"dir")
  if(all) {
    out <- list(valid,dats,ranged,bigs,massives,scalars,vecs,vec.files,vec.not.files,dats.not.files,numz,strz,fnames)
    names(out) <- c("valid","dats","ranged","bigs","massives","scalars","vecs","vec.files",
                    "vec.not.files","dats.not.files","numbers","strings","filenames")
    return(out)
  } else {
    return(valid)
  }
}


req.valid <- function(DT,plist=c(),grps=NA,dir=NULL,n.pcs=NA) {
  ## ensure that required fields are present and valid in a datatracker object
  if(!is.data.tracker(DT)) {
    warning("not a data.tracker object")
    return(NULL)
  }
  stdzr <- function(x) { tolower(gsub(".","",x,fixed=T)) } 
  mdt <- methods.data.tracker(T)
  numz <- mdt$numbers; strz <- mdt$strings; dnf <- mdt$dats.not.files
  if(is.null(dir)) { dir <- getSlot(DT,"dir") }
  is.val <- T
  if(length(plist)<1 | !is.character(plist)) 
  { warning("invalid list length or type"); return(NULL) }
  badlist <- plist[!stdzr(plist) %in% stdzr(mdt$valid)]
  if(length(badlist)>0)
  { warning(paste(badlist,collapse=",")," are not data.tracker methods") }
  plist <- plist[stdzr(plist) %in% stdzr(mdt$valid)] #only test valid methods
  for(cc in 1:length(plist)) {
    next.thing <- getSlot(DT,what=paste(plist[cc]),grps=grps,n.pcs=n.pcs)
    if(stdzr(plist[cc]) %in% stdzr(numz)) {
      nxt.val <- is.numeric(next.thing)
    } else {
      if(stdzr(plist[cc]) %in% stdzr(strz)) {
        nxt.val <- is.ch(next.thing)
      } else {
        if(stdzr(plist[cc]) %in% stdzr(dnf)) {
          nxt.val <- !is.null(next.thing)
        } else {
          if(is.null(next.thing) | length(next.thing)==0) {
            cat(" invalid DT object:",plist[cc],";",next.thing,"\n") } else {
          nxt.val <- is.file(next.thing,dir) }
        }
      }
    }
    if(!nxt.val) { warning(plist[cc]," not valid") }
    is.val <- (is.val & nxt.val)
  }
  return(is.val)
}


three.way.comparison <- function(clean.stats,sample.info,raw.fn="StatsPerSample",dir,
                                 batch.comps=c("plate","grp"),bxplot=T,pref="LRR",raw.stat.fn=NULL,...) 
{
  ## function assuming data is PC-clean and raw stats per sample summary files are
  # in the default location with default names. if so, makes nice visualisation
  dir <- validate.dir.for(dir,c("qc.pc"),warn=F)
  grpnumz <- as.numeric(unique(sample.info$grp))
  raw.stat.table <- combine.raw.stats(dir=dir,grpnumz=grpnumz,base.fn=raw.fn,pref=pref,raw.stat.fn=raw.stat.fn)  # raw.stat fn can come from DT
  sampqc.stat.table <- remove.qc.failers(raw.stat.table,sample.info,dir=dir)
  three.part.comp <- list(raw.stat.table,sampqc.stat.table,clean.stats)
  names(three.part.comp) <- c("Raw","Sample.QC","PC.Corrected")
  vars <- colnames(clean.stats)
  vars <- vars[vars %in% colnames(raw.stat.table)]
  if(length(vars)<2) { return(three.part.comp) } # don't bother w/ plot if not at least 2 stats 
  ### ENSURE EACH STAT TABLE HAS THE SAME COMMON COLUMN NAMES
  three.part.common <- three.part.comp
  for (jj in 1:3) {
    three.part.common[[jj]] <- 
      as.data.frame(three.part.common[[jj]][,match(vars,colnames(three.part.common[[jj]]))])
  }
    if(bxplot) {
    for (dd in 1:length(batch.comps)) {
      ofn <- cat.path(dir$pc,"ThreeStageComparison",suf=batch.comps[dd],ext=".pdf")
      pdf(ofn,...)
      par(mfcol=c(ncol(three.part.common[[3]]),length(three.part.common)))
      # some batch effects diagnostic plots to see how well the correction has worked
      for (cc in 1:length(three.part.common)) {
        batch.box.plots(three.part.common[[cc]],NULL,lookup=sample.info,batch=batch.comps[dd],
                        pref=paste(batch.comps[dd],names(three.part.common)[cc],sep="_"),
                        subtitle=paste("[",names(three.part.common)[cc],"]"),dir=dir$pc,to.file=F)
      }
      dev.off()
      cat("~wrote plots to:",ofn,"\n")
    }
  }
  return(three.part.comp)
}


combine.raw.stats <- function(dir,grpnumz=NA,base.fn="StatsPerSample",pref="LRR",raw.stat.fn=NULL) 
{
  # combine set of sample statistics from multiple files (cohorts) into 1 big file
  # or in the case of length grpnumz = 1, simply return the only file.
  # can deduce in 3 ways: i) user specified in run.PCA.correct; ii) from DT; iii) look for expected number of grps
  
  # automatically try to work out the preferences if file.spec.txt was there  #DEP?
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  ## an exact location of stats file might have been passed from the datatracker (or user specified in run.PCA.co..)
  if(is.null(raw.stat.fn)) {
    got.fn <- F 
  } else {
    if(is.file(raw.stat.fn,dir$qc.lrr,dir)) { got.fn <- T } else { got.fn <- F }
  }
  if(got.fn){
    ## file names passed down from data.tracker parameter 'stats'
    if(all(!is.na(grpnumz))) {
      if(length(grpnumz)!=length(raw.stat.fn)) { 
        warning(paste("number of files [",length(raw.stat.fn),"] didn't match group numbers [",length(grpnumz),"]\n",sep=""),
        "Note that if the files are correct, this function will work despite conflicting group count")
      }
    }
    stat.table <- list()
    for (kk in 1:length(raw.stat.fn)) {
      ifn <- find.file(raw.stat.fn[kk],dir$qc.lrr,dir)
      stat.table[[kk]] <- read.table(ifn)
    }
    stat.table <- do.call("rbind",stat.table)
  } else {
    if(all(is.na(grpnumz))) { warning("need parameter 'grpnumz' for combine.raw.stats(), failure likely") }
    if(length(grpnumz)==1) {
      ## only 1 group so looking for one file potentially without a suffix
      ifn <- cat.path(dir$qc.lrr,base.fn,pref=pref,ext=".tab")
      ifn2 <- cat.path(dir$qc.lrr,base.fn,pref=paste(pref,"1",sep=""),ext=".tab")
      if(!file.exists(ifn)) { 
        if(file.exists(ifn2)) {
          ifn <- ifn2
        } else {
          in.dir <- list.files(dir$qc.lrr)
          tryanything <- grep(base.fn,in.dir,fixed=T)
          if(length(tryanything)>0) {
            warning(paste("Didn't find expected sample-wise statistics text file in\n",dir$qc.lrr,"\n"),
            paste("but did find:",in.dir[tryanything][1],"\nattempting to use as stats file. \n"),
            paste("stats file with 1 grp should have been named:\n",ifn))
            ifn <- cat.path(dir$qc.lrr,in.dir[tryanything[1]])
          } else {
            Emsg=paste("Error: could not find any sample-wise stats file named *",base.fn,"* in:\n ",dir$qc.lrr,"\n",sep="")
            stop(Emsg)
          }
        }
      }
      stat.table <- read.table(ifn)
    } else {
      # multiple groups so looking for files with suffixes
      stat.table <- list()
      for (kk in grpnumz) {
        ifn <- cat.path(dir$qc.lrr,base.fn,pref=pref[kk],ext=".tab")
        if(file.exists(ifn)) {
          stat.table[[kk]] <- read.table(ifn)
        } else {
          warning(paste("Didn't find expected samplewise stats file:",ifn))
        }
      }
      stat.table <- do.call("rbind",stat.table)
    }
  }
  return(stat.table)
}


batch.box.plots <- function(plot.stats,pass.tab=NULL,lookup,batch="plate",pref="All",subtitle="",dir="",to.file=T)
{
  ## plot box plot of batch effect, using stats summary and looking up batch categories from table
  batch <- paste(batch[1])
  index <- find.id.col(lookup,ids=rownames(plot.stats))$index
  if(!(any(batch %in% colnames(lookup))))
  { cat("'lookup' expected to contain '",paste(batch,collapse=","),"' and 'id' named columns. No boxplot(s) produced.",sep=""); return(NULL)}
  batch <- batch[batch %in% colnames(lookup)] # at least 1 is present, remove any that aren't
  if(is.matrix(lookup)) { lookup <- as.data.frame(lookup) }
  if(is.matrix(plot.stats)) { plot.stats <- as.data.frame(plot.stats) }
  plts <- paste(unique(lookup[,batch]))
  pltnums <- 1:length(plts)
  samp.plts <- match(lookup[,batch],plts)
  plt.samp <- samp.plts[index]  #match(rownames(plot.stats),lookup[,"id"])]
  statz <- colnames(plot.stats) # e.g, c("Mean","DLRS","GCWave")
  ylz <- list( c(-.2,.2), c(0,.6) , c(-.2,.2) )
  spotz <- c("topleft","topright")
  if(to.file) {
    dir <- validate.dir.for(dir,c("qc.pl"),warn=F)
    ofn <- cat.path(dir$qc.pl,pre=batch,fn="BoxPlotsMDG",suf=pref,ext="pdf")
    pdf(ofn)
  }
  for (tt in 1:length(statz))
  { 
    stat <- paste(statz[tt])
    yl <- ylz[[tt]]
    spott <- spotz[tt]
    if(length(plot.stats[[stat]])>0 & length(plot.stats[[stat]])==length(factor(plt.samp))) 
    {
      boxplot(plot.stats[[stat]]~factor(plt.samp),xlab=paste(batch,"number"),
              ylab=paste("LRR",stat), lwd=1, pch=".", bty="l",
              main=paste("LRR",stat,"distribution"),sub=subtitle,ylim=yl,bty="l")
      # use segments here instead of abline to avoid overshadowing some boxes
      if(!is.null(pass.tab)) {
        abline(h=pass.tab[c("LB","UB"),stat],col="lightblue") #,lty="dotted")
        legend("top",legend=c("Lower and upper bounds"),
               lty="solid",col="lightblue",bty="n")
      }
    } else {
      warning(paste("skipped",stat,"in boxplot due to missing data"))
    }
  }
  if(to.file) {  dev.off() ; cat(paste("~wrote file:",ofn,"\n")) }
}


get.extreme.examples <- function(stat.table,failures,venn.lists)
{
 # find most extreme subject for each part of venn list
 # [NB: or for 'none' group - the least extreme ]
 last.n <- length(venn.lists)
 z.scores <- apply(stat.table,2,standardize)
 sample.dev.scores <- rowSums(abs(z.scores)*(failures)) # scores to get MOST extreme
 none.grp <- match(venn.lists[[last.n]],names(sample.dev.scores)) # group with no outliers
 sample.dev.scores[none.grp] <- rowSums(abs(z.scores))[none.grp] # scores to get LEAST extreme

 ex.id <- integer(last.n)
 # find most extreme subject for each part of venn list

 for (cc in 1:(last.n))
 {
   samp.nos <- match(venn.lists[[cc]],names(sample.dev.scores)) 
   # find max for 1-8, min for 9 [none group]
   if (length(samp.nos)>0) {
     if(cc<length(venn.lists))
     {
       mx <- which(sample.dev.scores[samp.nos]==max(sample.dev.scores[samp.nos],na.rm=T))
     } else {
       mx <- which(sample.dev.scores[samp.nos]==min(sample.dev.scores[samp.nos],na.rm=T))
     }
     ex.id[cc] <- names(sample.dev.scores)[samp.nos[mx[1]]]
   } else { ex.id[cc] <- NA }
 }
 names(ex.id) <- names(venn.lists)   
 return(ex.id)
}


col.extremes <- function(lat,selBlue="extr",selRed="extr",sdc=2)
{
 #if ="extr" then colour extreme values in a (pre) latex object
 #otherwise colour red/blue according to logical vars selBlue, selRed
 lat.num <- as.numeric(lat); dim(lat.num) <- dim(lat)
 colsds <- rep(apply(lat.num,2,sd,na.rm=T),each=nrow(lat))
 colmns <- rep(apply(lat.num,2,mean,na.rm=T),each=nrow(lat))
 dim(colsds) <- dim(colmns) <- dim(lat)
 Ff <- !is.na(lat.num)
 if (!is.logical(selBlue) ) {
   selBlue <- which(lat.num[Ff]<(colmns[Ff]-(sdc*colsds[Ff]))) }
 if (!is.logical(selRed)) {
   selRed <- which(lat.num[Ff]>=(colmns[Ff]+(sdc*colsds[Ff]))) }
 lat[Ff][selBlue] <- paste("\\textcolor{blue}{",lat[selBlue],"}",sep="")
 lat[Ff][selRed] <-  paste("\\textcolor{red}{",lat[selRed],"}",sep="")
 return(lat)
}


initialise.excl.files <- function(dir,reset.files=c("ChrAb.txt","DLRS.txt","GCWave.txt","Mean.txt","BadPlates.txt","BadReg.txt","CNVQC.txt"))
{
  # delete existing sample exclusion file if present
  dir <- validate.dir.for(dir,c("excl"),warn=F)
  cur.files <- list.files(dir$excl)
  if(length(reset.files)<1) { return(NULL) }
  for (cc in 1:length(reset.files))
  {
    if(reset.files[cc] %in% cur.files)
    {
      warning(" deleting existing file: ",reset.files[cc]," for initialisation purposes")
      unlink(cat.path(dir$excl,reset.files[cc]))
    }
  }
}


gs.bash.http <- function(repos="plumbCNV",script="getDataGS.sh") {
  # load package
  if(packages.loaded("RCurl")) {
    txt <- getURL(paste("https://raw.github.com/nicholasjcooper/",repos,"/master/",basename(script),sep=""), followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))
    return(txt)
  } else {
    warning("package RCurl not available, attempt to run bash script for data import failed")
    return(NULL)
  }
}





init.DATA.read <- function(dir,doLRR=T,doBAF=F,plink.imp=F,n.cores=1,scr.name="getDataGS.sh",scr.dir="",
                           hwe.thr=10^-8, callrate.samp.thr=.95, callrate.snp.thr=.95,
                           snp.info.sup="snpdata.map",genome.stud.file=F,combine.files=F) 
{
  dir <- validate.dir.for(dir,c("scr","baf.col","lrr.dat","base","col","sup","raw"))
  
  if(combine.files) {
    combin <- "CD"
  } else {
    combin <- ""
  }
  if(is.null(scr.dir) | all(scr.dir=="")) {
    scr.dir <- dir$scr
  }
  # make sure file.spec.txt is in place and if not, make a template ready to modify
  fspec <- cat.path(dir$lrr.dat,"file.spec.txt")
  if(!file.exists(fspec)) {
    cat("\n\n<<<<Error: plumbCNV requires setup file 'file.spec.txt' to run\n")
    cat("Will attempt to make a starting template now which can be\n")
    cat("easily modified to a correct setup file>>>>. \n\n")
    cat("\nRunning bash data import script in template mode\n")
    my.cmd <- paste(" -t -F '",dir$raw,"' -O '",dir$base,"'",sep="")
    scr.file <- cat.path(scr.dir,scr.name,must.exist=T)
    cmd <- paste(scr.file,my.cmd,collapse="",sep="")
    cat(cmd,"\n")
    system(command=cmd)
    #ie: ./getDataGS.sh -t -F '/Raw/GenomeStudio' -O '/MyCNV'
    cat("Please modify template in",dir$lrr.dat,"\n")
    stop()
  }
  if (genome.stud.file) { 
    gs <- "G" 
  } else { 
    gs <- "" 
    if(is.file(snp.info.sup,dir$ano,dir)) {
      nc <- nchar(snp.info.sup)
      suf <- substr(snp.info.sup,nc-3,nc)
      if(toupper(suf) %in% c(".BIM",".VCF",".MAP")) {
        cat("\nsnp file is",suf,"type, which should be supported\n")
      } else {
        warning(paste("not expecting snp file with type ",suf,", snp information may not be ",sep=""),
        "imported successfully. \ncheck results of this script carefully.")
      }
      snp.info.sup <- find.file(snp.info.sup,dir$ano,dir)
    } else {
      stop("Error: snp support file",snp.info.sup,"could not be found")
    }
  }
  ## Process LRR import ##
  if(doLRR) {
    if(!plink.imp) {
      #import LRR for QC in R:     
      cat("\nRunning bash data import script\n")
      my.cmd <- paste(" -lSM",gs,combin," -N ",n.cores," -T 'LRR' -F '",dir$raw,"' -O '",
                      dir$base,"' -m '",snp.info.sup,"'",sep="")
      scr.file <- cat.path(scr.dir,scr.name,must.exist=T)
      cmd <- paste(scr.file,my.cmd,collapse="",sep="")
      cat(cmd,"\n")
      system(command=cmd)
    } else {
      #import LRR and do plink QC:  
      cat("\nRunning bash data import script and calling plink for snp QC\n")
      my.cmd <- paste(" -SLMRPlf",gs,combin," -N ",n.cores,
                      " -T 'LRR' -F '",dir$raw,"' -O '",
                      dir$base,"' -m '",snp.info.sup,"'",
                      " -x ",callrate.samp.thr," -y ",callrate.snp.thr," -z ",format(hwe.thr,digits=9,scientific=F),
                      sep="")
      scr.file <- cat.path(scr.dir,scr.name,must.exist=T)
      cmd <- paste(scr.file,my.cmd,collapse="",sep="")
      cat(cmd,"\n")
      system(command=cmd)
    }
  }
  lrr.filez <- c(list.files(dir$col))
  ## Process BAF import ##
  if(doBAF) {
    cat("\nRunning bash data import script for BAF data\n")
    if(combin!="") { combin <- paste(" -",combin,sep="") }
    my.cmd <- paste("",combin," -N ",n.cores," -T 'BAF' -F '",dir$raw,"' -O '",dir$base,"'",sep="")
    scr.file <- cat.path(scr.dir,scr.name,must.exist=T)
    cmd <- paste(scr.file,my.cmd,collapse="",sep="")
    cat(cmd,"\n")
    system(command=cmd)
  }
  baf.filez <- list.files(dir$baf.col)
  cat("imports using bash scripts, complete!\n")
  out <- list(lrr.filez,baf.filez)
  names(out) <- c("lrr.files","baf.files")
  return(out)
}


run.SNP.qc <- function(DT=NULL, dir=NULL, import.plink=F, HD.mode=F, restore.mode=F,
                         callrate.samp.thr=.95, callrate.snp.thr=.95, hwe.thr=10^-7,
                         group.miss=F, grp.hwe.z.thr=4, grp.cr.thr=.001, min.snps=10,
                         het.lo=.1,het.hi=0.4,
                         tabulateCallRate=T, n.cores=1, low.ram=T, drawVennsOfHWEvsCR=T,
                         doDensityPlots=T, sample.info=NULL, snp.info=NULL,
                         sample.fn = NULL, snp.fn="snpNames.txt", snp.fields=NULL, 
                         geno.files = NULL,file.loc="snpMatLst",verbose=T,write.files=T,coverage.loss=T,
                         snp.info.fn="snpdata.map",anot="map3",ucsc="hg18",autosomes.only=T,...)
{
  # load in all the SNP-data need to do SNP-QC to exclude bad snps.
  # use snpStats functions to get callrates, etc, else can use plink results processed
  # at an earlier stage during 'getRawData.R' or 'getDataGS.sh' execution.
  # do SNP-qc and export lists of excluded samples, build sample.info and snp.info objects
  #
  load.all.libs() # load all main plumbcnv libraries
  must.use.package("snpStats",bioC=T) 
  if(n.cores>1) {  must.use.package("multicore",F) }
  # test for 
  if(restore.mode & import.plink) { cat("setting RESTORE=0, not compatible with PLINK=1\n"); restore.mode <- F }
  if(HD.mode & import.plink) { cat("setting HD=0 because HD=1 is not compatible with PLINK=1\n"); HD.mode <- F }
  if(is.data.tracker(DT)) {
    DTs.required <- c("genofile","ngrps","sample.info","snp.info","snps",
                      "subsamples","snp.fields")
    is.valid <- req.valid(DT,DTs.required)
    if(!is.valid) { 
      warning("file(s) required for 'run.SNP.qc' invalid in data.tracker")
      # possibly ressurected by additional input parameters
    } else {
      #sample.info <- getSlot(DT,"sample.info",ret.obj=T)
      #snp.info <- getSlot(DT,"snp.info",ret.obj=T)
      snp.fn <- getSlot(DT,"snps")
      sample.fn <- unlist(getSlot(DT,"subsamples"))
      snp.fields <- getSlot(DT,"snpfields")
      geno.files <- unlist(getSlot(DT,"genofile")) 
      while(is(snp.fields)[1]=="list") { snp.fields <- snp.fields[[1]]; cat("unlisting!\n") }
    }
    if(is.null(dir)) {
      dir <- getSlot(DT,"dir")
    }
  }
  dir <- validate.dir.for(dir,c("ano","cr","raw","col","cr.plk","ids","sort","sort2","excl","excl2"),warn=T)
  file.loc <- cat.path(dir$cr,file.loc,ext=".RData")

  ###DEP##snp.info.fn <- cat.path(dir$ano,snp.info.fn,must.exist=T)
  # read in snp list 
  snp.fn <- find.file(snp.fn,dir$ano,dir)
  if (file.exists(snp.fn)) {
    snpIDs <- readLines(cat.path(dir$ano,snp.fn))
  } else {
    linez <- character()
    linez[1] <- paste("Error: expecting file: '",snp.fn,"' in ",dir$ano,".",sep="")
    linez[2] <- "This file speeds up import by giving the snp order in the raw file"
    linez[3] <- "File should be generated using the bash command:"
    linez[4] <- paste("grep <anySampleid> <mydatafiles> | cut -f 2 | cat >",snp.fn)
    cat(linez,"\n")
    stop()
  }    
  # load and prepare a SNP info object if not using DT or not passed in as obj
  if(is.null(snp.info)) {
    snp.info <- make.snp.info(dir, snpIDs=snp.fn, anot=anot, ucsc=ucsc, 
                            snp.info.fn=snp.info.fn,sav=F,verbose=F,filt.by=NULL,
                            non.autosomes.excl.file=T)
  }
  if(is.null(sample.info)) { sample.info <- read.sample.info(dir) }
  if(n.cores>length(sample.fn)) { use.multi <- T } else { use.multi <- F }
  if(low.ram) { l.cores <- 1 } else { l.cores <- n.cores } # number of cores to use when ram limited
  if (!import.plink & !restore.mode) 
  { 
    # RESTORE = 0, PLINK = 0 #
    cat("\nImporting SNP data into R using snpStats package\n")
    snpMatLst <- import.snp.matrix.list(snpIDs,dir=dir,data=geno.files,samples=sample.fn,n.cores=n.cores,
                                        snp.fields=snp.fields,HD=HD.mode,multi=use.multi, verbose=verbose) 
    #print(paste("Length of snpMatLst:",length(snpMatLst))); print("temporary backup files to /backup")
    #system("cp /chiswick/data/ncooper/metabochipRunTest/CALLRATES/snpMat*RData /chiswick/data/ncooper/metabochipRunTest/CALLRATES/backup")
    if(!HD.mode) { 
      save(snpMatLst,file=file.loc); cat("~wrote snpMatLst to:",file.loc,"\n") 
      if(is.data.tracker(DT)) { DT <- setSlot(DT=DT,snp.stats=file.loc) } # list of all grp objs in one file
    } else {
      if(is.data.tracker(DT) & is.ch(snpMatLst)) { DT <- setSlot(DT=DT,snp.stats=snpMatLst) } # list of grp file names
    }
    snpMatLst <- sync.snpmat.with.info(snpMatLst,snp.info=snp.info,sample.info=NULL,dir=dir,n.cores=l.cores) 
  } else { 
    if(restore.mode)
    {
      if(!HD.mode & file.exists(file.loc)) {
        # RESTORE = 1, PLINK = 0, HD = 0, file.loc exists #
        cat(paste("\nRestoring SnpMatrix list from",file.loc,"\n"))
        snpMatLst <- get(paste(load(file.loc)))
        snpMatLst <- sync.snpmat.with.info(snpMatLst,snp.info=snp.info,sample.info=NULL,n.cores=l.cores)
      } else {
        # RESTORE = 1, PLINK = 0, HD = 1 or file.loc doesn't exist #
        #print(snpMatLst)
        if(HD.mode & file.exists(paste(dir$cr,"snpMat",1,".RData",sep=""))) {
          cat(paste("\nRestoring multiple SnpMatrix files from",paste("snpMat*.RData\n")))
          existo <- T; cnt <- 1; snpMatLst <- list()
          while(existo){
            nfn <- paste(dir$cr,"snpMat",cnt,".RData",sep="")
            if (file.exists(nfn)) 
            { snpMatLst[[cnt]] <- nfn ; cnt <- cnt+1 
            } else {
              existo <- F
            }
          }
          snpMatLst <- sync.snpmat.with.info(snpMatLst,snp.info=snp.info,sample.info=NULL,dir=dir,n.cores=l.cores)
        } else {
          # RESTORE = 1, PLINK = 0, no snpMat*RData files exist #
          stop("Restore mode failed, expected .RData files not present\n")
          snpMatLst <- NULL
        }
      }
    } else {
      # RESTORE = 0, PLINK = 1 #
      cat("\nUsing Plink QC files as input\n")
      snpMatLst <- NULL 
    }
  }
  # read in subject list
  if(is.null(sample.fn)) {
    subIDs.actual <- get.subIDs(dir,"combined")
  } else {
    if (is.file(sample.fn,dir$ids,dir)) {
      subs.list <- sapply(find.file(sample.fn,dir$ids,dir),readLines)
      subIDs.actual <- unlist(subs.list)
    } else {
      stop("Error: could not find sample IDs")
    }
  }
  # move inputted sample info 
  orig.sample.info <- sample.info; rm(sample.info)
  samp.qc.list <- doSampQC(dir=dir, subIDs.actual=subIDs.actual, plink=import.plink,
                           callrate.samp.thr=callrate.samp.thr, snpMatLst=snpMatLst, 
                           het.lo=het.lo,het.hi=het.hi,
                           n.cores=l.cores, proc=2)
  sample.info <- samp.qc.list$SAMPLE.INFO
  # Tabulate sample call rate stats
  if(tabulateCallRate) {  samp.result <- call.rate.summary(sample.info) }
  # write sample call rate failure lists to text file #
  ofn <- cat.path(dir$cr,"sampleqc.txt")
  write.table(sample.info,file=ofn,sep="\t",quote=F)
  cat(paste("~wrote file with samples call rate QC summary to:\n ",ofn,"\n"))
  ofn <- cat.path(dir$excl,"callrate.txt")
  writeLines(samp.qc.list$CR.EXCL,con=ofn)
  cat(paste("~wrote file with samples to exclude because of low call rate:\n ",ofn,"\n"))
  ofn <- cat.path(dir$excl,"hz.txt")
  writeLines(samp.qc.list$HZ.EXCL,con=ofn)
  cat(paste("~wrote file with samples to exclude because of heterozygosity outside (",het.lo,",",het.hi,"):\n ",ofn,"\n",sep=""))
  
  if(!import.plink) { snpMatLst <- sync.snpmat.with.info(snpMatLst,snp.info=NULL,
                                    sample.info[sample.info$QCfail==0,],dir=dir,n.cores=l.cores) } 
  
  orig.snp.info <- snp.info
  snp.qc.list <- doSnpQC(dir, plink=import.plink, n.cores=l.cores, proc=2,
                         callrate.snp.thr=callrate.snp.thr, hwe.thr=hwe.thr,
                         grp.hwe.z.thr=grp.hwe.z.thr, grp.cr.thr=grp.cr.thr,
                         snpMatLst=snpMatLst, subIDs.actual=subIDs.actual, group.miss=group.miss,
                         snp.info=snp.info, sample.info=sample.info, autosomes.only=autosomes.only)
  # still contain same core ranges/ids as original but if autosomes.only=T, then these now removed
  snp.info <- snp.qc.list$SNP.INFO
  
  # Tabulate snp call rate stats
  if(tabulateCallRate) { snp.result <- call.rate.summary(snp.info) }
  ofn <- cat.path(dir$cr,"snpqc.txt")
  rsnp.info <- snp.info
  # round before writing to file
  for (cc in 1:ncol(rsnp.info)) { if(is.numeric(rsnp.info[[cc]])) { (rsnp.info[[cc]] <- round(rsnp.info[[cc]],digits=5)) } }
  print(colnames(rsnp.info))
  sel.cols <- which(colnames(rsnp.info) %in% c("call.rate","P.hwe","Z.hwe",
                                   "grp.hwe.zmax","grp.miss.p","QCfail","het"))
  write.table(as.data.frame(rsnp.info[,sel.cols]),row.names=F,file=ofn,sep="\t",quote=F)
  cat(paste("~wrote file with snp call rate and HWE QC summary to:\n ",ofn,"\n"))
  ofn=cat.path(dir$excl2,"callrate.txt")
  writeLines(snp.qc.list$CR.EXCL,con=ofn)
  cat(paste("~wrote file with call rate failing snps to:\n ",ofn,"\n"))
  writeLines(snp.qc.list$HWE.EXCL,con=cat.path(dir$excl2,"HWE.txt"))
  cat("~wrote HWE exclusion list to annotion directory\n")
  if(group.miss & all(c("GRP.CR.EXCL","GRP.HWE.EXCL") %in% names(snp.qc.list))) {
    writeLines(snp.qc.list$GRP.CR.EXCL,con=cat.path(dir$excl2,"grp.missingness.txt"))
    cat("~wrote grp-wise missingness discrepancy list to annotion directory\n")
    writeLines(snp.qc.list$GRP.HWE.EXCL,con=cat.path(dir$excl2,"grp.HWE.txt"))
    cat("~wrote grp-wise HWE discrepancy list to annotion directory\n")
  }

  ## REMOVE MISSING VALUES, STORING ORIGINALS ##
  # note that graphing and tables below work poorly with missing
  if(nrow(snp.info[!is.na(snp.info[,1]),])>(length(snpIDs)/2))
  {
    # so long as there are a reasonable number of snps here, replace sort list 
    # in annotation directory with the sorted list just derived
    # (prevents replacement in most cases of crash/incorrect input)
    ex.fls <- list.files(dir$sort2)
    if(length(ex.fls)>0) { delete.file.list(cat.path(dir$sort2,ex.fls)) }
    cat("replacing any files in directory /ANNOTATION/SNP_SORT/ with new sorted snp list\n")
    ofn <- cat.path(dir$sort2,"snpsort.txt")
    writeLines(rownames(snp.info),con=ofn)
    cat("~wrote file:",ofn,"\n")
  } else {
    warning("number of valid snps in 'snp.info' object is less than half of original list, please check files")
  }
  
  # View overlap of call rate and HWE exclusion criteria for SNPs
  if(drawVennsOfHWEvsCR) {
    must.use.package("limma",T)
    cat(" generating venn diagrams\n")
    conds <- cbind(snp.qc.list$CR.cond,snp.qc.list$HWE.cond)
    conds[is.na(conds)] <- FALSE; vc <- vennCounts(conds) ; vc
    ofn <- paste(dir$cr,"venn_",(100*callrate.snp.thr),".pdf",sep="")
    pdf(ofn)
    vennDiagram(vc,names=c(paste("Call rate<.",callrate.snp.thr,sep=""),paste("HWE p<",hwe.thr,sep="")))
    dev.off()
    cat(paste("~wrote Venn figure to file:",ofn,"\n"))
  }
  
  # Visualize call rate stats
  if(doDensityPlots) {
    fn <- paste(dir$cr,"callRateDistributions.pdf",sep="")
    draw.density.plots(fn,sample.info,snp.info,samp.result=NULL,snp.result=NULL,
                       callrate.samp.thr=callrate.samp.thr, 
                       callrate.snp.thr=callrate.snp.thr)
    hwe.density.plots(Z.hwe=snp.info$Z.hwe,hwe.thr=hwe.thr,dir=dir)
    het.density.plots(data="sampleqc.txt",het.lo=.1,het.hi=0.4,zoom=T,
                                  het=NULL,dir=dir,fn="HZDistribution.pdf") 
    hwe.vs.callrate.plots(call.rate=snp.info$call.rate,Z.hwe=snp.info$Z.hwe,
                          callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,dir=dir)
    hz.vs.callrate.plots(data="sampleqc.txt",callrate.samp.thr=callrate.samp.thr,
                                    het.lo=.1,het.hi=0.4,
                                     zoom=T,full=T,dir=NULL,
                                     fn="HZvsCallrate.pdf")
  }
  
  if(coverage.loss) {
    # report on loss of genome coverage due to QC excluded snps #
    coverage.change(snp.info,cnv.sizes=c(50,100,200,500),dir=dir,min.snps=min.snps,do.plot=T) 
  }
  
  ## CHECK SNP.INFO AND SAMPLE.INFO AGAINST ORIGINALS ##
  if(!is.null(orig.snp.info)) {
    if(nrow(orig.snp.info)> nrow(snp.info)) {
      ##### SHOULD WE DO THIS??? - i think so, never hurts to have too many? #######
      ## if snp.info has reduced in length because of the snpQC removing autosomes,
      # just add the new columns to the original object instead of replacing it
      to.add <- c("call.rate","P.hwe","QCfail")
      #print(dim(orig.snp.info)) ; print(dim(snp.info)); print(head(orig.snp.info));print(head(snp.info))
      for (cc in 1:length(to.add)) {
      	selector <- (rownames(orig.snp.info) %in% rownames(snp.info))
        orig.snp.info[[to.add[cc]]] <- NA
        orig.snp.info[[to.add[cc]]][selector] <- snp.info[[to.add[cc]]][match(rownames(orig.snp.info)[selector],rownames(snp.info))]
      }
      snp.info <- orig.snp.info 
    } 
  }
  if(!is.null(orig.sample.info)) {
    ## if sample.info has reduced in length because of the QC somehow,
    # just add the new columns to the original object instead of replacing it
    new.sample.info <- sample.info
    sample.info <- add.to.sample.info(orig.sample.info,new.sample.info,verbose=verbose)
  }
  if(write.files | is.data.tracker(DT)) {
    f1 <- write.sample.info(sample.info,dir,T)
    f2 <- write.snp.info(snp.info,dir,verbose=T)
    cat("~wrote sampleinfo.tab, snpinfo.tab to:\n ",dir$ano,"\n")
  }
  n.snp.fail <- length(which(snp.info$QCfail!=0)); n.snp.pass <- nrow(snp.info)-n.snp.fail
  n.samp.fail <- length(which(sample.info$QCfail!=0)); n.smp.pass <- nrow(sample.info)-n.samp.fail
  cat("SNP-QC removed a total of",n.samp.fail,"samples and",n.snp.fail,"SNPs\n")
  if(is.data.tracker(DT)) {
    suppressWarnings(DT <- setSlot(DT,sample.info=f1, snp.info=f2,proc.done=2,warns=F))
    return(DT)
  }
  outlist <- list(sample.info,snp.info)
  names(outlist) <- c("sample.info","snp.info")
  return(outlist)
}


do.chromosomal.abberations <- function(dir,bigMat2,snp.info,pref="",
                                       colPlot=F,lob=2,hib=2.5,pctile.bound=0.01,plot=T,skip.chr.ab=F,...) {
  # detect, flag and plot samples with chromosomal abberations
  # better off passing a complete snp.info object, as will filter out any not in the data automatically
  # if skip.chr.ab = T then just write a blank exclusion file and continue
  if(skip.chr.ab) {
    ofn <- cat.path(dir$excl,"ChrAb.txt")
    writeLines("",con=ofn) 
    cat("\n~wrote empty chr.ab excluded subject list to:",ofn,"\n")
    return(NULL)
  }
  ## otherwise do chromosomal abberation detection
  if(!is(snp.info)[1]=="RangedData") { stop("snp.info was not a valid RangedData object") }
  snp.info <- snp.info[(rownames(snp.info) %in% rownames(bigMat2)),]  # snp.info object matching the bigmatrix rows exactly
  cat("\nDetecting chromosomal abberations\n")
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
  if(n.chr<2) { 
    cat(" only 1 chromosome in dataset so skipping detection of abberations\n\n")  
    out.list <- list(badcheckz=character(0),chrN=NULL,chrLab=character(0))
    return(out.list)
  }
  chr.stat <- get.chr.stats(bigMat2,snp.info,dir)
  chrWarns <- do.chr.ab.contrasts(chr.stat$chr.mean,lob,hib,nC=n.chr)
  outlist <- get.chr.ab.fail.subset(chrWarns, dir=dir, failerMode="NOTLRR") 
  chr.ab.report(chr.stat, chrWarns=chrWarns, dir=dir, 
                makeGraphs=F, writeExclList=T,append=T)
  if(plot) {
    plot.chr.ab.samples(dir=dir,bigMat2=bigMat2,chr.stat=chr.stat, chr.ab.samples.obj=outlist,
                        snp.info=snp.info,colPlot=colPlot,failerMode="NOTLRR",
                        lob=lob,hib=hib,pctile.bound=pctile.bound,pref=pref,...)
  }
  return(outlist)
}


do.plate.qc <- function(dir,IDlist=NULL,sample.info=NULL,stat.table=NULL,stat.sum=NULL,
                                      batch="plate",pref="",plot.and.report=T,badPlateThresh=.33,
                                      plate.fn="plate.lookup.txt") {
  ### PLATE STATS ###
  cat("\nPerforming plate QC\n")
  # ID.list is a subset (e.g, current bigMat subjs) to do these stats and plots for
  if(!is.character(batch)) { batch <- "plate" } # ensure valid value for important parameter
  if(length(batch)>1) { 
    warning("batch QC only possible for one category at a time")
    cat("Selecting:",batch[1],"and ignoring",paste(batch[-1],collapse=","),"\n")
    batch <- batch[1]
  }
  if(!is.null(sample.info)) {  
    other.cols <- batch[!tolower(batch) %in% c("well","plate")]
    if(length(other.cols)==0) { other.cols <- NULL }
    plate.lookup <- get.plate.from.sample.info(sample.info,other.cols=other.cols) 
  } else { 
    if(length(grep("plate",tolower(batch)))<1) { 
      warning(paste("'batch' was",batch[1],"but as sample.info was NULL, defaulting to 'plate'"))
      batch <- "plate"
    }
    plate.lookup <- get.plate.info(dir=dir,fn=plate.fn)[[1]]
    if(!is.null(plate.lookup)) { stop("Error: couldn't find plate.lookup.txt file\n")}
  }
  new.fr <- update.plate.bad.count.table(dir,plate.lookup,filt.list=IDlist,type=batch, incl.cr=T)
  ## yields here total count set for this plate only (because of filt.list above):
  this.fr <- new.fr[!is.na(new.fr$SIZE),]
  cat("\nTable of",toheader(batch),"failure counts:\n"); print(this.fr)
  ofn <- paste(dir$qc.pl,"CountsPer",toheader(batch[1]),pref,".tab",sep="")
  write.table(this.fr,sep="\t",col.names=T,row.names=F,quote=F,file=ofn)
  cat("written to",ofn,"\n")
  # the update.plate.... function....may also be useful when combining all file results
  
  ## PLATEWISE EXCLUSION - writes excluded samples to annotation directory
  # removes plates with more than one third samples excluded by LRR sample QC
  excl.info <- excl.bad.plates(this.fr,plate.lookup,dir=dir,badPlateThresh=badPlateThresh,
                               writeExclList=T,append=T,batch=batch, incl.cr=T)
  
  ### PLOTS / REPORTS ###
  if(plot.and.report & !is.null(stat.table) & !is.null(IDlist)) {
    ## function to do those scatter plots
    plate.lrr.stats <- round(get.plate.lrr.stats(plate.lookup,stat.table,IDlist),3)
    ofn <- paste(dir$qc.pl,"StatsPer",toheader(batch),pref,".tab",sep="")
    write.table(plate.lrr.stats,sep="\t",col.names=T,row.names=T,quote=F,file=ofn)
    cat("~table of LRR-statistics by ",toheader(batch),", written to:\n ",ofn,"\n",sep="")
    if(!is.null(stat.sum)) {
      batch.box.plots(stat.table,stat.sum,plate.lookup,batch=batch,pref=pref,dir=dir)
    }
  }
  return(excl.info)
}  


process.cohort.qc <- function(DT=NULL,grp,of,dir,sample.info,snp.info,restore.mode=F,
                              pref="",des.fn=paste(pref,"descrFile",sep=""),verbose=F,
                              med.chunk.out="medianStorage.RData",plate.fn="plate.lookup.txt",
                              dlrs=T,badPlateThresh=0.33,skip.chr.ab=F,lob=2,hib=2.5,pctile.bound=0.01,batch="plate",
                              plate.report=T,chr.ab.report=T,lrr.report=T,n.cores=1,
                              cohort.pc.correct=F, pc.to.keep=.13,num.pcs=9,nSD=3,
                              lo.cut=c("LB",NA,"LB"),hi.cut=c("UB","UB","UB"),ucsc="hg18",...)
{
  ### loop through once for each file
  ### passes arguments to lrr.sample.dists which passes to calc.chr.info
  # better off passing a complete snp.info object, as subfunctions will filter out any not in the data automatically
  dir <- validate.dir.for(dir,c("big","col","ano","qc.cab","qc.pl","qc.lrr","gc"),warn=T)
  load.all.libs() # load all plumb cnv libraries 
  #auto.mode <- auto.mode()
  # chromosomal aberattion z score boundary
  #hib <- 2.5 #2.575 #3.29 # define stringent z score boundary 
  #lob <- 2 #1.96 #2.575 # define lenient z score boundary
  # lo.cut=c("LB",NA,"LB"),hi.cut=c("UB","UB","UB") are which cutoff to use for the columns of the stat.table
  # which for this function should be: Mean, DLRS and GCWave respectively
  if(of==1) {
    pref <- "LRR" #  sample.fn <- "subIdsALL.txt"
  } else {
    pref <- paste("LRR",grp,sep="")
  }
  cat("\nProcessing data import and sample QC for cohort",grp,"of",of,"\n\n")
  ###des.fn <- paste(pref,"descrFile",sep="")
  med.chunk.out <- cat.path(dir$gc,pref=pref,fn=med.chunk.out)
  
  ## LOAD RAW DATA FROM SOURCE##
  if(file.exists(cat.path(dir$big,des.fn))) {
    cat("Loading big matrix from descriptor file",des.fn,":\n")
  } else { stop(paste("Error: big.matrix file",des.fn,"not found")) }
  
  ## FILTER ON SNP-QC, call rate analysis prepared earlier, etc ##
  #TEMPORARY FOR DEBUGGING!!! DEP
  #debug(big.exclude.sort)
  if(F) { R.descr <- cat.path(dir$big,"LRRSortdescrFile.RData") } else {
  R.descr <- big.exclude.sort(des.fn, dir=dir, pref=pref,verbose=verbose) }
  bigMat2 <- getBigMat(R.descr,dir)
  print.big.matrix(bigMat2,"Snp_QC_Mat")
  
  datlist <- lrr.sample.dists(bigMat2,snp.info=snp.info,dir=dir,n.cores=n.cores,
                              gc=T,dlrs=dlrs,pref=pref,med.chunk.fn=med.chunk.out,
                              ucsc=ucsc,restore.mode=restore.mode,
                              nSD=nSD,lo.cut=lo.cut,hi.cut=hi.cut,write.excl=T,
                              extended.report=lrr.report,...)

  ## do chromosomal abberations
  outlist <- do.chromosomal.abberations(dir=dir,bigMat2=bigMat2,snp.info=snp.info,skip.chr.ab=skip.chr.ab,pref=pref,
                                         colPlot=F,lob=lob,hib=hib,pctile.bound=pctile.bound,plot=chr.ab.report)

  ## do plate QC and plots and reports
  plate.list <- do.plate.qc(dir,IDlist=colnames(bigMat2),sample.info=sample.info,stat.table=datlist$stat.table,
              stat.sum=datlist$stat.summary,batch=batch,pref=pref,plot.and.report=plate.report,
              badPlateThresh=badPlateThresh,plate.fn=plate.fn)
  
  ## R.descr <- big.exclude.sort(descr, dir=dir, pref=paste("SampQC_LRR",grp,sep=""))
  ## want to return med.chunk.out somehow # return R.descr , datlist$tab.file  ##to pass to DT
  ## or could update DT here....
  ### would do extra PCA here if desired....
  if(cohort.pc.correct & of>1) {
    # do PCA and correction just for this cohort...
    cat("\nRunning extra cohort-wise PCA and correction on cohort",grp,"of",of,"\n")
    evs.filename <- cat.path(dir$pc,pref="Cohort",fn=grp,suf="PCsEVsFromPCA",ext="RData")
    filt.descr <- big.exclude.sort(R.descr, dir=dir, 
                                   pref=paste("SampQC_LRR",grp,sep=""),verbose=F)
    pc.descr <- run.PCA.correct(dir=dir,pc.to.keep=pc.to.keep,assoc=F,correct.sex=F,
                                num.pcs=num.pcs,n.store=max(50,num.pcs),add.int=F,
                                comparison=F,big.lrr.fn=filt.descr,n.cores=n.cores,
                                big.cor.fn=paste("pcaSubCorrected",grp,".RData",sep=""),
                                snp.sub.fn="pca.snp.subset.pre.txt",restore.mode=(grp!=1),pref.suf=paste(grp),
                                sub.pref="cohortPCA",cor.pref="cohortCorrect",
                                sub.fn=paste("pcaSubMat",grp,".RData",sep=""),
                                pcs.fn=evs.filename, med.fn="postpcmed.RData")
    if(is.data.tracker(DT)) { 
      DT <- setSlot(DT,pcas=pc.descr,eigens=evs.filename,grps=grp,warns=F) 
      cleanup <- T; if(cleanup) { unlink(filt.descr) }
    } else {
      warning("returning pc-corrected cohort subset, but data.tracker not in use.")
      return(pc.descr)
    } 
  }
  if(is.data.tracker(DT)) {
    gc.human.fn <- cat.path(dir$gc,"GCAv6.RData") # it always gets this name, see get.gc.human()
   ## print("big.filt empty (i hope)"); print(getSlot(DT,"big.filt",grps=grp)); print(R.descr)
    DT <- setSlot(DT,big.filt=R.descr,medians=med.chunk.out,stats=datlist$tab.file,
                        gc.genome=gc.human.fn,grps=grp,warns=T)
   ### print("big.filt filled (i hope)"); print(getSlot(DT,"big.filt",grps=grp))
    return(DT)
  } else {
    return(R.descr)  
  }
}


run.SAMPLE.qc <- function(DT=NULL,dir=NULL,init=T,big.lrr=NULL,sample.info=NULL,snp.info=NULL,ucsc="hg18",restore.mode=F,
                          med.chunk.out="medianStorage.RData",big.list.out="grpsBigMats.RData",verbose=F,
                          cmb.out="combinedBigMat.RData",badPlateThresh=0.33,skip.chr.ab=F,lob=2,hib=2.5,pctile.bound=0.01,
                          pc.to.keep=.13,num.pcs=9,nSD=3,mean.thr=c("LB","UB"),dlrs.thr=c(NA,"UB"),gc.thr=c("LB","UB"),
                          dlrs=T,batch="plate",lrr.report=T,chr.ab.report=T,plate.report=T,cohort.pc.correct=F,n.cores=1,...)
{
  # run the whole sample QC process from scratch including importing raw data
  load.all.libs() # load all main plumbcnv libraries
  # check validity of thresholds before proceeding and warn so changes can be made before time wasted
  valid.pars <- thresh.check.valid(nSD=nSD,mean.thr=mean.thr,dlrs.thr=dlrs.thr,gc.thr=gc.thr)
  if(!valid.pars[[3]]) { cat("*WARNING: detected invalid parameters for SAMPLE-QC thresholds.",
                   "Recommend fixing 'mean.thr/dlrs.thr/gc.thr' and restarting, else these will be set to NA\n") }
  lo.cut=valid.pars[[1]]; hi.cut=valid.pars[[2]]
  ## establish DT mode (or not)
  if(is.data.tracker(DT)) {
    DTs.required <- c("ngrps","sample.info","snp.info","big.lrr","subsamples")
    is.valid <- req.valid(DT,DTs.required)
    if(!is.valid) { 
      warning("file required for 'run.SAMPLE.qc' invalid in data.tracker")
    } else {
      # don't need these three now, but will need them in subroutines
      #sample.info <- getSlot(DT,"sample.info",ret.obj=T)
      #snp.info <- getSlot(DT,"snp.info",ret.obj=T)
      #sample.fn <- unlist(getSlot(DT,"subsamples"))
      big.lrr <- getSlot(DT,"big.lrr")
      ngrp <- getSlot(DT,"ngrp")
    }
    if(is.null(dir)) {
      dir <- getSlot(DT,"dir")
    }
  } else {
    # if filenames not inputted assume 1 group, else number of file names
    if(!is.null(big.lrr)) { ngrp <- length(big.lrr) } else { cat(" no file names inputted, assuming 1 group\n") ; ngrp <- 1 }
  }
  dir <- validate.dir.for(dir,c("big","res","excl"))
  # delete any existing sample exclusion files
  if(init) { initialise.excl.files(dir); cat(" initializing sample-QC exclusion files\n") }
  ## gather sample information (firstly loading from LoadInSNPs.R result)  
  if(is(snp.info)[1]!="RangedData") { snp.info <- read.snp.info(dir) }
  # if no snpsort file, assume snp.qc not run, so should create an ordered one
  sst.fn <- cat.path(dir$sort2,"snpsort.txt")
  if(!file.exists(sst.fn)) {
    writeLines(rownames(toGenomeOrder(snp.info)),con=sst.fn)
    cat(" created SNP sorting file based on ordered snp.info\n")
  }
  if(!is.data.frame(sample.info)) { sample.info <- read.sample.info(dir) }
  if(is.data.tracker(DT)) { sample.info <- make.grp.column(DT,sample.info,add=T) }
  sample.info <- validate.samp.info(sample.info,dir=dir,QC.update=T,verbose=F,proc=1)
  ## PERHAPS THIS ALL SHOULD BE FORCED EARLIER IN MAKE.SAMPLE.INFO - then get plate from sampleinfo in process.coho..
  ## make sure sample info is nice, in particular that 'grp' matches 'GRP' in 'file.spec.txt'
 # sample.info <- validate.samp.info(sample.info,dir) #DEP file.spec.txt
  # load plate info  #### ADD THIS FILE LOCATION TO DT OBJECT and to args!!!!! ####
  ####plate.lookup <- get.plate.info(dir,1,2,3,dup.action="trim")  # 1,4,5
  # add plate annotation to sample.info
  ##sample.info <- add.to.sample.info(sample.info,plate.lookup[[1]],to.add=c("plate","well"))
  ##file.info <- get.file.specs(dir)
  descr.list <- vector("list", ngrp)
  cat("\nRunning sample QC for",ngrp,"cohort(s)\n")
  for (grp in 1:ngrp) {
    cohort.qc <- process.cohort.qc(DT=DT,grp=grp,of=ngrp,dir=dir,sample.info=sample.info,snp.info=snp.info,verbose=verbose,
                                   pref="",des.fn=big.lrr[[grp]],med.chunk.out=med.chunk.out,plate.lookup=plate.lookup,
                                   badPlateThresh=badPlateThresh,skip.chr.ab=skip.chr.ab,lob=lob,hib=hib,pctile.bound=pctile.bound,
                                   batch=batch,plate.report=plate.report,chr.ab.report=chr.ab.report,lrr.report=lrr.report,
                                   dlrs=dlrs, pc.to.keep=pc.to.keep,num.pcs=num.pcs,restore.mode=restore.mode,
                                   cohort.pc.correct=cohort.pc.correct,nSD=nSD,lo.cut=lo.cut,hi.cut=hi.cut,n.cores=n.cores,...)
    if(is.data.tracker(DT)) {
      DT <- cohort.qc
      if(cohort.pc.correct) {
        descr.list[[grp]] <- getSlot(DT,"pcas",grps=grp)
      } else {
        descr.list[[grp]] <- getSlot(DT,"big.filt",grps=grp)
      }
    } else {
      descr.list[[grp]] <- cohort.qc
    }
  }
  descr.fn <- cat.path(dir$big,big.list.out)
  save(descr.list,file=descr.fn)
  new.descr <- exclude.combine.big.mats(descr.list,dir) 
  ofn <- cat.path(dir$big,cmb.out)
  save(new.descr,file=ofn)
  summary.list <- make.QC.summary.table(sample.list=rownames(sample.info),dir)
  passQCsamples <- summary.list$pass.samples
 ## sample.info$QCfail[rownames(sample.info) %in% get.all.samp.fails(dir)] <- 1
  sample.info <- validate.samp.info(sample.info,dir=dir,QC.update=T,verbose=F,proc=3) # updates fail flags
  ## update callrates to reflect current set of SNPs ###
  crs <- samp.cr.summary(getBigMat(new.descr,dir$big),print=T)
  sample.info[names(crs),"call.rate"] <- crs
  write.sample.info(sample.info,dir=dir,verbose=T,type="tab")
  if(is.data.tracker(DT)) {
    DT <- suppressWarnings(setSlot(DT,big.qc=ofn,proc.done=3,warns=F))
    return(DT)
  } else {
    return(new.descr)
  }
}


make.grp.column <- function(DT,sample.info,add=F) { 
    if(!is.data.tracker(DT)) { warning("DT not a data.tracker object, grp not returned") ; return(NULL) }
    if(!is.data.frame(sample.info)) { warning("sample.info not a data.frame object, grp not returned") ; return(NULL) }
    DTs.required <- c("ngrps","sample.info","subsamples")
    is.valid <- req.valid(DT,DTs.required)
    if(!is.valid) { warning("data.tracker object incomplete, no grps added to sample.info") }
    n.grp <- getSlot(DT,"ngrps")
    grp <- rep(1,times=nrow(sample.info)); acct.for <- 0
    for (cc in 1:n.grp) {
    	ID.list <- unlist(getSlot(DT,"sub.samples",ret.obj=T,grps=cc),recursive=F)
    	selection <- rownames(sample.info) %in% unlist(ID.list)
      grp[selection] <- cc
      acct.for <- acct.for+length(which(selection))
    }
    if(add) {
      cat(" ",acct.for,.Platform$file.sep,nrow(sample.info),
          " samples successfully assigned groups\n",sep="")
      sample.info[["grp"]] <- grp
      return(sample.info)
    } else {
      return(grp)
    }
}




dat.to.big.matrix <- function(dir, grp=NA, snp.fn="snpNames.txt", sample.fn=NULL, 
                              input.fn=NULL, pref="LRR", delete.existing=T,
                              ret.obj=F,DT=NULL,verbose=T)
{
  # import from a text (hopefully long format) datafile to a big.matrix
  # return bigmatrix description 'object' or name of description file according to 'ret.obj'
  ### CALIBRATE OPTIONS AND PERFORM CHECKS ###
  max.mem <- 4096 # stupidly large
  bafmode <- F # assume in LRR mode unless otherwise indicated by 'pref'
  if(length(grep("BAF",toupper(pref)))>0) {
    ## note these must match initial bash script that extracts data!
    dat.file.suf <- "BAF.dat"; bafmode <- T
  } else {
    dat.file.suf <- "LRR.dat"
  }
  ## Define data types for big.matrix
  dat.type <- double(1)
  dat.typeL <- "double"  # label
  dir <- validate.dir.for(dir,c("ano","big","col"),warn=F)
  #### GET THE SHAPED INPUT FILE NAME(S) ####
  if(all(is.null(input.fn))) {
    ## try to automatically detect datafiles if not specified
    if(is.data.tracker(DT)) {
      if(!bafmode) { input.fn <- getSlot(DT,"shape.lrr",grps=grp) 
      } else { input.fn <- getSlot(DT,"shape.baf",grps=grp) }
      if(is.list(input.fn)) { input.fn <- unlist(input.fn,recursive=F) } #cat("CHECK: removed 1 level of list\n") 
    } else {
      if(!bafmode) { input.fn <- get.lrr.files(dir,grp,suffix=dat.file.suf)    
      } else { input.fn <- get.lrr.files(dir,grp,suffix=dat.file.suf,baf=T) }
    }
  }
  # print(input.fn)
  #### GET THE CORRESPONDING SAMPLE LIST(S) ####
  if(all(!is.null(sample.fn))) {
    sample.fn <- find.file(sample.fn,dir$ids)  # sample id file name
    cat("Reading sample and snp lists from text files...\n")
    cat(paste(" samples entered manually, reading samples from",sample.fn,"\n"))
    ID.list <- lapply(sample.fn,readLines)  #readLines(sample.fn)
  } else {
    if(is.data.tracker(DT)) {
      ## GET samples using data.tracker object 'DT'
      if(all(!is.na(grp))) {
        if(verbose) { cat(" getting ID list(s) for group",grp,"\n") }
        # return file lists, listed by 'group'; separate lists for each file despite grp
        ID.list <- unlist(getSlot(DT,"sub.samples",grps=grp,ret.obj=T),recursive=F)
      } else {
        ng <- getSlot(DT,"ngrps")
        if(length(input.fn)==1 & ng==1) {
          # all samples in 1 set
          if(verbose) { cat(" 1 input file and 1 sample set only, getting IDs\n") }
          ID.list <- list(unlist(getSlot(DT,"samples",ret.obj=T))) 
        } else {
          # separate list for each file listed by group sets
          if(verbose) { cat(" getting all ID lists ignoring group [",ng," groups]\n",sep="") }
          ID.list <- unlist(getSlot(DT,"sub.samples",ret.obj=T),recursive=F) 
        }
      }
    } else {
      ## GET samples using 'file.spec.txt' # deprecated
      stop("if we need this code back, go to before march 19\n")
    }
  }
  ##print(headl(ID.list))
  cmb.ID.list <- paste(do.call("c",ID.list))
  ##print(length(ID.list[[1]]))
  numfls <- length(ID.list)
  ### GET THE ORDERED SNP LIST ###
  if(is.data.tracker(DT)) {
    snp.fn <- getSlot(DT,"snps")
  } else {
    snp.fn <- find.file(snp.fn,dir$ano,dir)  # snp annotation file name
  }
  #print(head(cmb.ID.list))
  cat(paste(" reading snps from",snp.fn,"\n"))
  snp.list <- readLines(snp.fn) ##; print(head(snp.list))
  ## Multiple possible input situations: 
  ##   FOR LRR: unlist all input files and subsamples for selected grp
  ##    n groups in study, 1 file each: should have input.fn length 1, id.list length 1
  ##    n groups in study, m_n files each: should have input.fn length m_n, id.list length m_n
  ##   FOR BAF: unlist all input files and subsamples for all grps
  ##    n groups in study, 1 file each: should have input.fn length n, id.list length n
  ##    n groups in study, m_n files each: should have input.fn length sum_i(m_n=i), id.list length sum_i(m_n=i)
  ## DEBUG print(length(input.fn)); print(numfls); print(input.fn)
  if(length(input.fn)>1) {
    if(length(input.fn)==numfls) {
      if(verbose) {
        warning(paste("reading a single cohort from",numfls,"source files. Edit file.spec.txt if this is unexpected"))
      }
    } else {
      stop("Error: when reading a single cohort from multiple source files, need same number of id files")
    }
  } else { if(numfls!=1) { warning(paste("length of ID list was",numfls,"but only 1 input file")) } } 
  #### DETERMINE FILE DIMENSIONS ####
  num.sub <- length(cmb.ID.list) #ID.list)
  smp.szs <- sapply(ID.list,length)
  fil.ofs <- c(0,cumsum(smp.szs)) #note last element is the end of the last file
  num.snp <- length(snp.list)
  cat(paste(" found",num.sub,"samples and",num.snp,"markers\n"))
  
  # use 'pref' as the name of the big.matrix backing files for this cohort
  bck.fn <- paste(pref,"bckfile",sep="")
  des.fn <- paste(pref,"descrFile",sep="")
  # attempt to find amount of free memory using 'top()'
  top1 <- top(CPU=F); if(!is.null(top1)) { rGb <- top1$RAM$Free } else { rGb <- NA }
  if(!is.numeric(rGb)) { rGb <- 2 }
  return(import.big.data(input.fn=input.fn, dir=dir, 
                         row.names=snp.list, col.names=ID.list, 
                         pref=pref, delete.existing=delete.existing, 
                         ret.obj=ret.obj, verbose=verbose, 
                         dat.type=dat.type, ram.gb=rGb, hd.gb=max.mem))
}



# old.dat.to.big.matrix <- function(dir, grp=NA, snp.fn="snpNames.txt", sample.fn=NULL, 
#            input.fn=NULL, pref="LRR", delete.existing=T,ret.obj=F,DT=NULL,verbose=T)
# {
#   # import from a text (hopefully long format) datafile to a big.matrix
#   # return bigmatrix description 'object' or name of description file according to 'ret.obj'
#   ### CALIBRATE OPTIONS AND PERFORM CHECKS ###
#   max.mem <- 4096 # stupidly large
#   bafmode <- F # assume in LRR mode unless otherwise indicated by 'pref'
#   if(length(grep("BAF",toupper(pref)))>0) {
#     ## note these must match initial bash script that extracts data!
#     dat.file.suf <- "BAF.dat"; bafmode <- T
#   } else {
#     dat.file.suf <- "LRR.dat"
#   }
#   ## Define data types for big.matrix
#   dat.type <- double(1)
#   dat.typeL <- "double"  # label
#   dir <- validate.dir.for(dir,c("ano","big","col"),warn=F)
#   #### GET THE SHAPED INPUT FILE NAME(S) ####
#   if(all(is.null(input.fn))) {
#     ## try to automatically detect datafiles if not specified
#     if(is.data.tracker(DT)) {
#       if(!bafmode) { input.fn <- getSlot(DT,"shape.lrr",grps=grp) 
#       } else { input.fn <- getSlot(DT,"shape.baf",grps=grp) }
#       if(is.list(input.fn)) { input.fn <- unlist(input.fn,recursive=F) } #cat("CHECK: removed 1 level of list\n") }
#     } else {
#       if(!bafmode) { input.fn <- get.lrr.files(dir,grp,suffix=dat.file.suf)    
#       } else { input.fn <- get.lrr.files(dir,grp,suffix=dat.file.suf,baf=T) }
#     }
#   }
# # print(input.fn)
#   #### GET THE CORRESPONDING SAMPLE LIST(S) ####
#   if(all(!is.null(sample.fn))) {
#     sample.fn <- find.file(sample.fn,dir$ids)  # sample id file name
#     cat("Reading sample and snp lists from text files...\n")
#     cat(paste(" samples entered manually, reading samples from",sample.fn,"\n"))
#     ID.list <- lapply(sample.fn,readLines)  #readLines(sample.fn)
#   } else {
#     if(is.data.tracker(DT)) {
#       ## GET samples using data.tracker object 'DT'
#       if(all(!is.na(grp))) {
#         if(verbose) { cat(" getting ID list(s) for group",grp,"\n") }
#         # return file lists, listed by 'group'; separate lists for each file despite grp
#         ID.list <- unlist(getSlot(DT,"sub.samples",grps=grp,ret.obj=T),recursive=F)
#       } else {
#         ng <- getSlot(DT,"ngrps")
#         if(length(input.fn)==1 & ng==1) {
#           # all samples in 1 set
#           if(verbose) { cat(" 1 input file and 1 sample set only, getting IDs\n") }
#           ID.list <- list(unlist(getSlot(DT,"samples",ret.obj=T))) 
#         } else {
#           # separate list for each file listed by group sets
#           if(verbose) { cat(" getting all ID lists ignoring group [",ng," groups]\n",sep="") }
#           ID.list <- unlist(getSlot(DT,"sub.samples",ret.obj=T),recursive=F) 
#         }
#       }
#     } else {
#       ## GET samples using 'file.spec.txt' # deprecated
#       stop("if we need this code back, go to before march 19\n")
#     }
#   }
#   ##print(headl(ID.list))
#   cmb.ID.list <- paste(do.call("c",ID.list))
#   ##print(length(ID.list[[1]]))
#   numfls <- length(ID.list)
#   ### GET THE ORDERED SNP LIST ###
#   if(is.data.tracker(DT)) {
#     snp.fn <- getSlot(DT,"snps")
#   } else {
#     snp.fn <- find.file(snp.fn,dir$ano,dir)  # snp annotation file name
#   }
#   #print(head(cmb.ID.list))
#   cat(paste(" reading snps from",snp.fn,"\n"))
#   snp.list <- readLines(snp.fn) ##; print(head(snp.list))
#   ## Multiple possible input situations: 
#   ##   FOR LRR: unlist all input files and subsamples for selected grp
#   ##    n groups in study, 1 file each: should have input.fn length 1, id.list length 1
#   ##    n groups in study, m_n files each: should have input.fn length m_n, id.list length m_n
#   ##   FOR BAF: unlist all input files and subsamples for all grps
#   ##    n groups in study, 1 file each: should have input.fn length n, id.list length n
#   ##    n groups in study, m_n files each: should have input.fn length sum_i(m_n=i), id.list length sum_i(m_n=i)
#   ## DEBUG print(length(input.fn)); print(numfls); print(input.fn)
#   if(length(input.fn)>1) {
#     if(length(input.fn)==numfls) {
#       if(verbose) {
#         warning(paste("reading a single cohort from",numfls,"source files. Edit file.spec.txt if this is unexpected"))
#       }
#     } else {
#       stop("Error: when reading a single cohort from multiple source files, need same number of id files")
#     }
#   } else { if(numfls!=1) { warning(paste("length of ID list was",numfls,"but only 1 input file")) } } 
#   #### DETERMINE FILE DIMENSIONS ####
#   num.sub <- length(cmb.ID.list) #ID.list)
#   smp.szs <- sapply(ID.list,length)
#   fil.ofs <- c(0,cumsum(smp.szs)) #note last element is the end of the last file
#   num.snp <- length(snp.list)
#   cat(paste(" found",num.sub,"samples and",num.snp,"markers\n"))
#   cls <- num.sub; rws <- num.snp
#   cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
#   memory.estimate <- as.double((as.double(rws)*as.double(cls))/cells.per.gb)
#   if (memory.estimate > max.mem) {
#     cat("\nError: Insufficient disk space availability expected for this import method\n")
#     cat("Please free up some disk space and try again\n")
#     stop()
#   }
#   # use 'pref' as the name of the big.matrix backing files for this cohort
#   bck.fn <- paste(pref,"bckfile",sep="")
#   des.fn <- paste(pref,"descrFile",sep="")
#   ### DELETE EXISTING FILE IF HAS SAME NAME ###
#   if ((!des.fn %in% list.files(dir$big)) | delete.existing )
#   {
#     if(delete.existing & (des.fn %in% list.files(dir$big)))
#     {
#       dfn <- cat.path(dir$big,des.fn)
#       cat("\n deleting",dfn,"\n")
#       unlink(dfn)
#     } else {
#       #all clear, no files already exist with same name
#     }
#   } else {
#     cat(paste("\nWarning: Big matrix description file",des.fn,"already exists in",dir$big,"\n"))
#     cat("You may wish to delete, rename or move this file, or use option 'delete.existing'=T, before re-running this script\n")
#     #stop()
#   }
#   ### MAIN IMPORT LOOP ###
#   cat("\nCreating big matrix object to store group data")
#   cat("\n predicted disk use: ",round(memory.estimate,1),"GB\n")
#   bigVar <- big.matrix(nrow=rws,ncol=cls, backingfile=bck.fn, dimnames=list(snp.list,cmb.ID.list),
#                        type=dat.typeL, backingpath=dir$big, descriptorfile=des.fn)
#   for(ff in 1:numfls) {
#     #ifn <- cat.path(dir$col,input.fn[ff],must.exist=T)
#     # create file name depending on whether in baf or lrr directory
#     if(!bafmode) {  
#       ifn <- cat.path(dir$col,input.fn,must.exist=T)
#     } else {  
#       ifn <- cat.path(dir$baf.col,input.fn,must.exist=T)  
#     }
#     #test file type if matrix
#     if(file.ncol(ifn[ff])>1) { input.is.vec <- F } else { input.is.vec <- T }
#     if(!input.is.vec) {
#       frm <- check.text.matrix.format(fn=ifn[ff],ncol=smp.szs[ff],header=ID.list[[ff]],row.names=snp.list)
#       if(frm$rname & !frm$match) { file.rn <- character(cls) } # recording rownames as we go
#     }
#     dat.file <- file(ifn[ff])
#     dir.ch <- .Platform$file.sep
#     open(con=dat.file,open="r")
#     cat(paste(" opening connection to ",c("matrix","long")[1+input.is.vec],
#               " format datafile (",ff,dir.ch,numfls,"): ",basename(ifn[ff]),"\n",sep=""))
#     cat("\nLoading text data into big matrix object:\n")
#     nxt.rng <- (fil.ofs[ff]+1):(fil.ofs[ff+1])
#     cls <- length(nxt.rng)
#     if(!input.is.vec)
#     {
#       ## read from matrix format tab file
#       twty.pc <- round(rws/5)
#       for (cc in 1:rws) {
#         if ((cc %% twty.pc)==0)  { fl.suc <- flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
#         loop.tracker(cc,rws)
#         next.line <- readLines(dat.file,n=1)
#         next.row <- strsplit(next.line,"\t",fixed=T)[[1]]
#         if (cc==1 & frm$header) { 
#             # need to go to next line to avoid reading header as data
#             next.line <- readLines(dat.file,n=1); next.row <- strsplit(next.line,"\t",fixed=T)[[1]]
#         }
#         if (frm$rnames) {
#           if(frm$match) {
#             lbl <- next.row[1]; bigVar[match(lbl,snp.list),nxt.rng] <- next.row[-1]
#           } else {
#             file.rn[nxt.rng[cc]] <- next.row[1]; bigVar[cc,nxt.rng] <- next.row[-1]
#           }
#         } else {
#           bigVar[cc,nxt.rng] <- next.row
#         }
#       }
#     } else {
#       ## read from (long) vector format tab file
#       twty.pc <- round(cls/5)
#       for (cc in 1:cls) {
#         loop.tracker(cc,cls)
#         if ((cc %% twty.pc)==0)  { fl.suc <- flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
#         bigVar[,(cc+fil.ofs[ff])] <- as(readLines(dat.file,n=rws),dat.typeL)
#       }
#     }
#     close(dat.file)
#   }
#   ### FINISH UP, RETURN BIGMATRIX DESCRIPTION ###
#   # in the case of different rownames found in matrix, then show following warning text:
#   if(!input.is.vec) {
#     if(frm$rname & !frm$match) {
#       ofn <- cat.path(dir$ano,pref=pref,"_file_rowname_list_check_this.txt")
#       warning(paste("rownames didn't match what was in file, check the list in the file at:\n ",ofn))
#       writeLines(paste(file.rn),con=ofn) ; cat("\n file preview:\n"); print(head(file.rn,10)); cat("\n")
#     }
#   }
#   cat("\n")
#   cat(paste(" created big.matrix description file:",des.fn,"\n"))
#   cat(paste(" created big.matrix backing file:",bck.fn,"\n"))
#   if(ret.obj) {
#     return(describe(bigVar))
#   } else {
#     return(des.fn)
#   }
#   cat("...complete!\n")
# }


make.dir <- function(dir.base="plumbCNV_LRRQC",dir.raw="LRRQC_Raw_Files",no.raw=F)
{
  #create a directory structure for plumbCNV to work upon starting with a rawdata and base directory
  #raw data:  dir.raw <- "/ipswich/data/Immunochip/FinalReports/"
  #location of support files directory dir.sup <- "/chiswick/data/store/metabochip/PLINK/"
  #location to write summaries of call rate: dir.base <- "/chiswick/data/ncooper/ImmunochipReplication/"
  #make sure each ends in a slash
  dir.force.slash <- reader:::dir.force.slash # use internal function from 'reader'
  paste.dr <- function(x,...) { dir.force.slash(paste(x,...)) }
  if(no.raw) { dir.raw <- NULL }
  dir.raw <- dir.force.slash(dir.raw)
  ##dir.sup <- dir.force.slash(dir.sup)
  dir.base <- dir.force.slash(dir.base)
  
  dir.res <- paste.dr(dir.base,"RESULTS",sep="")
  dir.cr <- paste.dr(dir.base,"SNPQC",sep="")
  dir.cr.plk <- paste.dr(dir.cr,"PLINK",sep="")
  #location of big matrix main data files (contains LRR data for all samplessnps)
  dir.big <- paste.dr(dir.base,"BIGMATRICES",sep="")
  #location of annotation files directory
  dir.ano <- paste.dr(dir.base,"ANNOTATION",sep="")
  # gc wave
  dir.gc <- paste.dr(dir.base,"GCDATA",sep="")
  #qc directories
  dir.qc <- paste.dr(dir.base,"SAMPLEQC",sep="")
  dir.qc.lrr <- paste.dr(dir.qc,"MEAN_DLRS_GC",sep="")
  dir.qc.cab <- paste.dr(dir.qc,"CHRAB",sep="")
  dir.qc.pl <- paste.dr(dir.qc,"PLATES",sep="")
  #dir.qc.pc <- paste(dir.qc,"PCCOMPARISON",sep="")
  dir.ind <- paste.dr(dir.qc,"EXAMPLES",sep="")
  #scripts directory
  dir.scr <- paste.dr(dir.base,"SCRIPTS",sep="")
  # sample exclude directory
  dir.excl <- paste.dr(dir.ano,"SAMPLE_EXCLUDE",sep="")
  dir.sort <- paste.dr(dir.ano,"SAMPLE_SORT",sep="")
  dir.excl2 <- paste.dr(dir.ano,"SNP_EXCLUDE",sep="")
  dir.sort2 <- paste.dr(dir.ano,"SNP_SORT",sep="")
  # principle components directories
  dir.pc <- paste.dr(dir.base,"PCA",sep="")
  #penn cnv files
  dir.cnv <- paste.dr(dir.base,"PENNCNV",sep="")
  dir.cnv.raw <- paste.dr(dir.cnv,"PENNRAWFILES",sep="")
  dir.cnv.pen <- paste.dr(dir.cnv,"PENNOUTPUT",sep="")
  dir.cnv.qc <- paste.dr(dir.base,"CNVQC",sep="")
  dir.cnv.qsub <- paste.dr(dir.cnv,"PENNCLUSTER",sep="")
  dir.lrr.dat <- paste.dr(dir.base,"LRRDATA",sep="")
  dir.baf <- paste.dr(dir.base,"BAFDATA",sep="")
  dir.col <- paste.dr(dir.lrr.dat,"RAWDATA",sep="")
  dir.baf.col <- paste.dr(dir.baf,"RAWDATA",sep="")
  dir.ids <- paste.dr(dir.lrr.dat,"SAMPLEFILES",sep="")
  ## now search the local name space for all dir.xxx
  if ("dir" %in% ls()) { rm(dir) }
  all.locs <- ls()[grep("dir.",ls())]
  nloc <- length(all.locs)
  dir <- vector("list", nloc)
  inval <- NULL
  for (cc in 1:nloc) { 
    nxt <- get(all.locs[cc])
    if(!is.function(nxt)) {
      dir[[cc]] <- nxt
    } else { inval <- c(inval,cc) }
  }
  if(!is.null(inval)) { dir <- dir[-inval] ; all.locs <- all.locs[-inval] }
  #dir <- dir.force.slash(dir) ; dir <- as.list(dir)
  names(dir) <- substr(all.locs,5,nchar(all.locs))
  return(dir)
}


lrr.dat.to.big.matrix  <- function(dir,ngrp=1,import.LRR="LRR.dat",
                                   sample.ids=NULL,snp.fn="snpNames.txt",
                                   pref="LRR",delete=T,DT=NULL) 
{
  # converts set of LRR long format/matrix *.dat files into big.matrix objects
  if(is.data.tracker(DT)) {
    dt.mode <- T 
    import.LRR <- NULL
  } else { 
    dt.mode <- F 
    if(!is.file(import.LRR,dir$col,dir)) {
      warning("did not find shaped lrr datafile from data tracker object, will auto search instead")
      import.LRR <- NULL # file not there, so allow for auto searching in plain.ve>>>ix()
    } else {
      import.LRR <- find.file(import.LRR,dir$col)
    }
  }
#  cat(" no pre-existing LRR big.matrices, or 'delete=T', importing now from text file\n")
  ## loop through for each group ##
  descr.list <- vector("list", ngrp)
  cat("\nRunning import for",ngrp,"cohort(s)\n")
  for (grp in 1:ngrp) {
    if(ngrp>1) { pre <- paste(pref,grp,sep="") } else { pre <- pref }
    if(dt.mode) {
      descr.list[[grp]] <- dat.to.big.matrix(dir, grp, DT=DT, input.fn=NULL, sample.fn=NULL,
                                             pref=pre, delete.existing=delete, ret.obj=F)
    } else {
      descr.list[[grp]] <- dat.to.big.matrix(dir, grp, input.fn=import.LRR, sample.fn=sample.ids, 
                                             pref=pre, delete.existing=delete, ret.obj=F)
    }
  }
  return(descr.list)
}


baf.dat.to.big.matrix <- function(dir,BAF.fn="BAFdescrFile",import.BAF="BAF.dat",
                                  sample.ids=NULL,snp.fn="snpNames.txt",
                                  delete=T,DT=NULL,pref="BAF") 
{
  # converts a BAF long format/matrix *.dat file(s) into a big.matrix object
  if( delete | !file.exists(cat.path(dir$big,BAF.fn))) {
    if(is.data.tracker(DT)) { dt.mode <- T } else { dt.mode <- F }
    if(dt.mode) {
      import.BAF <- as.character(unlist(getSlot(DT,"shape.baf")))
    }
    if(!is.file(import.BAF,dir$baf.col,dir)) {
      warning("did not find shaped baf datafile from data tracker object, will auto search instead")
      import.BAF <- NULL # file not there, so allow for auto searching in plain.ve>>>ix()
    } else {
      import.BAF <- find.file(import.BAF,dir$baf.col,dir)
    }
    cat(" no pre-existing BAF big.matrix, or 'delete=T', importing now from text file\n")
    if(dt.mode) {
      descr <- dat.to.big.matrix(dir, DT=DT, input.fn=import.BAF, sample.fn=NULL, snp.fn=snp.fn,
                     pref=pref, delete.existing=delete,ret.obj=F)
    } else {
      descr <- dat.to.big.matrix(dir, input.fn=import.BAF, sample.fn=sample.ids, snp.fn=snp.fn,
                                 pref=pref, delete.existing=delete,ret.obj=F) 
    }
  } else {
    cat(" using existing BAF big.matrix file:",BAF.fn,"\n")
    descr <- BAF.fn
  }
  return(descr)
}


import.DATA.big <- function(DT,...,dir=NULL) {
  ## If using the datatracker, allows import of baf and lrr shaped data to
  # big matrix format with minimal input requirements.
  if(!is.data.tracker(DT)) {
    stop("Error: 'import.DATA.big' requires 'DT' to be a data.tracker object\n")
  } else {
    dir <- getSlot(DT,"dir")
  }
  DTs.required <- c("shape.lrr","ngrps","shape.baf","snps","samples","subsamples")
  is.valid <- req.valid(DT,DTs.required)
  if(!is.valid) { stop("Error: file required for 'import.DATA.big' invalid in data.tracker\n")}
  BAF.fn <- baf.dat.to.big.matrix(dir,BAF.fn="BAFdescrFile",DT=DT,...) 
  DT <- setSlot(DT,big.baf=BAF.fn,warns=F)
  #write.data.tracker(DT,dir) # backup file in case of a crash or something
  ##sample.fn <- getSlot(DT,"shape.lrr")
  ngrp <- getSlot(DT,"ngrps")
  descr.list <- lrr.dat.to.big.matrix(dir,ngrp=ngrp,DT=DT,...)
  if(is.character(descr.list)) {
    ofn <- cat.path(dir$big,rmv.ext(descr.list),ext="RData")
    descr <- describe(getBigMat(descr.list))
    save(descr,file=ofn)
    descr.list <- ofn
  }
  DT <- setSlot(DT,big.lrr=descr.list,grps=1:ngrp,proc.done=1,warns=T)
  #write.data.tracker(DT,dir)
  return(DT)
}


confirm.corrected.filename <- function(dir,num.pcs=9)
{
  # expects PC corrected bigmatrix to have a file name of a certain format,
  # given the number of PCs corrected for. If this expected filename is not
  # found, then it will look for the nearest thing in the relevant directory
  LRR.fn <- paste("describePCcorrect",num.pcs,".RData",sep="")
  if(!LRR.fn %in% list.files(dir$big)) {
    anyfl <- grep("describePCcorrect",list.files(dir$big))
    if(length(anyfl)>0) {
      thefl <- list.files(dir$big)[anyfl][1]
      cat("\nWarning: couldn't find:\n ",LRR.fn,"but did find:\n ",thefl,"in:\n ",dir$big,"\n")
      cat("will continue using this file, if this is undesired please ctrl-C and \n")
      cat("review the file name passed to 'prepare.penn.cnv.samples()' before continuing.\n")
      return(thefl)
    } else {
      cat("No PC corrected big matrix 'R.Data' files found in:\n ",dir$big,"\n")
      cat("Please review the big matrix directory before re-running this function\n")
      stop()
    }
  } else {
    return(LRR.fn)
  }
}


prepare.one.sample <- function(lll, set.offset, baf.col.sel, lrr.col.sel,
                        cL, cB, lrr.txt, baf.txt, grand.rn, 
                        lrr.row.sel, baf.row.sel, bigLRR, bigBAF, next.dir) {
  # all the commands to process a single sample file.
  # mainly packaged into a function to tidy up the code in the main function
  loc <- set.offset+lll
  baf.sel <- baf.col.sel[loc]
  lrr.sel <- lrr.col.sel[loc]
  sampleid <- cL[lrr.sel]; s2 <- cB[baf.sel] 
  if(sampleid!=s2) { stop("Error: col indexes failed") }
  header <- paste("Name",paste(sampleid,lrr.txt,sep=""),
                  paste(sampleid,baf.txt,sep=""),sep="\t")
  to.write <- paste(grand.rn,bigLRR[lrr.row.sel,lrr.sel],
                    bigBAF[baf.row.sel,baf.sel],sep="\t")
  conn <- file(paste(next.dir,"/",sampleid,".txt",sep=""),"w")
  writeLines(header,con=conn)
  writeLines(to.write,con=conn)
  close(conn)
}


prepare.penncnv.samples <- function(dir,LRR.fn="",BAF.fn="BAFdescrFile",num.pcs=9,
                                    low.ram=T,n.cores=3,cores=3,verbose=F,add.name.files=T,
                                    relative=F,...) 
{
  # penn-cnv requires the LRR and BAF data for each sample to be in a separate file
  # this function imports the BAF data if not present, then writes individual LRR,BAF files
  # for all samples
  if(n.cores>1) { must.use.package("multicore") ; multi <- T } else { multi <- F }
  dir <- validate.dir.for(dir,c("cnv","big","cnv.raw"))
  # default PENNCNV header suffixes
  baf.txt <- ".B Allele Freq"  ;  lrr.txt <- ".Log R Ratio"
  cat("\nPreparing input files for PENN-CNV from LRR and BAF\n")
  if(is.null(LRR.fn)) { LRR.fn <- "" }
  if(LRR.fn=="") {  LRR.fn <- confirm.corrected.filename(dir,num.pcs) }
  if(!file.exists(cat.path(dir$big,LRR.fn))) {
    stop(paste("Error: could not find LRR big.matrix file:",LRR.fn))
  } else {
    cat(" using LRR big.matrix file:",LRR.fn,"\n")
  }
  if(!file.exists(cat.path(dir$big,BAF.fn))) {
    stop(paste("Error: could not find BAF big.matrix file:",BAF.fn))
  } else {
    cat(" using BAF big.matrix file:",BAF.fn,"\n")
  }
  bigBAF <- getBigMat(BAF.fn, dir$big)
  bigLRR <- getBigMat(LRR.fn, dir$big)
  print.big.matrix(bigBAF,"bigBAF") ; print.big.matrix(bigLRR,"bigLRR")
  rB <- rownames(bigBAF); cB <- colnames(bigBAF)
  rL <- rownames(bigLRR); cL <- colnames(bigLRR)
  n <- ncol(bigLRR)
  if(cores[1]>1) { ndirs <- cores[1] } else { ndirs <- ceiling(n/1000) }
  dir.sizes <- (n %/% ndirs) + c(rep(1,n %% ndirs),c(rep(0,ndirs-(n %% ndirs))))
  if(length(dir.sizes)>10) {
    sz.txt <- paste("N=",names(table(dir.sizes))," (",table(dir.sizes),")",sep="",collapse="; ")
  } else { sz.txt <- paste(dir.sizes,collapse=",") }
  cat("\nWill create",ndirs,"penn folders with sizes:",sz.txt,"\n")
  ##BAF.nms <- paste(cB,baf.txt,sep="") #used where??
  ##LRR.nms <- paste(cL,lrr.txt,sep="")  #used where??
  baf.col.sel <- match(cL,cB) ; cat(length(baf.col.sel[!is.na(baf.col.sel)]),"cols from LRR in BAF\n")
  baf.row.sel <- match(rL,rB) ; cat(length(baf.row.sel[!is.na(baf.row.sel)]),"rows from LRR in BAF\n")
  lrr.col.sel <- which(!is.na(baf.col.sel))
  mln <- length(which(is.na(baf.col.sel))); if(mln>0) { cat("*warning:") }
  cat(" BAF missing",mln,"samples in LRR file\n")
  lrr.row.sel <- which(!is.na(baf.row.sel))
  mln <- length(which(is.na(baf.row.sel))); if(mln>0) { cat("*warning:") }
  cat(" BAF missing",mln,"snps in LRR file\n")
  grand.rn <- rL[lrr.row.sel]
  ff <- FALSE # initialise failure flag
  if(!all(grand.rn==rB[baf.row.sel])) 
  { warning("row name selections don't match") ; ff <- T }
  tot <- n.to.process <- sum(dir.sizes)
  if(n.to.process!=length(lrr.col.sel)) {   
    ff <- T ;  warning("number of samples per directory and number of directories\n",
    " does not imply the number of common columns in the BAF,LRR sample matrices\n") }
  cat(" generating individual sample LRR/BAF files in Penn-CNV format:\n")
  if(multi) { 
    at.once <- n.cores
    i.nutin <- vector("list",at.once)
    #loop.tracker(1,tot) # write the percentages header as multi will skip to >'1'
  } else {
    at.once <- 1 
  }
  if(!ff){
    st.t <- proc.time()
    cat("minutes remaining: ")
    for (nn in 1:length(dir.sizes))
    {
      next.dir <- paste(dir$cnv.raw,"p",nn,sep="")
      dir.create(next.dir)
      set.offset <- sum(dir.sizes)-sum(dir.sizes[nn:ndirs])
      next.mm <- 1
      for (mm in 1:dir.sizes[nn])
      {
        if(multi) {
          if(mm < next.mm) { next } else { next.mm  <- min(dir.sizes[nn]+1,next.mm+at.once) }
          nutin <- i.nutin
          for(lll in mm:(next.mm-1)) {
            nutin[[lll]] <- multicore::parallel(prepare.one.sample(lll=lll, set.offset=set.offset, 
                        baf.col.sel=baf.col.sel, lrr.col.sel=lrr.col.sel,
                        cL=cL, cB=cB, lrr.txt=lrr.txt, baf.txt=baf.txt, grand.rn=grand.rn, 
                        lrr.row.sel=lrr.row.sel, baf.row.sel=baf.row.sel, 
                        bigLRR=bigLRR, bigBAF=bigBAF, 
                        next.dir=next.dir))  #
          }
          noob <- collect(nutin[mm:(next.mm-1)],wait=T)
        } else {
          prepare.one.sample(lll=mm, set.offset=set.offset, 
                      baf.col.sel=baf.col.sel, lrr.col.sel=lrr.col.sel,
                      cL=cL, cB=cB, lrr.txt=lrr.txt, baf.txt=baf.txt, grand.rn=grand.rn, 
                      lrr.row.sel=lrr.row.sel, baf.row.sel=baf.row.sel, 
                      bigLRR=bigLRR, bigBAF=bigBAF, 
                      next.dir=next.dir)
        }
        n.to.process <- max(0,(n.to.process-at.once))
        loop.tracker((tot-n.to.process),tot,st.t)
      }
      if(low.ram) {
        test <- runif(1) < ((1/500)*dir.sizes[nn])
        # flush and reload bigmatrices to clear RAM
        if(test) {
          rm(bigBAF); rm(bigLRR); gc()
          bigBAF <- getBigMat(BAF.fn, dir$big,verbose=F)
          bigLRR <- getBigMat(LRR.fn, dir$big,verbose=F)
        }
      }
    }
  }
  subj.in.each <- penn.name.files(dir,write.files=add.name.files,
                                  ret.fn=F,relative=relative)
  return(subj.in.each)
}


prepare.penncnv.markers <- function(snp.info,lrr.descr, baf.descr="BAFdescrFile", dir, ucsc="hg18", n.cores=1)
{
  #penn-cnv requires a pfb file and a gcm file to run. These show the average BAF and
  # the average GC percentage for every marker, as well as the locations and labels
  must.use.package(c("IRanges","genoset"),bioC=T)
  dir <- validate.dir.for(dir,c("cnv","big"),warn=F)
  baf.out.fn <- cat.path(dir$cnv,"BAF.pfb")
  gc.out.fn <- cat.path(dir$cnv,"marker.gcm")
  bigMat2 <- getBigMat(lrr.descr,dir)
  snp.list <- rownames(bigMat2)
  samp.list <- colnames(bigMat2)
  rm(bigMat2)
  gc.dat <- get.gc.markers(dir=dir,ret="bio",ucsc=ucsc,snp.info=snp.info, n.cores=n.cores) # snp.info can be made locally?
  row.selectSI <- match(snp.list,rownames(snp.info))
  row.selectGC <- match(snp.list,rownames(gc.dat))
  ## THIS IS WHERE I FOUND OUT THAT GC RETURNS MISSIN VALUES##
  #gc.file <- as.data.frame(gc.dat)[row.select,c("names","space","start","gc")]
  # NB: manual extraction as 'as.data.frame' method seems to fail
  gc1 <- rownames(snp.info); gc2 <- snp.info$space; gc3 <- start(snp.info)
  gc4 <- gc.dat$gc[match(rownames(snp.info),rownames(gc.dat))]
  gc.file <- data.frame(Name=gc1,Chr=gc2,Position=gc3,GC=gc4)
  gc.mean <- mean(gc.file$GC,na.rm=T)
  if(gc.mean<=1) {
    # convert to decimal percentage
    gc.file$GC <- round(gc.file$GC*100)
    gc.mean <- mean(gc.file$GC,na.rm=T)
  }
  gc.file$GC[is.na(gc.file$GC)] <- gc.mean
  gc.file <- gc.file[row.selectSI,]
  cat("~wrote SNP GC averages to file:\n ",gc.out.fn,"\n")
  # Name<TAB>Chr<TAB>Position<TAB>GC
  # markername<TAB>chromosome<TAB>position<TAB>meanGC
  write.table(gc.file,file=gc.out.fn,row.names=F,col.names=T,quote=F,sep="\t")
  marker.means <- get.baf.markers(baf.descr,samp.list,dir)
  v.select <- match(snp.list,names(marker.means))
  baf.file <- gc.file[,1:3]
  baf.file[[4]] <- marker.means[v.select]
  colnames(baf.file) <- c("Name","Chr","Position","PFB")
  cat(paste("~wrote file baf means for each SNP to file:\n ",baf.out.fn),"\n")
  # header: "Name","Chr","Position","PFB"              #
  # markername<TAB>chromosome<TAB>position<TAB>meanBAF #
  write.table(baf.file,file=baf.out.fn,row.names=F,col.names=T,quote=F,sep="\t")
  out <- list(gc.out.fn,baf.out.fn)
  names(out) <- c("gc","baf")
  return(out)
}


get.penn.cmd <- function(dir,ndir=1,gc.out.fn="marker.gcm",baf.out.fn="BAF.pfb",use.penn.gc=T,
                         penn.path="/usr/local/bin/penncnv/",pref="output",relative=F,hmm="hh550.hmm")
{
  # this function generates calls to runn penn-cnv for 'ndir' directories created by plumbcnv
  # can run with either absolute or relative file paths. marker.fn is a list
  #  containing elements 'baf' and 'gc' which are the files from prepare.penncnv.markers()
  dir <- validate.dir.for(dir,c("cnv","cnv.raw","cnv.pen"),warn=F)
  penn.cmd <- paste("perl",cat.path(penn.path,"detect_cnv.pl",must.exist=T),"-test")
  arg.nms <- c("hmm","pfb","list","log","out") ; num.to.chk <- c(1,2)
  if(use.penn.gc) { arg.nms <- c(arg.nms,"gcmodel"); num.to.chk <- c(num.to.chk,6) }
  penn.args.list <- vector("list",length(arg.nms))
  names(penn.args.list) <- arg.nms
  penn.args.list[[1]] <- paste(penn.path,"lib/",hmm,sep="")
  penn.args.list[[2]] <- find.file(baf.out.fn,dir$cnv,dir)
  if(use.penn.gc) { penn.args.list[[6]] <- find.file(gc.out.fn,dir$cnv,dir) }
  for (cc in num.to.chk) {
    nf <- penn.args.list[[cc]]
    if(!file.exists(nf)) { warning(paste("file",nf,"from PennCNV command not found")) }
  }
  #print(ndir)
  penn.call <- character(ndir)
  for (dd in 1:ndir) {
    penn.args.list[[3]] <- cat.path(dir$cnv.raw,paste("p",dd,"names.txt",sep=""))
    penn.args.list[[4]] <- cat.path(dir$cnv.pen,paste(pref,"p",dd,".log",sep=""))
    penn.args.list[[5]] <- cat.path(dir$cnv.pen,paste(pref,"p",dd,".cnv",sep=""))
    penn.args <- NULL
    for (cc in 1:length(penn.args.list)) {
      penn.args <- paste(penn.args,paste("-",names(penn.args.list)[cc]," ",penn.args.list[[cc]]," ",sep=""),sep="")
    }
    penn.call[dd] <- paste(penn.cmd,penn.args)
  }
  if(relative) {
    # remove dir$cnv to make dirs relative rather than absolute (neater)
    penn.call <- gsub(dir$cnv,"./",penn.call)
  }
  return(penn.call)
}  


process.penn.results <- function(dir,sample.info=NULL,penn.path="/usr/local/bin/penncnv/",ucsc="hg18",min.sites=10,rare.olp=.5,rare.pc=1,
                                 baf.fn="BAF.pfb",cnv.qc=T,rare.qc=T,plate.qc=T,pval=0.05,del.rate=0.4,dup.rate=0.18,thr.sd=3,plate.thr=3,rmv.low.plates=F,
                                 raw.main="raw",plink.out="plink.cnv",rmv.bad.reg=T,bad.reg.fn="badRegions.txt",hide.plink.out=T,verbose=F) 
{
  # This function runs the penn cnv output through a series of filters from penn scripts,
  # plink commands, both utilising bash commands, and also R-functions, to end up with
  # the final QC-ed set of deletions, duplications, CNVs, both common and rare.
  dir <- validate.dir.for(dir,c("cnv.plk","cnv.qc","excl"))
  ##script names##
  nn <- vector("list",10) # to store ignored plink output
  scn.scrpt <- "scan_region.pl"; cln.scrpt <- "clean_cnv.pl"
  conv.scrpt <- "penncnv_to_plink.pl"; force.scripts <- c(scn.scrpt,cln.scrpt)
  if(!file.exists(penn.path)) { stop("Error: PennCNV path incorrect or software not installed") }
  penn.available <- list.files(penn.path)
  missing.ones <- !force.scripts %in% penn.available
  if(any(missing.ones)) {
    cat("PennCNV script:",force.scripts[missing.ones],
               "not found in:",penn.path,". Please download from PennCNV website\n")   
    stop("PennCNV script file(s) missing!")
  }
  setwd(dir$cnv.qc)
  commands <- character(28)
  # annotation to help remove centromeric, telomeric and immunoglobic overlap regions
  make.bad.region.file(dir,fn=bad.reg.fn,ucsc=ucsc)
  # make plink family file for whole dataset 
  if(is.null(sample.info)) { sample.info <- read.sample.info(dir) }
  sample.info <- validate.samp.info(sample.info,dir=dir,QC.update=T,verbose=F)
  nsamps <- length(which(sample.info$QCfail==0)) # used below for filtering
  # DEFINE RARE #
  if (rare.olp>1) { rare.olp <- rare.olp/100 } # in case entered as % convert to decimal
  if (rare.pc>=1) { rare.pc <- rare.pc/100 } # % of incidence below which is rare
  rare.cut <- round(rare.pc*nsamps) # number of samps this corresponds to
  # make plink family file
  make.fam.file(sample.info[sample.info$QCfail==0,],dir=dir,out.fn=cat.path("",plink.out,ext="fam"))
  raw.fn <- cat.path("",raw.main,ext="cnv")
  cat.arg <- cat.path(dir$cnv.pen,"outputp*",ext="cnv")
  penn.cmb <- paste("cat",cat.arg,">",raw.fn)
  commands[1] <- paste(penn.cmb)
  # get pure IDs by removing directory and file extension
  #rmv.dir.penn.cnv.file(raw.fn,append="",ext=T,verbose=T) # remove directory garbage from penn.cnv file
  cmd <- paste("perl ",penn.path,cln.scrpt,sep="")
  baf.fn <- cat.path(dir$cnv,baf.fn)
  cur.pen <- cat.path(fn=raw.fn,suf=".merge1",ext="cnv")
  args <- paste("combineseg --signalfile ",baf.fn," --fraction 0.2 --bp ",raw.fn," > ",cur.pen,sep="")
  commands[2] <- paste(cmd,args)
  prev.pen <- cur.pen; cur.pen <- cat.path(fn=raw.fn,suf=".merge2",ext="cnv")
  args <- paste("combineseg --signalfile ",baf.fn," --fraction 0.2 --bp ",prev.pen," > ",cur.pen,sep="")
  commands[3] <- paste(cmd,args)
  cmd <- paste("perl ",penn.path,scn.scrpt,sep="")
  cnv.in.bad <- cat.path("",bad.reg.fn,pre="cnvs_in_")
  args <- paste(cur.pen,bad.reg.fn,"-minqueryfrac 0.5 > ",cnv.in.bad)
  commands[4] <- paste(cmd,args)
  prev.pen <- cur.pen; cur.pen <- cat.path("",cur.pen,suf=".rmv.bad",ext="cnv")
  cat(" merging adjacent cnvs with <20% gaps between them ('clean_cnv.pl combineseg')..\n")
  #j0 <- proc.time()[3]
  #j1 <- proc.time()[3]; print(j1-j0)
  # plumbCNV() will remove any existing cnv exclusion file from SAMPLE_EXCLUDE prior to this
  bfn <- cat.path(dir$excl,"BadReg.txt"); brf <- FALSE # default no bad found
  if(rmv.bad.reg) {
    cat(" removing CNVs with >50% overlap with immunoglobin, telomeric and centromeric regions\n")
    commands[5] <- paste("fgrep -v -f",cnv.in.bad,prev.pen,">",cur.pen)
  } else {
    cat(" NB: did not remove CNVs in known low quality call regions (immunoglobin, centromere, telomere)\n")
    commands[5] <- paste("cat",prev.pen,">",cur.pen)
  } 
  # need to run these commands before processing the sample ids excluded
  nn[[1]] <- sapply(commands[1:5],system,intern=hide.plink.out, ignore.stderr=hide.plink.out) 
  if(rmv.bad.reg) {
    tmp <- rmv.dir.penn.cnv.file(cnv.in.bad,append=".tmp",ext=T,verbose=F)
    p.file <- read.penn.cnv.file(tmp,readtable=T)
    if(nrow(p.file)>0) {
      pcoln <- find.id.col(p.file,rownames(sample.info),ret="col")$col
      bad.reg.ids <- p.file[,pcoln]
      brf <- T
    }
    unlink(tmp) # delete temporary file
  } 
  # convert to plink format
  cat(". done\n")
  if(conv.scrpt %in% list.files(penn.path)) {
    cmd <- paste("perl ",penn.path,conv.scrpt,sep="")
    args <- paste("-i ",cur.pen," -o ",plink.out," -c 1 ",sep="")
    commands[6] <- paste(cmd,args)
    nn[[2]] <- system(commands[6],intern=hide.plink.out)
  } else {
    # do it in R instead (slower)
    cat(" converting penn-cnv file to plink format\n")
    convert.penn.to.plink(penn.in=cur.pen,plink.out=plink.out)
    commands[6] <- "# conversion to plink format was done in R"
  }
  # remove directory garbage from plink file
  rmv.dir.plink.file(plink.out,append="",ext=T,verbose=verbose)
  if(brf) {
    # bad region CNVs were removed, check if any samples should be removed as result
    still.cnv.ids <- read.plink.file(plink.out,readtable=T)$IID
    bad.reg.ids <- bad.reg.ids[!bad.reg.ids %in% still.cnv.ids]
    writeLines(paste(bad.reg.ids),con=bfn)
    cat(" wrote sample ids with no remaining CNVs after filtering bad regions to:",bfn,"\n")
  }
  #j2 <- proc.time()[3]; print(j2-j1)
  commands[7] <- "# directories and file extensions removed in R"
  commands[8] <- "# plink .fam file created in R"
  # make map file in plink
  commands[9] <- paste("plink --cnv-list ",plink.out," --cnv-make-map --out ",rmv.ext(plink.out),sep="")
  # check for overlapping regions (duplicates)
  commands[10] <- paste("plink --cfile",rmv.ext(plink.out),"--cnv-check-no-overlap --allow-no-sex")
  nn[[3]] <- sapply(commands[9:10],system,intern=hide.plink.out,USE.NAMES=F)
  #j3 <- proc.time()[3]; print(j3-j2)
  ovlp.fn <- cat.path(dir$cnv.qc,paste(rmv.ext(plink.out),".cnv.overlap",sep=""))
  if(file.exists(ovlp.fn)) {
    nr <- file.nrow(ovlp.fn)
    if(nr>1) { warning("Found ",nr-1," overlapping CNV segments within same sample(s). See file\n",
                       ovlp.fn,"\nWARNING! This is potentially a critical problem, please check the data for a systematic error")}
  }
  filter.hifreq.cnv.samples(plink.out, nsamps=nsamps, rem.file=plink.out, rare="",
                            dir=dir, thr.sd=thr.sd, QC.ON=cnv.qc, pval=pval,
                            ndels=NULL, ndups=NULL, del.rate=del.rate, dup.rate=dup.rate)
  #j4 <- proc.time()[3]; print(j4-j3)
  commands[11] <- "# remove samples with too many CNVs in R"
  suf.for.rmv <- "_CNVcnt"
  plink.rmv <- cat.path(dir="",plink.out,suf=suf.for.rmv,ext="cnv")
  dupR.file <- "rareDUP"
  delR.file <- "rareDEL"
  # rare dels
  commands[12] <- paste("plink --cfile",rmv.ext(plink.rmv),
                        "--cnv-del --cnv-sites",min.sites,"--cnv-drop-no-segment --cnv-freq-exclude-above",rare.cut,"--cnv-overlap",rare.olp,"--cnv-write --out",delR.file,"--allow-no-sex")
  # rare dups
  commands[13] <-paste("plink --cfile",rmv.ext(plink.rmv),
                       "--cnv-dup --cnv-sites",min.sites,"--cnv-drop-no-segment --cnv-freq-exclude-above",rare.cut,"--cnv-overlap",rare.olp,"--cnv-write --out",dupR.file,"--allow-no-sex")
  nn[[4]] <- sapply(commands[12:13],system,intern=hide.plink.out,USE.NAMES=F)
  
  suf.for.rmv <- "_DUPrcnt"
  plink.rmv2 <- cat.path(dir="",plink.rmv,suf=suf.for.rmv,ext="cnv")
  filter.hifreq.cnv.samples(dupR.file, nsamps, rem.file=plink.rmv, rare="DUP", dir=dir, 
                            thr.sd=thr.sd, ndels=NULL, ndups=NULL, 
                            del.rate=del.rate, dup.rate=dup.rate, QC.ON=rare.qc, pval=pval)
  commands[14] <- "# remove samples with too many rare deletions or duplications in R"
  suf.for.rmv <- "_DELrcnt"
  plink.rmv3 <- cat.path(dir="",plink.rmv2,suf=suf.for.rmv,ext="cnv")
  filter.hifreq.cnv.samples(delR.file, nsamps, rem.file=plink.rmv2, rare="DEL", dir=dir, 
                            thr.sd=thr.sd, ndels=NULL, ndups=NULL, 
                            del.rate=del.rate, dup.rate=dup.rate, QC.ON=rare.qc, pval=pval)
  suf.for.rmv <- "_PLATEcnt"
  plink.rmv4 <- cat.path(dir="",plink.rmv3,suf=suf.for.rmv,ext="cnv")
  filter.hifreq.plates(plink.rmv3,sample.info=sample.info,rem.file=plink.rmv3,dir=dir,rare=F,
                       psd=plate.thr,QC.ON=plate.qc,rmv.low.plates=rmv.low.plates,suf.for.rmv=suf.for.rmv)

  suf.for.rmv <- "_PLATErcnt"
  plink.rmv5 <- cat.path(dir="",plink.rmv4,suf=suf.for.rmv,ext="cnv")
  filter.hifreq.plates(c(delR.file,dupR.file),sample.info=sample.info,rem.file=plink.rmv4,dir=dir,rare=T,
                       psd=plate.thr,QC.ON=plate.qc,rmv.low.plates=rmv.low.plates,suf.for.rmv=suf.for.rmv)
  
  #f1 <- cat.path(dir="",rmv.ext(plink.out),suf="allFilt",ext="cnv"); g1 <- plink.rmv5
  #f2 <- cat.path(dir="",rmv.ext(plink.out),suf="allFilt",ext="cnv.map"); g2 <- paste(plink.rmv5,".map",sep="")
  #f3 <- cat.path(dir="",rmv.ext(plink.out),suf="allFilt",ext="fam"); g3 <- cat.path(dir="",rmv.ext(plink.rmv5),ext="fam")
  #system(paste("cp",g1,f1),intern=hide.plink.out)
  #system(paste("cp",g2,f2),intern=hide.plink.out)
  #system(paste("cp",g3,f3),intern=hide.plink.out)
  
  commands[15] <- "# remove samples with too many rare or total CNVs per plate in R"
  # all cnvs
  commands[16] <- paste("plink --cfile",rmv.ext(plink.rmv5),"--cnv-drop-no-segment --cnv-write --out allCNV --allow-no-sex")
  # make map file in plink
  commands[17] <- paste("plink --cnv-list allCNV.cnv --cnv-make-map --out allCNV")
  # all dels
  commands[18] <- paste("plink --cfile allCNV --cnv-del --cnv-drop-no-segment --cnv-write --out allDel --allow-no-sex")
  # all dups
  commands[19] <- paste("plink --cfile allCNV --cnv-dup --cnv-drop-no-segment --cnv-write --out allDup --allow-no-sex")
  # rare dups
  commands[20] <- paste("plink --cfile allCNV --cnv-dup --cnv-sites",min.sites,"--cnv-drop-no-segment --cnv-freq-exclude-above",rare.cut,"--cnv-overlap",rare.olp,"--cnv-write --out",dupR.file,"--allow-no-sex")
  # make a map file for each
  commands[21] <- paste("plink --cnv-list ",paste(delR.file,".cnv",sep="")," --cnv-make-map --out ",delR.file,sep="")
  commands[22] <- paste("plink --cnv-list ",paste(dupR.file,".cnv",sep="")," --cnv-make-map --out ",dupR.file,sep="")
  #nn[[5]] <- sapply(commands[16:19],system,intern=hide.plink.out,USE.NAMES=F)
  #nn[[6]] <- system(commands[20],intern=hide.plink.out)
  #nn[[7]] <- sapply(commands[21:22],system,intern=hide.plink.out,USE.NAMES=F)
  # rare dels
  #print(min.sites); print(rare.cut); print(rare.olp); print(delR.file);
  commands[23] <- paste("plink --cfile allCNV --cnv-del --cnv-sites",min.sites[1],"--cnv-drop-no-segment --cnv-freq-exclude-above",rare.cut[1],"--cnv-overlap",rare.olp[1],"--cnv-write --out",delR.file[1],"--allow-no-sex")
  commands[24] <- paste("wc -l *.cnv | grep -vi '.map.cnv' | grep -vi 'total' | sort -bgr > cnv_qc_summary.txt",sep="")  # count the number of cnvs at various stages

  ## JUST FOR RARE ##
  # create segment files (in each 1mb count cases vs controls, make pretty pic) # use this to graph
  commands[25] <- paste("plink --cfile",delR.file,"--cnv-seglist  --out",delR.file,"--allow-no-sex")
  commands[26] <- paste("plink --cfile",dupR.file,"--cnv-seglist  --out",dupR.file,"--allow-no-sex")
  # create CNVR regions file(s)
  commands[27] <- paste("plink --cfile",delR.file,"--segment-group --out ",delR.file," --cnv-overlap 0.5 --allow-no-sex")
  commands[28] <- paste("plink --cfile",dupR.file,"--segment-group --out ",dupR.file," --cnv-overlap 0.5 --allow-no-sex")
  
  nn[[5]] <- suppressWarnings(sapply(commands[16:28],system,intern=hide.plink.out,USE.NAMES=F))
  #nn[[9]] <- system(commands[28],intern=hide.plink.out)
  #nn[[6]] <- sapply(commands[28],system,intern=hide.plink.out,USE.NAMES=F)
  #now nn contains all suppressed plink output
  sample.info <- samp.update.qc.fail(sample.info,dir=dir,verbose=F,proc=6)
  write.sample.info(sample.info,dir)
  return(commands)
}



q.submit <- function(penncmd,cc,dest,cluster.fn="q.cmd",grid.name="all.q",
                     filepref="qsub_done_",sub.dir.nm="complete_flags/") {
  # make the sh file
  # make the qsub call
  # submit to queue
  sub.dir.path <- paste(dir.force.slash(dest),sub.dir.nm,sep="")
  if(!file.exists(sub.dir.path)) {
    dir.create(paste(sub.dir.path),showWarnings=F)
  }
  extraline <- paste("touch ",sub.dir.path,filepref,cc,sep="")
  nxt.fn <- cat.path(dest,"penn",suf=cc,ext="sh")
  next.content <- paste(c(penncmd,extraline))
  writeLines(next.content,con=nxt.fn)
  #q.call <- paste("qsub -q ",grid.name," -o ",dest," -j y ",nxt.fn,sep="")
  if(is.function(cluster.fn) | (is.character(cluster.fn))) {
    q.call <- do.call(cluster.fn,args=list(file.name=nxt.fn,output.dir=dest,id=grid.name))
  }
  call.out <- system(q.call,intern=T)
  return(call.out)
}



q.test.done <- function(dest,filepref="qsub_done_",sub.dir.nm="complete_flags/") {
  # check the flag files and return the numbers of any completed
  sub.dir.path <- paste(dir.force.slash(dest),sub.dir.nm,sep="")
  if(file.exists(sub.dir.path)) {
    indir <- list.files(sub.dir.path)
    if(length(indir)>0) {
      present <- as.numeric(gsub(filepref,"",indir))
      if(length(narm(present))>0) {
        return(narm(present))
      }
    } 
  } else {
    warning("qsub path doesn't exist!")
  }
  return(numeric(0))
}


do.penn.qsub <- function(penn.calls, dir, hrs.guess=1, grid.name="all.q", cluster.fn="q.cmd") {
  cat("\nSubmitting PennCNV calls to grid...\n")
  cat("\n expect HMM processing via the cluster to take roughly ",round(hrs.guess,2),"hrs ...\n",sep="")
  cat(" note that if restarting previous run, existing penn-cnv jobs should be deleted to avoid conflict\n")
  ## remove existing files
  #print("can remove this temporary assignment for next full run")
  #dir[["cnv.qsub"]] <- paste(dir$cnv,"PENNCLUSTER/",sep="")
  #print(dir$cnv.qsub)
  if(length(list.files(dir$cnv.qsub))>0) {
    delete.file.list(cat.path(dir$cnv.qsub,list.files(dir$cnv.qsub)),verbose=F) # delete existing files in dir
    warning("deleted previous contents of directory: ",dir$cnv.qsub)
  }
  
  for (tt in 1:length(penn.calls)) {
    # make the sh file
    # make the qsub call
    # q.submit(penncmd,cc,destination,dir)
    q.submit(penncmd=penn.calls[tt],tt,dest=dir$cnv.qsub,cluster.fn=cluster.fn,grid.name=grid.name) 
  }
  cat(" submitted",length(penn.calls),"jobs to",grid.name,"[qsub]\n")
  # check whether it's actually running
  wait(45,"s"); files.so.far <- list.files(dir$cnv.qsub)
  if(length(files.so.far)==0) {
    warning("qsub command seems to have failed. You have 30 seconds to cancel and fix the issue, or\n",
    "the PennCNV commands will be run serially (can be slow)")
    wait(30,"s")
    return(F)
  } else {
    not.complete <- T
    kk <- proc.time()[3]
    done.now <- 0; no.act <- 0; notes.per.hr <- 2; hrs.so.far <- 0
    while(not.complete) {
      # check whether it's finished and update progress occassionally
      cat("."); wait(60,"s") # only bother checking every 60 seconds (show a dot per minute)
      which.done <- q.test.done(dest=dir$cnv.qsub)
      if(length(which.done)>length(done.now)) {
        cat(" completed ",length(which.done),"/",length(penn.calls)," qsub jobs\n",sep="")
        done.now <- which.done; no.act <- 0
        if(length(done.now)>=length(penn.calls)) {
          not.complete <- F
        }
      } else {
        tmp <- cat.path(dir$cnv.qsub,"temp.temp")
        system(paste("qstat | cat >",tmp)) # write queue status update to file
        temp <- readLines(tmp)
        unlink(tmp)
        if(length(temp)<1) { no.act <- no.act + 1 } 
        if(no.act>10) { cat(" several minutes of no qsub activity detected - suggest cancelling [ctrl-c]\n") }
      }
      jj <- proc.time()[3]
      latest.hrs <- ((jj-kk)/3600); 
      if(floor(latest.hrs*notes.per.hr)>(hrs.so.far*notes.per.hr)) { 
        hrs.so.far <- floor(latest.hrs*notes.per.hr)/notes.per.hr
        cat("PennCNV qsub now has been running for",hrs.so.far,"hours\n")
      }
    }
    jj <- proc.time()[3]
    cat(" PennCNV qsub commands processing completed in",round(((jj-kk)/3600),2),"hours\n")
  }
  return(T)
}


run.PENN.cnv <- function(DT=NULL,dir=NULL,num.pcs=NA,LRR.fn,BAF.fn="BAFdescrFile",
            n.cores=3,q.cores=NA,grid.id="all.q",cluster.fn="q.cmd",relative=T,run.manual=F,low.ram=T,
            penn.path="/usr/local/bin/penncnv/",ucsc="hg18",sample.info=NULL,snp.info=NULL,
            restore.mode=F,print.cmds=T,hide.penn.out=T,hmm="hh550.hmm",use.penn.gc=T)
{
  #takes PC corrected data, generates the prerequisite PENN-CNV input files
  # and runs PennCNV to detect CNVs (or gives the commands to run manually)
  # will import and process the BAF data/file during this process
  load.all.libs() # load all main plumbcnv libraries
  if(!is.character(grid.id)) { grid.id <- NA }
  if(!check.linux.install("qsub")) { q.cores <- NA }
  if(!is.numeric(q.cores[1]) | (as.numeric(q.cores[1])<1) | is.na(grid.id[1]) ) { q.cores <- NA }
  if(n.cores>1) { multi <- T } 
  marker.fn <- NULL
  if(!is.na(q.cores)) { qsub <- T ; multi <- F ; relative <- F } else { qsub <- F } # run penn on qsub
  if(is.data.tracker(DT)) {
    DTs.required <- c("ngrps","sample.info","snp.info","big.pcc","big.baf")
    is.valid <- req.valid(DT,DTs.required,n.pcs=num.pcs)
    if(!is.valid) { 
      ## if invalid let it try to work via loading in separate function parameters
      warning("file required for 'run.PENN.cnv' invalid in data.tracker")
      if(!req.valid(DT,"big.pcc",n.pcs=num.pcs)) { 
        cat("NB: if more than one set of PCs have been calculated, must specify num.pcs\n") 
      }
    }
    if(is.null(dir)) {
       dir <- getSlot(DT,"dir")
    }
    LRR.fn <- getSlot(DT,"big.pcc",n.pcs=num.pcs)
    BAF.fn <- getSlot(DT,"big.baf")
  }
  cat("\nPreparing for the running of PennCNV.\n")
  if(!file.exists(penn.path)) {
    cat("Error: Make sure this open-source perl script is correctly installed to:",penn.path,"\n")
    cat("[NB: Perl must also be installed, which should be fairly standard on linux]\n")
    stop("Find PennCNV at: http://www.openbioinformatics.org/penncnv/penncnv_download.html")
  }
  dir <- validate.dir.for(dir,c("cnv","cnv.raw","cnv.pen","cnv.qsub"),warn=F)
  if(is(snp.info)[1]!="RangedData") {  snp.info <- read.snp.info(dir) }
  #if(!is.data.frame(sample.info)) {  sample.info <- read.sample.info(dir)  }
  ##sample.info <- validate.samp.info(sample.info,QC.update=T,verbose=F) only used by process.penn.results 
  # which validates it anyway
  if(is.null(LRR.fn)) { LRR.fn <- "" }
  if(LRR.fn=="") {  LRR.fn <- confirm.corrected.filename(dir,num.pcs) }
  ndir <- count.penn.dirs(dir) 
  ndir.test <- (ndir=={ if(qsub) q.cores else n.cores })
  if(!restore.mode | !ndir.test) {
    ## must delete any existing files in penn raw subdirectories! ##
    if(length(ndir)>0 & length(list.files(dir$cnv.raw))>0) { delete.file.list(cat.path(dir$cnv.raw,list.files(dir$cnv.raw)),verbose=F);
                              cat(" deleting old penn-cnv raw subject files from:\n",dir$cnv.raw,"\n") }
    if(qsub) { cores <- q.cores } else { cores <- n.cores }
    dir.contents <- prepare.penncnv.samples(dir,LRR.fn=LRR.fn,BAF.fn=BAF.fn,num.pcs=num.pcs,
                                       low.ram=low.ram,ucsc=ucsc,n.cores=n.cores,cores=cores)
    lenz <- unlist(sapply(dir.contents,length))
    if(length(lenz)>10) {
      sz.txt <- paste("N=",names(table(lenz))," (",table(lenz),")",sep="",collapse="; ")
    } else { sz.txt <- paste(lenz,collapse=",") }
    cat(" number of samples written per PENN-CNV directory:",sz.txt,"\n")
  } else {
    # assume there are already files in the penn directories #TEMP for debug
    dir.contents <- penn.name.files(dir,write.files=F,ret.fn=F,relative=relative) 
    samp.fnd <- length(unlist(dir.contents))
    if(samp.fnd<10) 
      { stop("Error: option to use existing PennCNV files didn't find enough sample files") }
    ## restoring from existing files so assume we already have a .gcm and .pfb file
    cat(" detected existing set of ",samp.fnd," valid penn-cnv raw files (",ndir," subsets)\n",sep="")
    gc.out.fn <- list.files(dir$cnv,pattern=".gcm")
    baf.out.fn <- list.files(dir$cnv,pattern=".pfb")
    if(length(gc.out.fn)>0 & length(baf.out.fn)>0) {
      marker.fn <- list(cat.path(dir$cnv,gc.out.fn[1]),cat.path(dir$cnv,baf.out.fn[1]))
      names(marker.fn) <- c("gc","baf")
    }
  }
  if(is.null(marker.fn)) { marker.fn <- prepare.penncnv.markers(snp.info,
                                  lrr.descr=LRR.fn,baf.descr=BAF.fn,dir=dir,ucsc=ucsc,n.cores=n.cores) }
  #print(marker.fn)
  ##  need to make sure that any existing penn-cnv output files are deleted
  ext.files <- cat.path(dir$cnv.pen,list.files(dir$cnv.pen)); 
  if(length(list.files(dir$cnv.pen))>0) { delete.file.list(ext.files,verbose=F);
                cat(" deleting old penn-cnv output files from:\n",dir$cnv.pen,"\n") }
  
  if(!run.manual) { 
    cat("\nRunning the PennCNV HMM - expect this to take several hours\n")
    cat(" NB: to run with family data, press Ctrl-C and run PennCNV manually with appropriate options\n")
  }
  if(!multi & !run.manual & !qsub) cat(" may also be faster to run PennCNV manually on several computers in parallel\n")
  ndir <- count.penn.dirs(dir) 
  penn.calls <- get.penn.cmd(dir,ndir=ndir,gc.out.fn=marker.fn$gc,baf.out.fn=marker.fn$baf,
                             relative=relative,penn.path=penn.path,hmm=hmm,use.penn.gc=use.penn.gc)
  n.calls <- length(penn.calls)
  if(multi) {  must.use.package("multicore",F) } # run in parallel
  if(relative) { 
    rel.cd <- paste("cd",dir$cnv) 
    cur.dir <- getwd(); system(rel.cd)  # use relative paths
  }
  ## count files against number of parallel processes to estimate proc time ##
  tot.to.proc <- sum(sapply(dir.contents,length))
  if(multi) { tot.to.proc <- (tot.to.proc/(n.calls*.75)) } # faster with parallel proc
  if(qsub) { tot.to.proc <- (tot.to.proc/(n.calls*.9)) }
  per.hr <- 250; hrs.guess <- tot.to.proc/per.hr
  # create multi parallel dir structure and name files for penn
  subj.in.each <- penn.name.files(dir,write.files=T,ret.fn=F,relative=relative)
  if(print.cmds) { cat("\nCalls sent to PennCNV:\n\n"); print(penn.calls) }
  if(run.manual) { 
    cat("\n\nCOPY AND PASTE COMMANDS BELOW TO RUN MANUALLY IN TERMINAL\n\n")
    if(relative) { cat(rel.cd,"\n\n") }
    cat(paste(penn.calls,"\n\n"),sep="") ; cat("\n\n")
    return(penn.calls) 
  }
 
  ## RUN PENN CNV ##
  if(relative) {  setwd(dir$cnv) } # only works in 'relative mode' when calls are made from here
  if(!multi) {
    if(qsub){
      qsub <- do.penn.qsub(penn.calls=penn.calls, dir=dir, hrs.guess=hrs.guess, 
                           grid.name=grid.id, cluster.fn=cluster.fn)
    } 
    if(!qsub) {
      # run one after the other
      options(warn = -1)
      time.per.it <- (hrs.guess/n.calls)*60; total.time <- 0
      for (tt in 1:n.calls) {
        cat(" expect HMM process ",tt,"/",n.calls," to take roughly ",round(time.per.it,2),"minutes\n",sep="")
        kk <- proc.time()
        killme <- system(penn.calls[tt],intern=hide.penn.out, ignore.stderr=hide.penn.out)
        jj <- proc.time()
        time.per.it <- round((jj[3]-kk[3])/60)
        total.time <- total.time+time.per.it
        cat(" PennCNV processing for the first ",tt," calls has taken",time.per.it,"minutes so far\n",sep="")
      }
      options(warn = 0)
    }
  } else {
    kk <- proc.time()
    nullList <- list()
    cat(paste("\n parallel PennCNV processing for",n.calls,"files initiated: "))
    for (tt in 1:n.calls) {
      cat(paste(tt,"..",sep=""))
      nullList[[tt]] <- multicore::parallel(system(penn.calls[tt],intern=hide.penn.out, ignore.stderr=hide.penn.out))
    }  
    cat("\n expect HMM processing to take roughly ",round(hrs.guess,2),"hrs ...",sep="")
   # print(sapply(nullList,is))
   # print(penn.calls[tt])
    nullList <- collect(nullList,wait=T)
    jj <- proc.time()
    cat(" parallel PennCNV took",round((jj[3]-kk[3])/60),"minutes\n")
  }
  if(relative) { setwd(cur.dir) } # go back to original directory
  # run CNV-QC on processed files
  if(is.data.tracker(DT)) {
    DT <- setSlot(DT,proc.done=5,gc.marker=cat.path(dir$gc,"marker.gc.RData"),warns=F)
    return(DT)
  } else {
    return(penn.calls)
  }
}  
  

run.CNV.qc <- function(DT=NULL,dir=NULL,num.pcs=NA,penn.path="/usr/local/bin/penncnv/",
                         ucsc="hg18",sample.info=NULL,snp.info=NULL,
                         out.format="Ranges",result.pref="cnvResults",
                         cnv.qc=T,rare.qc=T,plate.qc=T,restore.mode=F,
                         pval=0.05,del.rate=0.4,dup.rate=0.18,thr.sd=3,plate.thr=3,rmv.low.plates=F,
                         min.sites=10,rare.olp=.5,rare.pc=1,rmv.bad.reg=T,hide.plink.out=T,verbose=F)
{
  #takes PC corrected data, generates the prerequisite PENN-CNV input files
  # and runs PennCNV to detect CNVs (or gives the commands to run manually)
  # will import and process the BAF data/file during this process
  load.all.libs(more.bio=c("BSgenome","Rsamtools","rtracklayer","biomaRt","gage",
                           "graph","multtest","GenomicFeatures","AnnotationDbi")) # load all main plumbcnv libraries
  if(is.data.tracker(DT)) {
    DTs.required <- c("ngrps","sample.info","snp.info")
    is.valid <- req.valid(DT,DTs.required,n.pcs=num.pcs)
    if(!is.valid) { 
      ## if invalid let it try to work via loading in separate function parameters
      warning("file required for 'run.CNV.qc' invalid in data.tracker")
      if(!req.valid(DT,"big.pcc",n.pcs=num.pcs)) { 
        cat("NB: if more than one set of PCs have been calculated, must specify num.pcs\n") 
      }
    } else {
      if(is.null(dir)) {
        dir <- getSlot(DT,"dir")
      }
    }
  }
  dir <- validate.dir.for(dir,c("cnv","cnv.qc","gc","res"))
  cat("Running CNV-QC on processed PennCNV files\n")
  if(!restore.mode) {
    warning("deleting existing plink files from directory:",dir$cnv.qc)
    if(length(list.files(dir$cnv.qc))>0) { delete.file.list(cat.path(dir$cnv.qc,list.files(dir$cnv.qc)),verbose=F) }
  } else {
    warning("restore mode selected, any existing CNVQC files won't be deleted",
            " before this function proceeds; results may be inaccurate")
  }
  my.out <- process.penn.results(dir,penn.path=penn.path,sample.info=sample.info,
                   ucsc=ucsc,hide.plink.out=hide.plink.out,
                   cnv.qc=cnv.qc,rare.qc=rare.qc,plate.qc=plate.qc,pval=pval,
                   del.rate=del.rate,dup.rate=dup.rate,thr.sd=thr.sd,plate.thr=plate.thr,
                   min.sites=min.sites,rare.olp=rare.olp,rare.pc=rare.pc,rmv.low.plates=rmv.low.plates,
                   rmv.bad.reg=rmv.bad.reg,verbose=verbose) 
  cat("\nCNV-QC complete\n")
  penn.cmds <- my.out
  if(!out.format %in% c("cnvGSA","Ranges")) { out.format <- "Ranges" }
  cnvResults <- import.called.cnvs(dir=dir,out.format=out.format)
  cnv.res.fn <- (cat.path(dir$res,result.pref,ext="RData",suf=num.pcs))
  save(cnvResults,file=cnv.res.fn)
  cat("~wrote CNVs as",out.format,"object, saved to:\n ",cnv.res.fn,"\n")
  if(is.data.tracker(DT)) {
    DT <- setSlot(DT,cnv.result=cnv.res.fn,n.pcs=num.pcs,proc.done=6,
                              gc.marker=cat.path(dir$gc,"marker.gc.RData"),warns=F)
    return(DT)
  } else {
    return(cnvResults)
  }
}


extract.val.penn <- function(penn.col)
{
  # penn cnv files have format varname=var.value; extract just the value
  split.list <- strsplit(penn.col,"=",fixed=T)
  copies <- sapply(split.list,tail,1)
  return(copies)
}


count.penn.dirs <- function(dir) {
  # looks at list of files in dir$cnv.raw and finds the highest suffix p-xx
  # and assumes this corresponds to the number of penn-cnv directories
  dir <- validate.dir.for(dir,"cnv.raw")
  allf <- list.files(dir$cnv.raw)
  splitz <- strsplit(allf,"p",fixed=T)
  nz <- sapply(splitz,tail,1)
  suppressWarnings(n_z <- narm(as.numeric(paste(nz))))
  if(length(n_z)<1) { warning(paste("no pXX directories detected in:\n ",dir$cnv.raw)); return(0) }
  highest <- max(n_z)
  all.there <- (1:highest %in% n_z)
  if(length(which(all.there))<highest) {
    warning(paste("some pXX directories missing:",paste(which(!all.there),collapse=",")))
  }
  return(highest)
}


remove.for.step <- function(DT,step=0,n.pcs=NA) {
  ## helps reduce file space by deleting large files no longer
  ## needed from previous steps of the pipeline. can delete after complete step 'n'
  if(!all(step %in% c(0:6))) {
    warning("step(s) must be entered as numbers from 0:6. No files removed")
    return(F)
  }
  if(!is.data.tracker(DT)) {
    warning("DT must be a data.tracker object. No files removed")
    return(F)
  }
  steps.complete <- getSlot(DT,"proc.done") # see which steps have been done in DT
  if((1 %in% step) & any(steps.complete>0)){
    cat("Deleting the long format text files produced in import from Genome Studio file\n")
    del.files <- c(getSlot(DT,"shape.lrr"),getSlot(DT,"shape.baf"))
    delete.file.list(del.files)
  }
  if(2 %in% step & any(steps.complete>1)){
    cat("Deleting snpMatrix files produced in snp-QC\n")
    # just the snpMatrix files
    del.files <- list.files(dir$cr,pattern=".RData")
    delete.file.list(cat.path(dir$cr,del.files))
  }
  if(3 %in% step & any(steps.complete>2)){
    cat("Deleting big.matrix files produced in import to big.matrix\n")
    first.bigs <- unlist(getSlot(DT,"big.lrr"))
    for (cc in 1:length(first.bigs)) { delete.big.matrix.files(first.bigs[cc],dir) }
    cat("Deleting big.matrix files produced in sample-QC\n")
    first.bigs <- c(unlist(getSlot(DT,"big.filt")),unlist(getSlot(DT,"pcas")))
    for (cc in 1:length(first.bigs)) { delete.big.matrix.files(first.bigs[cc],dir) }
  }
  if(4 %in% step & any(steps.complete>3)){
    cat("Deleting big.matrix files produced for PCA, and combined sample-QC file\n")
    bigfls <- list.files(dir$big)
    find.sub.mat <- grep("PCAMatrixSort",bigfls) # try to find pca sub matrix big files
    if(length(find.sub.mat)>0) { unlink(cat.path(dir$big,bigfls[find.sub.mat])) }
    first.bigs <- c(unlist(getSlot(DT,"big.qc")))
    for (cc in 1:length(first.bigs)) { delete.big.matrix.files(first.bigs[cc],dir) }
  }
  if(5 %in% step & any(steps.complete>4)){
    cat("Deleting PC-corrected big matrix file\n")
    first.bigs <- unlist(getSlot(DT,"big.pcc",n.pcs=n.pcs))
    for (cc in 1:length(first.bigs)) { delete.big.matrix.files(first.bigs[cc],dir) }
    cat("Deleting Penn CNV raw subject LRR/BAF text files and directories\n")
    del.list <- cat.path(dir$cnv,list.files(dir$cnv.raw))
    if(length(list.files(dir$cnv.raw))>0) { delete.file.list(del.list) }
  }
  if(6 %in% step & any(steps.complete>5)){
    cat("No files yet set for deletion after set 6 (files are all small)\n")
  }
  return(T)
}


delete.big.matrix.files <- function(bigMat,dir) {
  # delete the backing, description and optionally 'RData' file for a specific
  # bigMatrix. RData only deleted if bigMat is the name of this file
  dir <- validate.dir.for(dir,"big")
  to.kill <- NULL
  if(is.character(bigMat)) {
    if(is.file(bigMat,dir$big)) {
      if(length(grep("RData",bigMat))>0) {
          to.kill <- bigMat
      }
    } else {
      return(NULL) # was like a filename but couldn't find the file
    }
  }
  bm <- getBigMat(bigMat,dir$big)
  ott <- describe(bm)
  bckFileName <- attr(ott,"description")$filename
  rm(bm); rm(ott); gc() # remove the objects in case deleting the backing files causes problems
  to.kill <- c(to.kill,cat.path(dir$big,bckFileName))
  all.files.dir <- cat.path(dir$big,list.files(dir$big))
  # this bit is a hack - read the first 3 lines of all files in directory to
  # see whether they are description files, then later whether they point to the
  # big matrix backing we want to delete
  all.files.hds <- suppressWarnings(lapply(all.files.dir,function(x) { readLines(x,n=3) }))
  bigDescs <- suppressWarnings(sapply(all.files.hds,
                    function(x) { length(grep("big.matrix.descriptor",paste(x)))>0 }))
  if(length(which(bigDescs))>0) {
    all.descr.files <- all.files.dir[bigDescs]
    all.descr.heads <- all.files.hds[bigDescs]
    for (cc in 1:length(all.descr.files)) {
      if(length(grep(paste(bckFileName),all.descr.heads[cc]))>0) {
        to.kill <- c(to.kill,all.descr.files[cc])
      }
    }
  }
  delete.file.list(cat.path(dir$big,to.kill))
  return(to.kill)
}


delete.file.list <- function(files,verbose=T) {
  ## delete a list of files
  if(is.ch(files))  { files <- unlist(files) } else { cat("invalid file list"); return(NULL) }
  if(length(files)<1) { return(NULL) }
  for (cc in 1:length(files)) {
    if(file.exists(files[cc])) {
      if(verbose) { cat(" deleting:",files[cc],"\n") }
      unlink(files[cc],recursive=T)
    } else {
      cat(" couldn't find file:",files[cc],"\n")
    }
  }
}


delete.all.files <- function(dir,rescue=NULL,to="~/Documents/") {
  if(!is.list(dir)) { stop("'dir' should be a list of directories specific to plumbCNV")}
  if(length(rescue)>0 & is.ch(rescue)) {
    if(!file.exists(dirname(to))) { stop("stopped delete as rescue path did not exist") }
    rescue <- unlist(rescue)
    for (cc in 1:length(rescue)) {
      if(is.file(rescue[cc],dir)) {
        dest <- cat.path(to,basename(rescue[cc]))
        cat("\n~saving a copy of:\n ",rescue[cc],"to:\n ",dest,"\n")
        file.copy(from=rescue[cc], to=dest, overwrite = T,
                  recursive = T,copy.mode = T)
      }
    }
  } else {
    if(!is.null(rescue)) { stop("stopped delete as a non-character 'rescue' list was entered")}
  }
  cat("\n\n====== DELETING ALL FILES IN:",dir$base,"=======\n\n")
  init.dirs.fn(dir,overwrite=T,silent=F,update.bash=F,plate.info=NULL,pheno=NULL,file.spec=NULL)
}


print.big.matrix <- function(bigMat,dir="",name=NULL,dat=T,descr=NULL,bck=NULL,mem=F,row=3,col=2,
                          rcap="SNPs",ccap="Samples",rlab="SNP-id",clab="Sample IDs",...) {
  print.big.matrix(bigMat=bigMat,dir=dir,name=name,dat=dat,descr=descr,bck=bck,
                    mem=mem,row=row,col=col,
                    rcap=rcap,ccap=ccap,rlab=rlab,clab=clab)
}


#dummylst <- Rfile.index("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")

#rm(dummylst)


####### QUALITY SCORE FUNCTIONS ####

do.roc <- function(y=numeric(),x=numeric(),fn="ROC.pdf",plot=F) {
  myform <- as.formula("y ~ x")
  dd <- predict(glm(myform,family=binomial("logit")),type="response")
  if(plot) { pdf(cat.path(dir$res,fn)) }
  X <- ROC(test=dd,stat=y)
  if(plot) { dev.off() }
  mx <- max(X$res[, 1] + X$res[, 2])
  mhv <- which((X$res[, 1] + X$res[, 2]) == mx)
  fv <- X$res[,5]
  mxf <- fv[mhv]
  cat(paste("Sens: ", formatC(100 * X$res[mhv, 1], digits = 1, format = "f"), 
            "%\n", "Spec: ", formatC(100 * X$res[mhv, 2], digits = 1, format = "f")),"\n")
  return(c(mxf=mxf))
}


qscore.cnvs <- function(cnv.ranges,DT=NULL,file="all.ranges.pdf",dir="",pc.flank=5,snp.info=NULL,
                        col1="black",scheme="mono",LRR=T,BAF=F,bafOverlay=F,hzOverlay=F,
                        PREPOSTPC=F,baf.file="big.baf",lrr.file="big.pcc",n.cores=1,...) {
  must.use.package(c("multicore","bigmemory","biganalytics")); must.use.package("genoset",T)
  if(is(cnv.ranges)[1]!="RangedData") { warning("not a RangedData object"); return(NULL) }
  if(any(!c("id") %in% colnames(cnv.ranges))) { warning("cnv.ranges must contain id"); return(NULL) }
  idz <- cnv.ranges$id; stz <- start(cnv.ranges); enz <- end(cnv.ranges); chrz <- chr(cnv.ranges); wz <- width(cnv.ranges)
  dir <- validate.dir.for(dir,c("big","res"))
  if(!is.data.tracker(DT) & all(dir!="")) { DT <- read.data.tracker(dir,warn.only=T) }
  if((!is.character(lrr.file) | !is.character(baf.file))) { lrr.file <- baf.file <- "";
                                                            warning("lrr.file/baf.file must be character locations, defaults will be used from tracker") }
  XX <- get.tracker.for.plot(DT=DT,dir=dir,PREPOSTPC=PREPOSTPC,baf.file=baf.file,lrr.file=lrr.file,
                             LRR=LRR,BAF=BAF,samples=idz)
  if(all(c("lrr.file","baf.file","lrr.file.raw") %in% names(XX))) {
    lrr.file <- XX$lrr.file; baf.file <- XX$baf.file; lrr.file.raw <- XX$lrr.file.raw
    if(!PREPOSTPC) {
      # load big objects to save loading each time for many ranges
      if(LRR) { lrr.file <- getBigMat(lrr.file,dir$big) }
      if(BAF) { baf.file <- getBigMat(baf.file,dir$big) }
    } else {
      warning("PREPOSTPC==T is slow for multiple ranges as data must be reloaded each iteration")
    }
  } else {
    warning("lookup of baf.file and/or lrr.file file(s) failed - will try to proceed but failure is likely")
  }
  ofn <- cat.path(dir$res,file)
  n.to.plot <- nrow(cnv.ranges)
  if(!is(snp.info)[1]=="RangedData") { snp.info <- read.snp.info(dir) }
  if(!is.na(scheme) & ("color" %in% colnames(snp.info))) { snp.info <- snp.info[,-match("color",colnames(snp.info))] } 
  if(!all(c("QCfail","gindx") %in% colnames(snp.info))) {
    snp.info <- add.gindx.to.Ranges(snp.info)
    snp.info <- snp.update.qc.fail(snp.info,dir=dir)
  }
  cL <- get.chr.lens(dir)
  if(n.to.plot<1) { return(NULL) } else { cat(n.to.plot,"CNVs to plot:\n")}
  #pdf(ofn)
  if(n.to.plot>0) {
    pdf(ofn)
    for(cc in 1:n.to.plot) {
      poz <- c( max(1,(stz[cc]-(pc.flank*wz[cc]))), min((enz[cc]+(pc.flank*wz[cc])),cL[as.numeric(chrz[cc])]) ) #; print(poz)
      #  suppressWarnings({
      #   cnv.plot(DT=DT,samples=idz[cc],Chr=chrz[cc],Pos=poz,Cnv=c(stz[cc],enz[cc]),
      #           dir=dir,snp.info=snp.info,
      #          lrr.file=lrr.file,baf.file=baf.file,PREPOSTPC=PREPOSTPC,
      #         LRR=LRR,BAF=BAF
      
      Pos <- force.chr.pos(Pos,Chr,snp.info)
      
      }
    dev.off()
    }
  #dev.off()
  cat("~wrote file:",ofn)
  }


ratio.flanks <- function(pos,chr,ratio=5,bp=NA,chr.lens=get.chr.lens()) {
  stz <- pos[1]; enz <- pos[2]; wz <- enz-stz+1
  if(is.na(bp)) { widen.by <- ratio*wz } else { widen.by <- bp }
  poz1 <- c(max(1,(stz-(widen.by))), pos[1]-1)
  poz2 <- c(pos[2]+1,min((enz+(widen.by)),chr.lens[as.numeric(chr)]) )
  return(c(poz1,poz2))
}


get.ratio.set <- function(ranged,ratio=5,bp=NA) {
  stz <- start(ranged)    ; enz <- end(ranged); chr <- chr(ranged)
  chr.lens <- get.chr.lens()
  tt <- matrix(nrow=length(stz),ncol=4)
  for (cc in 1:length(stz)) { tt[cc,] <- ratio.flanks(c(stz[cc],enz[cc]),chr[cc],ratio=ratio,bp=bp,chr.lens=chr.lens) }
  return(tt)
}


big.extract.snp.ranges <- function(first.last.snps,samples,bigMat,snp.info=NULL) {
  if(is.null(dim(first.last.snps))) { st <- first.last.snps[1]; en <- first.last.snps[2] 
  } else { st <- first.last.snps[,1]; en <- first.last.snps[,2] }
  if(is.null(snp.info)) {
    sorter <- 1:nrow(bigMat)
  } else {
    sorter <- match(rownames(snp.info),rownames(bigMat))
    if(any(is.na(sorter))) { warning("non matches between snp.info and bigMat data") }
    sorter <- narm(sorter)
  }
  idz <- samples; #print(length(idz)); print(dim(first.last.snps))
  rN <- rownames(bigMat)[sorter]; cN <- colnames(bigMat)
  st.rw <- match(st,rN); en.rw <- match(en,rN); nrng <- length(st.rw)
  samp.cl <- match(idz,cN)
  invalids <- which(is.na(st.rw) | is.na(en.rw) | is.na(samp.cl))
  out.list <- vector("list",nrng)
  for(cc in 1:nrng) {
    if(cc %in% invalids) {
      out.list[[cc]] <- NA 
    } else {
      out.list[[cc]] <- as.numeric(bigMat[sorter[st.rw[cc]:en.rw[cc]],samp.cl[cc]])
    }
  }
  #print(length(out.list))
  names(out.list) <- idz
  return(out.list)
}


prob.norm <- function(x,mu=0,sig=1) {
  X <- (x-mu)/(max(10^-100,abs(sig)))
  return(1-pnorm(abs(X)))
}


get.adj.nsnp <- function(snp.info,ranged,nsnp=10) {
  snp.info <- toGenomeOrder(snp.info); rw.cnt <- 1
  all.fl <- matrix(ncol=4, nrow=0)
  for(cc in chr.nums.of.info(ranged)) {
    nxt.nm <- rownames(snp.info[paste(cc)]); pos <- start(snp.info[paste(cc)])
    rng <- ranged[paste(cc)]
    st.en.snp <- range.snp(snp.info,ranged=rng)
    fl <- matrix(ncol=4, nrow=nrow(st.en.snp))
    fl[,2] <- start(rng); fl[,3] <- end(rng);
    for(dd in 1:nrow(st.en.snp)) {
      x1 <- pos[max(1,match(st.en.snp[dd,1],nxt.nm)-nsnp)]
      x2 <- pos[min(length(nxt.nm),match(st.en.snp[dd,2],nxt.nm)+nsnp)]
      #print(x1); print(x2); print(length(x1)); print(length(x2))
      fl[dd,1] <- x1[1]
      fl[dd,4] <- x2[1]
    }
    all.fl <- rbind(all.fl,fl)
  }
  fl <- (all.fl)
  fl[fl[,1]>fl[,2],1] <- fl[fl[,1]>fl[,2],2]
  fl[fl[,3]>fl[,4],1] <- fl[fl[,3]>fl[,4],4]
  return(fl)
}


get.flanks.from.big.mat <- function(ranged,bigMat,ratio=5,bp=NA,nsnp=NA,snp.info=NULL,L=T,R=T) {
  rownames(ranged) <- rN <- paste("CNV",1:nrow(ranged),sep="")
  if(!is.na(nsnp) & is(snp.info)[1]=="RangedData") {
    fl <- get.adj.nsnp(snp.info,ranged,nsnp)
    fl2 <- get.ratio.set(ranged,ratio=ratio,bp=bp)
    lenz <- (fl[,2]-fl[,1]) + (fl[,4]-fl[,3]); 
    lenz2 <- (fl2[,2]-fl2[,1]) + (fl2[,4]-fl2[,3]) ;
    print(length(which(lenz>lenz2)))
    fl3 <- fl; fl3[lenz>lenz2,] <- fl2[lenz>lenz2,]; fl <- fl3
  } else {
    fl <- get.ratio.set(ranged,ratio=ratio,bp=bp)
  }
  idz <- (ranged$id);  #print(length(idz)) #print(dim(ranged))
  if(!L & !R) { warning("L and R both unselected, returning nothing"); return(NULL) }
  if(L) {
    left <- RangedData(ranges=IRanges(start=fl[,1],end=fl[,2],names=rownames(ranged)),space=chr(ranged),id=idz)
    #rownames(left) <- rownames(ranged)
    flanking_1 <- range.snp(snp.info,toGenomeOrder(left,strict=T))
    if(nrow(flanking_1)!=length(idz)) { stop("unequal lengths, or id column missing from ranged") }
    if(!is.null(rownames(flanking_1))) { flanking_1 <- flanking_1[rN,] } # ensures same order
  }
  if(R) {
    right <- RangedData(ranges=IRanges(start=fl[,3],end=fl[,4],names=rownames(ranged)),space=chr(ranged),id=idz)
    #    rownames(right) <- rownames(ranged)
    flanking_2 <- range.snp(snp.info,toGenomeOrder(right,strict=T))
    if(nrow(flanking_2)!=length(idz)) { stop("unequal lengths, or id column missing from ranged") }
    if(!is.null(rownames(flanking_2))) { flanking_2 <- flanking_2[rN,] } # ensures same order
  }
  if(L & R) {
    flanking <- lapply(list(flanking_1,flanking_2),big.extract.snp.ranges,samples=idz,bigMat=bigMat,snp.info=snp.info)
    flanking <- mapply(FUN=c,flanking[[1]],flanking[[2]])
  } else {
    if(L) {
      #left
      flanking <- big.extract.snp.ranges(flanking_1,samples=idz,bigMat=bigMat,snp.info=snp.info)
    } else {
      #right
      flanking <- big.extract.snp.ranges(flanking_2,samples=idz,bigMat=bigMat,snp.info=snp.info)
    }
  }
  return(flanking)
}



kf.dist.test <- function(shrt,lng,k=10,prior=0.5) {
  itsz <- max(1,round(length(shrt)/k))
  minif <- function(shrt,itsz,reg.mn,reg.sd,prior=0.5) {
    excl <- sample(length(shrt),itsz)
    t.out <- suppressWarnings(most.likely.copy(lrr=shrt[-excl],
                              test=NA,prior=prior,reg.mn=reg.mn,reg.sd=reg.sd))
    return(c(t.out[1],t.out[2]))
  }
  reg.mn <- mean(lng,na.rm=T) ; reg.sd <- sd(lng,na.rm=T)
  state.ps <- replicate(k,minif(shrt,itsz=itsz,prior=prior,reg.mn=reg.mn,reg.sd=reg.sd))
  return(state.ps)
}


cnv.sanity.test <- function(shrt,lng,DEL=T,mag=1) {
  lS <- length(shrt); lL <- length(lng)
  if(DEL) { 
    P1a <- length(shrt[shrt < -.05])/lS
    P1b <- length(shrt[shrt < 0.01])/lS
    P2a <- length(lng[lng > -.25])/lL
    P2b <- length(lng[lng > -.5])/lL
    p <- P1a*(P1b^sqrt(lS))*sqrt(P2a*(P2b^sqrt(lL))) # sqrt(P2) to de-emphasise in case CNV not extended enough
    p <- P1a*P1b*sqrt(P2a*P2b) # sqrt(P2) to de-emphasise in case CNV not extended enough
  } else {
    P1a <- length(shrt[shrt > .01])/lS # allows more errors for dup
    P1b <- length(shrt[shrt > -.05])/lS # allows more errors for dup
    P2a <- length(lng[lng < .2])/lL
    P2b <- length(lng[lng < .4])/lL
    #p <- sqrt(P1a*(P1b^sqrt(lS))*sqrt(P2a*(P2b^sqrt(lL)))) # sqrt(as DUP mean is harder to discriminate, helps scores be more comparable)
    p <- sqrt(P1a*P1b*sqrt(P2a*P2b)) # sqrt(as DUP mean is harder to discriminate, helps scores be more comparable)
  }
  return(p)
}




t.kf.test <- function(shrt,lng,k=10) {
  itsz <- max(1,round(length(shrt)/k))
  minif <- function(shrt,itsz,lng) {
    excl <- sample(length(shrt),itsz)
    t.out <- t.test(x=shrt[-excl],y=lng,var.equal=F)
    return(c(t.out$statistic,p=t.out$p.value,t.out$estimate[1],t.out$estimate[2]))
  }
  ts <- replicate(k,minif(shrt,itsz=itsz,lng=lng))
  return(ts)
}

t.kf.test.fix <- function(cnv,fix,k=10,...) {
  itsz <- max(1,round(length(cnv)/k))
  minif <- function(cnv,itsz,fix) {
    excl <- sample(length(cnv),itsz)
    t.out <- t.test(x=cnv[-excl],mu=fix,...)
    return(c(t.out$statistic,p=t.out$p.value,t.out$estimate))
  }
  ts <- replicate(k,minif(cnv,itsz=itsz,fix=fix))
  return(ts)
}


test.lrr.as <- function(lrr, adj, copy=1, prior=0.5, k=10) {
  state.ps <- kf.dist.test(shrt=lrr,lng=adj,k=k,prior=prior)
  if(length(which(!is.na(state.ps)))>1) {
    state.ps <- state.ps[,!is.na(state.ps[1,])]
    state.ps <- matrix(state.ps,nrow=2)
    k.ps <- numeric(ncol(state.ps))
    k.ps[state.ps[1,]==copy] <- state.ps[2,][state.ps[1,]==copy]
    k.ps[state.ps[1,]!=copy] <- 1-state.ps[2,][state.ps[1,]!=copy]
    result <- mean(k.ps,na.rm=T)
    error <- sd(k.ps,na.rm=T)
    SE <- error/sqrt(k)
    CI <- (c(-1.96,1.96)*SE) + result
  } else {
    result <- SE <- NA; CI <- c(NA,NA)
  }
  return(c(p=result,SE=SE,lower=CI[1],upper=CI[2]))
}


most.likely.copy <- function(baf=NULL,lrr=NULL,test=NA,lik.mat=NULL,prior=0.5,warn.only=F,...) {
  # returns most likely copy number for a baf sample, with confidence versus next best option
  # if using test=x, this tests a specific copy-x against the most likely alternative
  no.result <- c(copy=NA,conf=NA)
  if(!is.numeric(baf) & !is.numeric(lrr)) { warning("one of baf or lrr must be a vector of numeric values"); return(no.result) }
  if(!is.null(baf) & !is.null(lrr)) { warning("only one of baf or lrr should be used, defaulting to BAF") }
  if(is.numeric(baf)) { baf.on <- T } else { baf.on <- F }
  if(is.null(lik.mat) & baf.on) { lik.mat <- get.baf.lik(n=100000,noise=2,recalc=F) }
  n.points <- max(c(length(baf),length(lrr)),na.rm=T)
  if(n.points<2) { warning("insufficient data") ; return(no.result) }
  if(baf.on) {
    csl <- copy.state.baf.lik(baf,lik.mat=lik.mat,...)
  } else {
    csl <- copy.state.lrr.lik(lrr,...)
  }
  csl <- as.numeric(csl)
  if(length(csl)==0) { warning("copy state data empty") ; return(no.result) }
  if(is.na(test)) {
    max.ll <- which(csl==max(csl))
  } else { max.ll <- test+1 }
  next.ll <- (c(1:length(csl))[-max.ll])[which(csl[-max.ll]==max(csl[-max.ll]))]
  x1 <- csl[max.ll]; x2 <- csl[next.ll]; 
  cnd <- x1>=x2; if(length(cnd)<1)  { return(no.result) }
  if(cnd) {
    ratio <- exp(abs(x2-x1)) # ratio of most likely prob vs next most prob
    if(is.infinite(ratio)) { ratio <- 10^100 }
    # prior <- .5 chance of genuine CNV [correct state inference] vs not
    conf <- (ratio*prior)/((ratio*prior) + ((1-prior)*1))
    ##conf <- pchisq(D.next,2)  #probability value of most likely (1 - null-p)
  } else {
    ratio <- exp(-abs(x1-x2)) ##this should happen only when !is.na(test)
    conf <- (ratio*(1-prior))/((ratio*(1-prior)) + ((prior)*1))
  }
  result <- max.ll-1
  if(n.points < 3)  {
    if(warn.only) {
      warning("copy likelihood test run for only ",n.points," points") 
    } else {
      # not enough points for a valid test
      return(c(copy=result,conf=NA))
    }
  }
  return(c(copy=result,conf=conf))
}


#max.n <- function(x,n) { rev(sort(x))[1:n] }

copy.state.baf.lik <- function(baf,lik.mat=NULL,state.priors=c(.001,.33,.33,.33,.001),ext.filt=F) {
  if(is.null(lik.mat)) { lik.mat <- get.baf.lik(n=100000,noise=2,recalc=F) }
  min.n <- function(x,n) { sort(x)[1:n] }
  baf <- round(baf,2)
  obs <- table(baf)
  obs.mat <- matrix(rep(obs,ncol(lik.mat)),ncol=ncol(lik.mat))
  the.mat <- lik.mat[names(obs),]*obs.mat
  nm <- nrow(the.mat); if(is.null(nm)) { nm <- 1 }
  if(ext.filt & length(baf)>9 & nm>6) {
    colmaxs <- apply(the.mat,2,max,na.rm=T)
    if(length(baf)<30 | nm<8) {
      colmins <- apply(the.mat,2,min,na.rm=T)
    } else {
      # remove percentage of points when many snps
      colmins <- apply(the.mat,2,min.n,n=ceiling(length(baf)/20))
      if(dim(colmins)[2]!=length(colmaxs[2])) {
        if(dim(colmins)[1]==length(colmaxs[2])) { colmins <- t(colmins) ; print("yes") }
      } else { print("no") }
      colmins <- colSums(colmins)
    }
  } else { colmins <- colmaxs <- 0 }
  #print(colSums(the.mat)); print(colmins); print(colmaxs)
  liks <- colSums(the.mat) - colmins - colmaxs
  if(is.numeric(state.priors)) { if(length(state.priors)==5) { liks <- liks + log(state.priors) } }
  names(liks) <- paste(0:4)
  return(liks)
}



copy.state.lrr.lik <- function(lrr,min.p=10^-6,reg.mn=NA,reg.sd=NA,state.priors=c(.001,.33,.33,.33,.001)) {
  mnz <- lrr.by.copy(); sdz <- lrr.by.copy(SD=T)
  if(is.numeric(reg.mn) & is.numeric(reg.sd)) {
    if(!is.na(reg.mn) & !is.na(reg.sd)) {
      mnz[3] <- reg.mn; sdz[3] <- reg.sd # punch in custom normal copy mean & sd
    } }
  lik.mat <- matrix(ncol=length(mnz),nrow=length(lrr))
  for (j in 1:length(mnz)) { lik.mat[,j] <- prob.norm(lrr,mu=mnz[j],sig=sdz[j]) }
  lik.mat[lik.mat==0] <- min.p
  liks <- colSums(log(lik.mat))
  if(is.numeric(state.priors)) { if(length(state.priors)==5) { liks <- liks + log(state.priors) } }
  names(liks) <- paste(0:4)
  return(liks)
}


est.state.priors <- function(in.cnv=T,near.cnv=F,cnv.copy=1,accuracy=.5,power=1) {
  # if it 's in the cnv then likely to be that copy state
  # near cnv only used when in.cnv =F; in.cnv set to F if cnv.copy=2
  a <- accuracy; b <- 1-a; if(cnv.copy==2) { in.cnv <- F }
  cnv.copy <- min(round(cnv.copy),4)
  if(!cnv.copy %in% 0:4){ stop("invalid cnv.copy") }
  if(in.cnv) {
    if(cnv.copy==0) { prior <- c(a,.65*b,.349*b,.00099*b,.00001*b) }
    if(cnv.copy==1) { prior <- c(.002*b,a,.987*b,.01*b,.001*b) }
    if(cnv.copy==3) { prior <- c(.001*b,.002*b,.987*b,a,.01*b) }
    if(cnv.copy==4) { prior <- c(.00001*b,.00099*b,.349*b,.65*b,a) }  
  } else {
    if(near.cnv) {
      if(cnv.copy==0) { prior <- c(.75*b,.24998*b,a,.00001*b,.00001*b) }
      if(cnv.copy==1) { prior <- c(.00998*b,.99*b,a,.00001*b,.00001*b) }
      if(cnv.copy==2) { prior <- c(.001*b,.499*b,a,.499*b,.001*b) }
      if(cnv.copy==3) { prior <- c(.00001*b,.00001*b,a,.99*b,.00998*b) }
      if(cnv.copy==4) { prior <- c(.00001*b,.00001*b,a,.24998*b,.75*b) }
    } else {
      # away from cnv, should really be copy=2 or copy=1 due to LOH [p's from Penn HMM: hh550.hmm]
      prior <- c(0.000001064,0.001165429,0.998795591,0.000012479,0.000000912) 
    }
  }
  prior <- (prior^power)/sum((prior^power)) # strengthen or weaken extreme priors (ensure sums to 1)
  names(prior) <- paste(0:4)
  return(prior)
}



lrr.by.copy <- function(SD=F,penn.file=F,...) {
  # the expected LRR mean for each copy number
  if(penn.file) { 
    if(!SD) {
      lrr.by.copy <- read.penn.hmm.pars(...)$mean
    } else {
      lrr.by.copy <- read.penn.hmm.pars(...)$sd
    }
    if(is.null(lrr.by.copy)) {
      penn.file <- F
    }
  }
  if(!penn.file) {
    if(!SD) {
      #Mean
      lrr.by.copy <- c(-3.527211,-0.664184,0.000000,0.395621,0.678345)
    } else {
      #SD
      lrr.by.copy <- c(1.329152,0.284338,0.159645,0.209089,0.191579)
    }
  }
  names(lrr.by.copy) <- paste(0:4)
  return(lrr.by.copy)
}


read.penn.hmm.pars <- function(file="/usr/local/bin/penncnv64/lib/hh550.hmm",mn=T,sd=T) {
  if(!file.exists(file)) { 
    mnz <- NULL; sdz <- NULL 
  } else {
    hmm <- readLines(file)
    mn.ln <- grep("B1_mean",hmm)+1
    sd.ln <- grep("B1_sd",hmm)+1
    mnz <- strsplit(hmm[mn.ln]," ")[[1]][-4]  # -4 is to remove LOH state
    sdz <- strsplit(hmm[sd.ln]," ")[[1]][-4]
    mnz <- as.numeric(mnz); sdz <- as.numeric(sdz)
  }
  if(mn & sd) { return(list(mean=mnz,sd=sdz)) }
  if(mn) { return(mnz) } else { return(sdz) }
}


calc.quality.stats <- function(DEL, bigBAF, bigPCC, snp.info, pr.acc=0.5, 
                               k=10, power=1, short.ratio=2, skip.sanity=F, nsnp=NA, do.chr=T) {
  ## get the raw data
  cat("\nDeriving CNV quality scores\n")
  cat(" extracting data for CNVs and adjacent regions\n")
  cnv <- big.extract.snp.ranges(range.snp(snp.info,DEL),samples=DEL$id,bigPCC)
  cnvbaf <- big.extract.snp.ranges(range.snp(snp.info,DEL),samples=DEL$id,bigBAF,snp.info=snp.info)
  if(do.chr) {
    #skipping chromosome QS reduces need to grab large chunks of data, reduces RAM req'mnt
    chrN <- chrIndices(snp.info)
    chrSnps <- cbind(rownames(snp.info)[chrN[,1]],rownames(snp.info)[chrN[,2]])
    chrDat <- big.extract.snp.ranges(chrSnps[as.numeric(chr(DEL)),],samples=DEL$id,bigPCC)
  }
  flanking1mb <- get.flanks.from.big.mat(DEL,bigPCC,bp=10^6,snp.info=snp.info, nsnp=NA)
  #print(head(flanking1mb))
  flanking2 <- get.flanks.from.big.mat(DEL,bigPCC,ratio=short.ratio,snp.info=snp.info, nsnp=nsnp)
  #print(head(flanking2))
  #flanking5 <- get.flanks.from.big.mat(DEL,bigPCC,ratio=5,snp.info=snp.info)
  fl <- get.ratio.set(DEL,ratio=short.ratio)
  FLK <- RangedData(ranges=IRanges(start=fl[,1],end=fl[,4]),space=chr(DEL),id=DEL$id)
  FLK <- toGenomeOrder(FLK,strict=T)
  bafDat2R <- get.flanks.from.big.mat(DEL,bigBAF,ratio=short.ratio,snp.info=snp.info,L=F, nsnp=nsnp)
  #print(head(bafDat2R))
  bafDat2L <- get.flanks.from.big.mat(DEL,bigBAF,ratio=short.ratio,snp.info=snp.info,R=F, nsnp=nsnp)
  bafDatExtAdj <- get.flanks.from.big.mat(FLK,bigBAF,ratio=1,snp.info=snp.info, nsnp=nsnp)
  #print(head(bafDatExtAdj))
  ## analyse the data
  cat(" analyzing called CNVs and adjacent regions\n")
  LRRs <- REGs <- BAFs <- BAF_REGs <- LENs <- list(); MONO <- logical()
  n.cnvs <- nrow(DEL)
  for (cnt in 1:n.cnvs) {
    #cnt <- 1
    copynum <- DEL$cn
    if(copynum[cnt]<2) { alt <- "greater" } else { alt <- "less" }
    near.test <- test.lrr.as(lrr=cnv[[cnt]], adj=flanking2[[cnt]], copy=as.numeric(copynum[cnt]), prior=pr.acc, k=k)  
    if(skip.sanity) {
      near.2.test <- 1
    } else {
      near.2.test <- cnv.sanity.test(shrt=cnv[[cnt]],lng=flanking2[[cnt]],DEL=(as.numeric(copynum[cnt])<2))
    }
    
    far.test <- test.lrr.as(lrr=cnv[[cnt]], adj=flanking1mb[[cnt]], copy=as.numeric(copynum[cnt]), prior=pr.acc, k=k) 
    if(do.chr) {
      chr.test <- test.lrr.as(lrr=cnv[[cnt]], adj=chrDat[[cnt]], copy=as.numeric(copynum[cnt]), prior=pr.acc, k=k) 
      far.vs.chr.test <- test.lrr.as(lrr=flanking1mb[[cnt]], adj=chrDat[[cnt]], copy=2, prior=pr.acc, k=k)
    }
    near.vs.far.test <- test.lrr.as(lrr=flanking2[[cnt]], adj=flanking1mb[[cnt]], copy=2, prior=pr.acc, k=k)    
    lrr.exp <- lrr.by.copy()[paste(copynum[cnt])]
    # lrr.const.test <- t.kf.test.fix(cnv[[cnt]],lrr.exp,alternative=alt)
    cnv.esp <- est.state.priors(T,T,copynum[cnt],accuracy=pr.acc,power=power)
    copy.prob <- suppressWarnings(most.likely.copy(baf=cnvbaf[[cnt]],test=as.numeric(copynum[cnt]),
                                  state.priors=cnv.esp,prior=pr.acc))
    mono <- (all(narm(c(cnvbaf[[cnt]], bafDat2L[[cnt]], bafDat2R[[cnt]]))<.02)) #monomorphic
    adj.esp <- est.state.priors(F,T,copynum[cnt],accuracy=pr.acc,power=power)
    adj.prob2L <- suppressWarnings(most.likely.copy(baf=bafDat2L[[cnt]],test=2,
                                   state.priors=adj.esp,prior=pr.acc))
    adj.prob2R <- suppressWarnings(most.likely.copy(baf=bafDat2R[[cnt]],test=2,
                                   state.priors=adj.esp,prior=pr.acc))
    far.esp <- est.state.priors(F,F,copynum[cnt],accuracy=pr.acc,power=power)
    adj.prob5 <- suppressWarnings(most.likely.copy(baf=bafDatExtAdj[[cnt]],test=2,
                                  state.priors=far.esp,prior=pr.acc))
    LRRs[[cnt]] <- c(near.test["p"]*near.2.test, far.test["p"], if(do.chr) { chr.test["p"] } else { NULL } )
    REGs[[cnt]] <- c(near.vs.far.test["p"],{ if(do.chr) { far.vs.chr.test["p"] } else { NULL }})
    BAFs[[cnt]] <- c(copy.prob[2])
    MONO[cnt] <- mono # is it monomorphic based on BAF
    BAF_REGs[[cnt]] <- c(adj.prob2L[2], adj.prob2R[2], adj.prob5[2])
    loop.tracker(cnt,n.cnvs)
  }
  LENs <- cbind(sapply(cnv,length),sapply(cnvbaf,length),sapply(flanking2,length),
                sapply(flanking1mb,length),{if(do.chr) { sapply(chrDat,length)} else { NULL }}, 
                sapply(bafDat2L,length), sapply(bafDat2R,length), sapply(bafDatExtAdj,length))
  colnames(LENs) <- c("cnv","cnvbaf","fl2","flmb",{if(do.chr) { "chr"} else { NULL}},"lf","rt","outer")
  return(list(lrr=LRRs, lrr.reg=REGs, baf=BAFs, baf.reg=BAF_REGs, len=LENs, mono=MONO))
}


#LRRs; REGs; BAFs; BAF_REGs
#



one.cnv.report <- function(cnt) {
  cat("\nCNV",cnt,"\n")
  # check whether BAF of CNV and surrounding area checks out
  if(all(!is.na(copy.prob[2]))) { if(copy.prob[2]<.5) { more.lik <- most.likely.copy(cnvbaf[[cnt]],state.priors=cnv.esp,prior=pr.acc); 
                                                        warning("In cnv #",cnt,", for cnv region, copy ",round(more.lik[1]),
                                                                " had ",round(more.lik[2],4)," probability") 
  } else {
    cat("CNV BAF matched copy",round(copy.prob[1]),"with likelihood",round(100*copy.prob[2],2),"%\n")
  } }
  if(all(!is.na(adj.prob2L[2]))) { if(adj.prob2L[2]<.5) { DD <- bafDat2L[[cnt]]; more.lik <- most.likely.copy(DD,state.priors=adj.esp,prior=pr.acc) ;
                                                          warning("In cnv #",cnt,", for left immediate region (",length(DD)," SNPs), copy ",round(more.lik[1]),
                                                                  " had ",round(more.lik[2],4)," probability (may be LOH if copy=1)") 
  } else {
    cat("In cnv #",cnt,"left immediate region BAF matched copy 2 with likelihood",round(100*adj.prob2L[2],2),"%\n")
  } }
  if(all(!is.na(adj.prob2R[2]))) { if(adj.prob2R[2]<.5) { DD <- bafDat2R[[cnt]]; more.lik <- most.likely.copy(DD,state.priors=adj.esp,prior=pr.acc) ;
                                                          warning("In cnv #",cnt,", for right immediate region (",length(DD)," SNPs), copy ",round(more.lik[1]),
                                                                  " had ",round(more.lik[2],4)," probability (may be LOH if copy=1)") 
  } else {
    cat("In cnv #",cnt,"right immediate BAF matched copy 2 with likelihood",round(100*adj.prob2R[2],2),"%\n")
  } }
  if(all(!is.na(adj.prob5[2]))) { if(adj.prob5[2]<.5) { DD <- bafDatExtAdj[[cnt]]; more.lik <- most.likely.copy(DD,state.priors=far.esp,prior=pr.acc) ;
                                                        warning("In cnv #",cnt,", for expanded adjacent region (",length(DD)," SNPs), copy ",round(more.lik[1]),
                                                                " had ",round(more.lik[2],4)," probability (may be LOH if copy=1)") 
  } else {
    cat("In cnv #",cnt,"expanded adjacent region BAF matched copy 2 with likelihood",round(100*adj.prob5[2],2),"%\n")
  } }
  # cat("\nmax test results CNV vs expected value of",lrr.exp,"at copy",copynum[cnt]," (NS) \n")
  # print((lrr.const.test))    #2=max p value
  cat("\nmean, CI test results CNV vs close surrounding region\n")
  print((near.test))
  cat("\nmean, CI test results CNV vs 1Mb surrounding region\n")
  print((far.test))
  cat("\nmean, CI test results CNV vs whole chromosome LRR\n")
  print((chr.test))
  cat("\nmean, CI test results close surrounding region LRR vs 1Mb surrounding region LRR (NS)\n")
  print(round((near.vs.far.test),5))
  cat("\nmean, CI test results 1Mb surrounding region LRR vs whole chromosome LRR (NS)\n")
  print(round((far.vs.chr.test),5))
}



get.quality.scores <- function(ranges,dir,n.pcs=NA,...) {
  dir <- validate.dir.for(dir,"big")
  DT <- read.data.tracker(dir)
  DEL <- ranges
  snp.info <- read.snp.info(dir)
  snp.info <- snp.info[snp.info$QCfail==0,]
  snp.info <- toGenomeOrder(select.autosomes(snp.info))
  bigPCC <- getSlot(DT,"big.pcc",ret.obj=T,n.pcs=n.pcs)
  bigBAF <- getSlot(DT,"big.baf",ret.obj=T,n.pcs=n.pcs)
  #print(dim(bigPCC)); print(dim(bigBAF))
  # YY <- calc.quality.stats(DUP, bigBAF, bigPCC, snp.info, pr.acc=0.5, k=10, power=1, short.ratio=3)
  XX <- calc.quality.stats(DEL, bigBAF, bigPCC, snp.info, pr.acc=0.5, k=10, power=1, short.ratio=2.5,...)
  # XX$lrr XX$lrr.reg XX$baf XX$baf.reg XX$len #[("cnv","cnvbaf","fl2","flmb","chr","lf","rt","outer")] XX$mono
  return(XX)
}


make.qs.table <- function(XX) {
  LRR.QS <- sapply(XX$lrr,mean)
  BAF.QS <- unlist(XX$baf) #^((XX$len[,"cnv"])^.25)
  CNV.QS <- rowMeans(cbind(LRR.QS,BAF.QS),na.rm=T)
  #CNV.QS[as.logical(XX$mono)] <- CNV.QS[as.logical(XX$mono)]*.98  # if monomorphic and BAF is high, ignore BAF
  CNV.QS[as.logical(XX$mono) & (BAF.QS>.95)] <- LRR.QS[as.logical(XX$mono) & (BAF.QS>.95)]  # if monomorphic and BAF is high, ignore BAF 
  #  BAF.QS[as.logical(XX$mono)] <- NA
  thr <- .8; nn <- length(XX$lrr)
  longer.than.tag <- in.roh.region <- local.bias <- rep(0,nn) 
  
  # likely CNV is longer than what is tagged
  cnd1 <- (sapply(XX$baf.reg,"[",3)<thr)
  cnd2 <- ((sapply(XX$baf.reg,"[",2)<thr) | (sapply(XX$baf.reg,"[",1)<thr))
  cnd3 <- ((sapply(XX$baf.reg,"[",2)<thr) & (sapply(XX$baf.reg,"[",1)<thr))
  longer.than.tag[ which(!cnd1 & cnd2 & !cnd3) ]  <- 1 
  longer.than.tag[ which(!cnd1 & cnd2 & !cnd3 & CNV.QS>.95) ]  <- 2
  cnd1 <- which(sapply(XX$lrr,"[",1)<thr)
  longer.than.tag[cnd1] <- longer.than.tag[cnd1] + 1
  cnd2 <- which(sapply(XX$lrr.reg,"[",1)<thr)
  longer.than.tag[cnd2] <- longer.than.tag[cnd2] + 1
  cnd3 <- sqrt(XX$len[,"cnv"])<8
  longer.than.tag[cnd3] <- longer.than.tag[cnd3] - 1
  longer.than.tag <- pmin(3,pmax(longer.than.tag,0))
  
  # likely CNV is in ROH region 
  cnd1 <- which(sapply(XX$baf.reg,"[",1)<thr)
  cnd2 <- which(sapply(XX$baf.reg,"[",2)<thr)
  in.roh.region[cnd1] <- 1 
  in.roh.region[cnd2] <- in.roh.region[cnd2]+1
  cnd <- which(sapply(XX$baf.reg,"[",3)<thr)
  in.roh.region[cnd] <- in.roh.region[cnd]+1
  in.roh.region <- pmin(3,pmax(in.roh.region,0))
  
  # likely there is local bias (wave interference) in the CNVR
  cnd <- which((sapply(XX$lrr,"[",2)<thr) | (sapply(XX$lrr,"[",3)<thr))
  local.bias[cnd] <- 2 
  cnd <- which((sapply(XX$lrr.reg,"[",2)<thr) | (sapply(XX$lrr.reg,"[",3)<thr))
  local.bias[cnd] <- local.bias[cnd] + 2
  local.bias <- local.bias - ceiling(longer.than.tag/2); local.bias <- pmin(3,pmax(local.bias,0))
  fs <- c("unlikely", "possible", "likely", "probable")
  
  # quality summary for all DELs/DUPs
  QS.results.matrix <- data.frame(cnvqs=CNV.QS, lrrqs=LRR.QS, bafqs=BAF.QS, n.snps=XX$len[,"cnv"],
                                  trunc=fs[longer.than.tag+1], roh=fs[in.roh.region+1], wave=fs[local.bias+1])
  return(QS.results.matrix)
}


########### DEFINE USER FUNCTIONS ############

boundary.stats <- function(lower=F,upper=F,nSD=3) {
  # return all stats, or potential lower/upper bounds
  # lower/upper bounds go from least extreme to most extreme (least deleterious) cutoff
  if(!is.numeric(nSD)) { nSD <- 3 }
  nsd <- paste(c("-","+"),paste(round(nSD,2),"SD",sep=""),sep="")
  low <- c("Q1","-2SD","Low1%","LB",nsd[1],"Min")
  upp <- c("Q3","+2SD","Hi1%", "UB",nsd[2],"Max")
  all <- c("Mean","SD","Min","Q1","Median","Q3","Max",
           "LB","UB","-2SD","+2SD",nsd[1],nsd[2],"Low1%","Hi1%")
  if(lower & !upper) { return(low) }
  if(!lower & upper) { return(upp) }
  return(all)
}


select.autosomes <- function(snp.info) {
  # select autosomes only in RangedData object
  must.use.package("genoset",bioC=T)
  if(!is(snp.info)[1]=="RangedData") { warning("not a RangedData object"); return(snp.info) }
  if(length(unique(space(snp.info))) < length(levels(space(snp.info)))) {
    # this fixes the problem when a subset of a ranges object with less 
    #  chromosomes still has empty chr slots from previous object
    snp.info <- snp.info[as.numeric(unique(space(snp.info)))]
  }
  select <- which(rownames(chrInfo(snp.info)) %in% paste(1:22))
  ok.chrs <- rownames(chrInfo(snp.info))[select]
  return(snp.info[ok.chrs])
}


chr.nums.of.info <- function(snp.info,warn=F) {
  must.use.package("genoset",bioC=T)
  if(!is(snp.info)[1]=="RangedData") { warning("not a RangedData object"); return(NULL) }
  txt <- rownames(chrIndices(snp.info))
  nums <- suppressWarnings(as.numeric(txt))
  num.na <- length(nums[is.na(nums)])
  if(num.na>0) { 
    if(warn) { warning(paste("chr.nums requested for non-autosomes, will assign numbers >23 to letters",
                  paste(txt[is.na(nums)],collapse=","))) }
    nums[is.na(nums)] <- 24:(24+num.na-1)
  }
  return(nums)
}


calc.cov <- function (snp.info, targ.int, min.snps) {
  # function used by chip.coverage() to calculate coverage gaps
  cur.chr <- chr(snp.info)
  if(nrow(snp.info)<10) { warning("not enough snps to calculate coverage"); return(NA) }
  nsnp <- nrow(snp.info)-min.snps;
  beg <- fin <- Chr <- integer(nsnp) 
  lenz <- dif <- numeric(nsnp)
  startz <- start(snp.info)
  skipnum <- 0
  for(cc in 1:nsnp) { 
    # loop through each snp, looking at that snp pos, and the pos 'min.snps' further along
    if( skipnum<=1) {
      lenz[cc] <- startz[cc+min.snps]-startz[cc] # distance between the positions
      if(lenz[cc]>targ.int)
      {
        # if the 'min.snps' distance is longer than 'targ.int'
        dif[cc] <- 0
        gap <- lenz[cc]-targ.int # difference between interval and targ.int
        while(dif[cc]<gap) {
          # see how many snps we need to skip ahead until the gap to 10th [min.snp] snp is <= targ.int
          # record this distance as a gap where a CNV can't be detected, then start again from the 'skipnum' snp
          skipnum <- skipnum+1
          dif[cc] <- startz[cc+skipnum] - startz[cc]   ## old idea, don't think it works : dif2 <- startz[cc+min.snps+skipnum+1] - startz[min.snps+cc]
        }
        beg[cc] <- startz[cc]; fin[cc] <- startz[cc+skipnum]; Chr[cc] <- cur.chr[cc]
      }
    } else {
      # this counter helps to skip snps which are too far apart to call a targ.int size CNV with min.snps
      skipnum <- skipnum-1
    }
  }
  ## if using the option to export the gaps as ranges
  sel <- !is.na(beg) & !is.na(fin) & beg!=0
  beg <- beg[sel] ; fin <- fin[sel] ; Chr <- Chr[sel]
  rd <- RangedData(IRanges(start=beg,end=fin),space=Chr)
  rd <- toGenomeOrder(rd,strict=T)
  return(rd)
}


chip.coverage <- function(snp.info,targ.int=100000,by.chr=F,min.snps=10,
                          full.chr.lens=T,verbose=T,dir="",ucsc="hg18",ranges=F) {
  # calculate % of chromosomes/genome covered by a set of microarray markers in terms of
  # possible CNVs detectable
  # targ.int #100000  ## size of windows to find
  # min.snps ## minimum number of snps per window
  # can return the gaps - they'll be unmerged...
  must.use.package(c("genoset","IRanges"),T)
  if(!is(snp.info)[1]=="RangedData") { warning("snp.info wasn't RangedData") ; return(NULL) }
  snp.info <- select.autosomes(snp.info)
  snp.info <- toGenomeOrder(snp.info,strict=T)
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(snp.info)
  # single core -  multicore doesn't really help much (like 5 seconds saving for 22 cores for a large dataset)
  rd <- calc.cov(snp.info, targ.int, min.snps)
  rd <- reduce(rd) # combine adjacent ranges
  # calculate total length as total length of all chromosomes
  if(is.null(dir)) { dir <- "" }
  chr.lens.full <- (get.chr.lens(dir,ucsc=ucsc)[chr.set])-(targ.int)
  # calculate total  length using the first and last pos on each chromosome [allows calc for subranges]
  chr.lens.range <- as.numeric(sapply(snp.info,function(x) { diff(range(start(x))) })) #- (targ.int)
  if (full.chr.lens) { 
    # add gaps before the first snp and after the last snp for each chr
    chr.lens <- chr.lens.full 
    stz <- sapply(snp.info,function(x) { min(start(x)) }) #- (targ.int)
    enz  <- (chr.lens.full - sapply(snp.info,function(x) { rev(sort(start(x)))[min.snps-1] })) - (targ.int)
    rd2 <- RangedData(IRanges(start=rep(1,n.chr)[stz>0],end=stz[stz>0]),space=chr.set[stz>0])
    rd3 <- RangedData(IRanges(start=(chr.lens.full-enz-targ.int)[enz>0],end=(chr.lens.full-targ.int)[enz>0]),space=chr.set[enz>0])
    rd  <- rbind(toGenomeOrder(rd,strict=T),toGenomeOrder(rd2,strict=T),toGenomeOrder(rd3,strict=T))
  } else { 
    chr.lens <- chr.lens.range  
  }
  rd <- toGenomeOrder(rd,strict=T)
  gappoz <- sum(as.numeric(width(rd)))
  if(by.chr) {
    ind.chrs <- 1-(sapply(rd,function(x) { sum(width(x)) })/chr.lens)
  }  else { ind.chrs <- c() }
  total.chr <- sum(chr.lens)
  cov <- 1-(gappoz/total.chr)
  if(verbose) {
    if(full.chr.lens) { onl <- ", " } else { onl <- "\n" }
    cat("found ",nrow(rd)," gaps >",targ.int,"K, not covered by at least ",min.snps," snps\n",sep="")
    cat("combined length of gaps: ",round(gappoz/10^6,2),"MB\n",sep="")
    cat("overall proportion of autosomes with gaps: ",gappoz/(total.chr),onl,"total autosomal genome size: ",round(total.chr/10^6,2),"MB",sep="")
    if(full.chr.lens) { cat("\n") } else { cat("*\n *ignoring regions before and after first and last SNPs of each chromosome\n")}
    cat("total CNV coverage:",cov,"\n")
  }
  outlist <- list(total=cov,chr=ind.chrs,gaps=rd)
  if(!by.chr) { outlist[["chr"]] <- NULL }
  if(!ranges) { outlist[["gaps"]] <- NULL }
  if(length(outlist)==1) { outlist <- unlist(outlist) } # simplify if only overall coverage is desired
  return(outlist)
}


coverage.change <- function(snp.info,cnv.sizes=c(100,500),min.snps=10,dir=NULL,
                            plot.fn="coverageChangesSnp",do.plot=T) {
  # report on loss of genome coverage due to QC excluded snps #
  must.use.package(c("genoset","IRanges"),T)
  if(is(snp.info)[1]!="RangedData") { warning("snp.info is wrong type") ; return(NULL) }
  snp.info <- select.autosomes(snp.info)
  n.snps.rmv <- paste(length(which(snp.info$QCfail!=0)),"/",nrow(snp.info),sep="")
  nL <- length(cnv.sizes); all.cov <- qc.cov <- vector("list",nL)
  cat(" after SNP-QC removed",n.snps.rmv,"SNPs, coverage of the genome has decreased from: \n")
  chr.lst <- chr.nums.of.info(snp.info)
  plot.fn <- cat.path(dir$cr,plot.fn,suf=min.snps,ext="pdf")
  if(do.plot & !is.null(dir)) { pdf(plot.fn) }
  for (j in 1:nL) {
    all.cov[[j]] <- chip.coverage(snp.info,cnv.sizes[j]*10^3,min.snps=min.snps,by.chr=T,
                                  full.chr.lens=T,verbose=F,dir=dir)
    qc.cov[[j]] <- chip.coverage(snp.info[snp.info$QCfail==0,],cnv.sizes[j]*10^3,min.snps=min.snps,by.chr=T,
                                 full.chr.lens=T,verbose=F,dir=dir)
    cat("     ",round(100*all.cov[[j]][[1]],2),"% to ",round(100*qc.cov[[j]][[1]],2),"%",
        " for CNVs of size ",cnv.sizes[j],"K\n",sep="")
    yl <- min(1.1,max(all.cov[[j]][[2]])*1.5) # make y-limits big enough for legend
    barplot(rbind(qc.cov[[j]][[2]],(all.cov[[j]][[2]]-qc.cov[[j]][[2]])),
            xlab="chromosome",ylab="% coverage",ylim=c(0,yl),names.arg=paste(chr.lst),
            main=paste("Coverage for",cnv.sizes[j],"K CNVs, pre vs post SNP-QC\n[filtered",n.snps.rmv,"SNPs]"))
    legend("top",legend=paste("coverage",c("pre","post"),"QC"),pch=22,pt.cex=2,pt.bg=c("grey90","grey32"),bty="n",ncol=2)
  }  
  if(do.plot & !is.null(dir)) {
    dev.off()
    cat("~wrote plot",plot.fn,"\n")
  }
}


write.snp.info <- function(snp.info,dir,...) { 
  write.sample.info(dir=dir,sample.info=snp.info,type="bin",fn="snpinfo",...) 
}


write.data.tracker <- function(DT,dir=NULL,fn="datatracker",...) { 
  if(is.data.tracker(DT)) {
    if(is.null(dir)) { dir <- getSlot(DT,"dir") }
    write.sample.info(dir=dir,sample.info=DT,type="bin",fn=fn,...) 
  } else {
    stop("Error: couldn't write DT, not a data.tracker object")
  }
}


read.snp.info <- function(dir,...) { 
  read.sample.info(dir=dir,type="bin",fn="snpinfo",...) 
}


read.data.tracker <- function(dir,fn="datatracker",...,warn.only=F) { 
  DT <- read.sample.info(dir=dir,type="bin",fn=fn,...) 
  if(is.data.tracker(DT)) { return(DT) } else { 
    if(warn.only) {
      warning("returning NULL as file did not contain a datatracker object")
      return(NULL)
    } else {
      stop("Error: file did not contain a datatracker object")
    }
  }
}


write.sample.info <- function(sample.info,dir,verbose=F,type="tab",fn="sampleinfo",nprev=0) 
{
  # write standard 'sample.info' object to file
  dir <- validate.dir.for(dir,"ano")
  ofn <- cat.path(dir$ano,fn,ext=type)
  if(toupper(type) %in% c("RDATA","BIN","BINARY")) {
    ofn <- cat.path(dir$ano,fn,ext="RData")
    support <- sample.info # this line only here to make object name in file more appropriate
    save(support,file=ofn)
  } else {
    ofn <- cat.path(dir$ano,fn,ext="tab")
    write.table(sample.info,file=ofn,quote=F)
  }
  if(verbose) {
    cat(paste("==> full info object written to:",ofn,"\n"))
  }
  if(is.numeric(nprev)) {
    if(nprev>0) {
      cat("\nFirst",nprev,"records excerpt:\n")
      print(head(sample.info,nprev))
    } }
  return(ofn)
}


read.sample.info <- function(dir,type="tab",fn="sampleinfo",verbose=F,remove.fail=F,nprev=0) 
{
  # read standard 'sample.info' object from file
  dir <- validate.dir.for(dir,"ano")
  if(toupper(type) %in% c("RDATA","BIN","BINARY")) {
    ofn <- cat.path(dir$ano,fn,ext="RData")
    if(!file.exists(ofn)) { warning("No info binary file was found") ; return(NULL) }
    sample.info <- get(paste(load(ofn)))
  } else {
    ofn <- cat.path(dir$ano,fn,ext="tab")
    if(!file.exists(ofn)) { warning("No info tab file was found") ; return(NULL) }
    sample.info <- read.table(file=ofn)
  }
  if(verbose) {
    cat(paste("sample.info object read from:",ofn,"\n"))
  }
  if(is.numeric(nprev)) {
    if(nprev>0) {
      cat("\nFirst",nprev,"records excerpt:\n")
      print(head(sample.info,nprev))
    } }
  if(remove.fail) {
    qc.col <- which(toupper(colnames(sample.info)) %in% "QCFAIL")
    if(length(qc.col)>0) { sample.info <- sample.info[sample.info[,qc.col]!=1,] }
  }
  return(sample.info)
}


plink.to.Ranges <- function(plink.cnv) {
  # read in a plink cnv file and return Ranged data object
  if(get.ext(plink.cnv)!="cnv") {
    warning(paste("this function is meant for importing plink *.cnv files so using a file with",
                  "extension",get.ext(plink.cnv),"may cause unpredictable results\n"))
  }
  dat <- read.table(plink.cnv,header=T)
  # Plink format:
  #FID  IID  CHR  BP1  BP2	TYPE	SCORE	SITES
  outData <- RangedData(ranges=IRanges(start=dat$BP1,end=dat$BP2),id=dat$IID,space=dat$CHR,
                        cn=dat$TYPE,numSnps=dat$SITES,fid=dat$FID,score=dat$SCORE)
  outData <- toGenomeOrder(outData,strict=T)
  # note the IDs can't be rownames due to duplicates
  return(outData)
}


Ranges.to.cnvgsa <- function(dat) {
  # convert an Ranged data object to a cnv-gsa object
  # for use with the package cnvGSA
  if(!is(dat)[1]=="RangedData") { stop("This function only accepts objects of type 'RangedData") }
  # cnvGSA format below:
  #SampleID ID assigned to the subject’s DNA sample in which the CNV was found
  #Chr Chromosome on which the CNV is located
  #Coord_i Start position of the CNV on the chromosome
  #Coord_f End position of the CNV on the chromosome
  #Type CNV type (e.g. "DEL" or "DUP")
  #Genes (empty string)
  #CnvID ID assigned to the CNV
  types <- character(nrow(dat))
  types[dat$cn>2] <- "DUP"
  types[dat$cn<2] <- "DEL"  
  nmz <- rownames(dat); if(is.null(nmz)) { nmz <- paste(1:nrow(dat)) }
  cnvGSA <- data.frame(SampleID=dat$id,Chr=space(dat),Coord_i=start(dat),Coord_f=end(dat),
                       Type=types,Genes="",CnvID=nmz)
  cat(" cnvGSA format generated; this data.frame can go into the 'cnvData' slot of a CnvGSAInput object\n")
  cat(" use function full.cnvGSA() to make a full CnvGSAInput object using this data.frame\n")
  return(cnvGSA)
}


add.genes.to.GSA <- function(cnv.frame, delim=";", ucsc="hg18", genemap=NULL) {
  ## use cnvGSA function 'getCnvGenes' to add gene annotation to a cnv slot object for cnvGSA
  gsa.cols <- c("SampleID","Chr","Coord_i","Coord_f","Type","Genes","CnvID")
  if(!all(colnames(cnv.frame)==gsa.cols))
  { cat(" invalid frame, genes not added\n") ; return(cnv.frame) }
  must.use.package("cnvGSA",T)
  if(is.null(genemap)) {
    genemap <- get.gene.annot(ucsc=ucsc,range.out=F)
    gga.cn <- c("gene","chr","start","end","band")
    if(all(colnames(genemap)==gga.cn)) {
      genemap <- genemap[,c(2,3,4,1)]
      colnames(genemap) <- c("Chr","Coord_i","Coord_f","GeneID") # cnvGSA needs these names
    } else {
      cat("Error: colnames from annotation were:",paste(colnames(genemap),collapse=","))
      cat(" but needed:",paste(gga.cn,collapse=",")); stop()      
    }
    genes <- getCnvGenes(cnv.frame, genemap, delim=delim)
  }
  cnv.frame$Genes <- genes
  return(cnv.frame)
}


get.geneData.obj <- function() {
  # get object needed for cnvGSA package
  # contains lookup of SYMBOL vs Entrez ID, and descriptions of genes
  must.use.package("org.Hs.eg.db",T)
  ann <- list( gene2sy = character(0), gene2name = character(0) )
  x <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(x)
  ann$gene2sy <- unlist( as.list( x[mapped_genes] ) )
  x <- org.Hs.egGENENAME
  mapped_genes <- mappedkeys(x)
  ann$gene2name <- unlist( as.list( x[mapped_genes] ) )
  geneData_demo <- list(ann)
  rm( ann, x, mapped_genes )
  return(geneData)
}


full.cnvGSA <- function(do.assoc=F, cnv.frame=NULL, gmt.file=NULL, grandtotals.mode=c("all","cnv","cnvGen"),
                        sample.classes=c("case","ctrl"),fdr.iter=2,ext.report=200,boxplots=F,
                        limit.type=c("DEL","DUP"),limits.size=NULL,rem.genes=NULL,
                        delim=";", ucsc="hg18", genemap=NULL) {
  # create a cnvGSA object using cnv.frame or RangedData object generated by plumb cnv
  must.use.package("cnvGSA",T)
  if(is(cnv.frame)[1]=="RangedData") { 
    uv <- tolower(universe(cnv.frane)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
    cnv.frame <- Ranges.to.cnvgsa(cnv.frame) 
  } else {
    if(is.data.frame(cnv.frame)) {
      if(!all(colnames(cnv.frame)==gsa.cols))
      { cat(" invalid cnv.frame frame, cnvGSA object not produced\n") ; return(NULL) }
    } else {
      stop(" cnv.frame should be RangedData or a data.frame")
    }
  }
  ## cnv data.frame
  cnv.frame <- add.genes.to.GSA(cnv.frame, delim=";", ucsc=ucsc, genemap=NULL)
  ## gene data
  geneData <- get.geneData.obj()
  ## gmt file
  if(is.null(gmt.file)) {
    cat(" Parameter 'gmt.file' must be specified. See the following links to download/describe the GMT format\n")
    cat("http://www.geneontology.org/\n")
    cat("http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data formats\n")
    stop()
  } else {
    if(file.exists(gmt.file)) {  gsData <- readGMT(gmt.file) } else { stop(" gmt.file not found") }
  }
  ## params 
  if(is.character(grandtotals.mode) & is.character(sample.classes) &
    is.numeric(fdr.iter) & is.numeric(ext.report) & is.logical(boxplots)) {
    params <- list(grandtotals_mode=grandtotals.mode[1],sample_classes=sample.classes,
                   fdr_iter=fdr.iter[1],extended_report=ext.report[1],
                   filters=list(limits_type=limit.type,limits_size=limits.size,
                                rem_genes=rem.genes),boxplots_PDFs=boxplots[1]) 
  } else {
    warning("invalid format for params. Will leave params object as NULL. Set it using params()<-")
    params <- NULL
  }
  input <- CnvGSAInput(cnvData=cnv.frame, gsData=gsData, geneData=geneData,  params=params)
  if(do.assoc) {
    anal <- cnvGSAFisher(input); return(anal)
  } else {
    cat("CnvGSAInput object created. This can now be analysed using cnvGSAFisher()\n")
    return(input)
  }
}


plink.to.cnvGSA <- function(cnv.data) {
  # read in a plink cnv file and return cnvGSA object
  if(get.ext(cnv.data)!="cnv") {
    warning(paste("this function is meant for importing plink *.cnv files so using a file with",
                  "extension",get.ext(cnv.data),"may cause unpredictable results"))
  }
  dat <- plink.to.Ranges(cnv.data)
  dat <- Ranges.to.cnvGSA(dat)
  return(dat)
}


## export!
get.PCA.subset <- function(dir,pc.to.keep=.13,assoc=F,autosomes=T,big.fn="combinedBigMat.RData",
                           snp.sub.fn="pca.snp.subset.txt",use.current=F,pref="PCAMatrix",n.cores=1,
                           descr.fn="pcaSubMat.RData",nprev=0,snp.info=NULL,sample.info=NULL) 
{
  ## extract LRR matrix with subset of SNPs, ready for PCA analysis
  dir <- validate.dir.for(dir,c("ano","big"))
  #load.all.libs() # load all main libraries used by plumbCNV
  #if(add.pheno) {
  #  sample.info <- read.sample.info(dir)
  #  phenotype <- read.table(file=cat.path(dir$ano,"pheno.lookup.txt"))
  #  sample.info <- add.to.sample.info(sample.info,phenotype,"phenotype")
  #  write.sample.info(sample.info,dir)
  #}
  if(!is.data.frame(sample.info)) { sample.info <- read.sample.info(dir,nprev=nprev) }
  if(is(snp.info)[1]!="RangedData") { snp.info <- read.snp.info(dir,nprev=nprev) }
  #sample.info <- validate.samp.info(sample.info,QC.update=F,verbose=F) #this done done later anyway
  #samp.fn <- "combined.samples.txt"
  if(use.current & is.file(snp.sub.fn,dir$ano,dir)) {
    snps.to.keep <- get.vec.multi.type(snp.sub.fn,dir=dir)
  } else {
    snps.to.keep <- extract.snp.subset(snp.info,sample.info,pc.to.keep=pc.to.keep,assoc=assoc,autosomes=autosomes,
                                       writeResultsToFile=T,big.fn=big.fn,out.fn=snp.sub.fn,dir=dir, n.cores=n.cores)
  }
  ###bigMat <- getBigMat(big.fn,dir)
  if(length(snps.to.keep)>100) {
    ##writeLines(colnames(bigMat),paste(dir$ano,samp.fn,sep=""))
    subset.descr <- big.exclude.sort(big.fn,dir=dir,T,tranMode=1,pref=pref,f.snp=snps.to.keep,verbose=F)
  } else {
    stop("Error: list of snps to keep is too small - trying running again with a higher pc.to.keep\n")
  }
  #if(descr.fn!="") {
  #  save(subset.descr,file=cat.path(dir$big,descr.fn))
  #} else { warning("submatrix description returned but not saved\n")}
  print(subset.descr)
  return(subset.descr)
}



random.spacing.snp.select <- function(snp.info,pc.to.keep=.05,dir,autosomes.only=T) {  
  # this assumes roughly whole genome coverage. breaks down if this is not roughly true
  # uses the chr, pos, label of each snp, stored in a genoset RangedData object
  # to yield a % (e.g, 5%) subset of snps evenly spaced throughout the genome
  cat(" choosing spaced",round(pc.to.keep*100,1),"% subset of SNPs\n")
  # to use in PCA [only for Chr 6, exclude the MHC region]
  ## MHC chr6:30,018,000-33,606,563; build 37 only slighty earlier
  must.use.package("genoset",T)
  if(is(snp.info)[1]!="RangedData") { stop("Error: snp.info must be a RangedData object") }
  if(autosomes.only) { snp.info <- select.autosomes(snp.info); cat(" selecting only autosomes\n") }
  skip.chr.mhc <- 6; chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
  mhc <- c(29500000,34000000) # ok for build 36/37
  cat("snp.info contains data for chromsomes:",paste(chr.set,collapse=","),"\n")
  ratio = min(1,(1.15*pc.to.keep))  # upward adjust to counter empty bits and exclusions
  ngPCA <- (nrow(snp.info)*ratio) # number of intervals in genome
  av.int <- (chrInfo(snp.info)[length(snp.info),"stop"])/ngPCA # mean interval length
  chr.uniq <- list()
  sumo <- 0 # counter for non-unique mappings
  for (cc in 1:n.chr)
  {
    rr <- range(start(snp.info[(cc)])) # start and end position of chromosome
    l.out <- (rr[2]-rr[1])/av.int
    # create a sequence of evenly spaced locations in the chromosome
    lociz <- seq(from=((av.int/2)+rr[1]),to=(rr[2]-(av.int/2)),length.out=l.out)
    if (chr.set[cc]==skip.chr.mhc) {
      nskip <- length(which(lociz>mhc[1] & lociz<mhc[2]))
      lociz <- lociz[lociz<mhc[1] | lociz>mhc[2]]   
      cat(" skipped",nskip,"locations in MHC region\n")
    }
    # find closest SNPs to each location (lociz)
    lc <- (matchpt(lociz, start(snp.info[(cc)])))
    # calculate number of locations not matching uniquely
    not.uniq <- length(lc[,1])-length(unique(lc[,1]))
    ##cat("",paste("chr",cc,": for",not.uniq,"locations the most proximal snp was already used\n"))
    sumo <- sumo + not.uniq
    chr.uniq[[cc]] <- unique(lc[,1])
  }
  kpt <- sum(sapply(chr.uniq,length))
  
  cat(paste(" total number of SNPs kept:",kpt,"[",round(kpt/nrow(snp.info)*100,2),"%]\n"))
  cat(paste(" total of",sumo,"removed because they did not uniquely map to 1 location\n"))
  
  # combine snps from each chromosome get the full list of SNP IDs for PCA
  snp.to.pca <- character()
  for (cc in 1:n.chr) { snp.to.pca <- c( snp.to.pca,rownames(snp.info[cc])[chr.uniq[[cc]]] ) }
  return(snp.to.pca)
}


table2d <- function(...,col,row,rn=NULL,cn=NULL) {
  inp <- list(...); bad1 <- F 
  if(length(inp)<2) { bad1 <- T } else { 
    if(length(inp[[1]])!=length(inp[[2]])) { bad1 <- T 
    } else { if(length(inp[[1]])<1) { bad1 <- T } }
  }
  if(bad1) {
    warning("at least 2 arguments of equal length required for table2d, passing to table")
    return(table(...)) 
  }
  rr <- inp[[1]]; cc <- inp[[2]]
  # add 1 of each possible cell to the 2 vectors, to ensure each table cell is represented
  rr <- c(rep(sort(row),each=length(col)),rr)
  cc <- c(rep(sort(col),length(row)),cc)
  nmz <- sapply(match.call(expand.dots=TRUE)[-1], deparse) # get ... arg names
  TT <- table(rr,cc,dnn=nmz[1:2])-1
  if(is.character(rn)) { if(length(rn)==length(row)) { rownames(TT) <- rn } }
  if(is.character(cn)) { if(length(cn)==length(col)) { colnames(TT) <- cn } }
  return(TT)    
}


print.snp.sample.summary <- function(dir) {
  ## make simple summary of number of passing/failing SNPs at different processing stages
  sample.info <- read.sample.info(dir)
  snp.info <- read.snp.info(dir)
  samp.tab <- table(sample.info$QCfail); snp.tab <- table(snp.info$QCfail)
  cat("\n=== SNP SUMMARY ===\n")
  cat(nrow(snp.info),"SNPs in starting dataset\n")
  if(!is.na(snp.tab["2"])) { cat(snp.tab["2"],"SNPs failed SNP-QC\n") }
  if(!is.na(snp.tab["1"])) { cat(snp.tab["1"],"SNPs failed due to other criteria\n") }
  if(!is.na(snp.tab["0"])) { cat(snp.tab["0"],"SNPs passed QC\n") }
  cat("\n=== SAMPLE SUMMARY ===\n")  
  cat(nrow(sample.info),"samples in starting dataset\n")
  if(!is.na(samp.tab["2"])) { cat(samp.tab["2"],"samples failed SNP-QC\n") }
  if(!is.na(samp.tab["3"])) { cat(samp.tab["3"],"samples failed Sample-QC\n") }
  if(!is.na(samp.tab["4"])) { cat(samp.tab["4"],"samples failed during batch correction\n") }
  if(!is.na(samp.tab["6"])) { cat(samp.tab["6"],"samples failed CNV-QC\n") }
  if(!is.na(samp.tab["1"])) { cat(samp.tab["1"],"samples failed due to other criteria\n") }
  if(!is.na(samp.tab["0"])) { cat(samp.tab["0"],"samples passed QC\n") }
  pf <- sample.info$QCfail; pf[pf>0] <- 1; cn <- colnames(sample.info)
  
  if("grp" %in% cn) { gr <- sample.info$grp; if(length(unique(gr))>0) {
    cat("\n=== COHORT SUMMARY ===\n")
    grp.tab <- table2d(pf,gr,row=c(0,1),col=unique(gr),rn=c("pass","fail"),cn=paste("grp",unique(gr)))
    #rownames(grp.tab) <- c("pass","fail") ; colnames(grp.tab) <- paste("grp",unique(gr))
    pc.rmv <- (t(grp.tab)/colSums(grp.tab))[,2]
    print(grp.tab); cat(paste(round(100*pc.rmv,2),"% failed in grp=",
                 unique(gr),"\n",sep=""),sep="")
  }}
  if("phenotype" %in% cn) { ph <- sample.info$phenotype; if(length(unique(ph))>0) {
    cat("\n=== PHENOTYPE SUMMARY ===\n")
    #ph.tab <- table(pf,ph)
    ph.tab <- table2d(pf,ph,row=c(0,1),col=unique(ph),rn=c("pass","fail"),cn=paste("pheno",unique(ph)))
    #rownames(ph.tab) <- c("pass","fail"); colnames(ph.tab) <- paste("pheno",unique(ph))
    pc.rmv <- (t(ph.tab)/colSums(ph.tab))[,2]
    print(ph.tab)
    cat(paste(round(100*pc.rmv,2),"% failed in phenotype=",
                    unique(ph),"\n",sep=""),sep="")
  }}
}


calculate.gc.for.samples <- function(bigMat,snp.info,dir,med.chunk.fn="",restore.mode=F,ucsc="hg18",n.cores=1)
{
  # calculate GC-wave stats for each sample
  must.use.package(c("BiocGenerics","IRanges","Biobase"),bioC=T)
  ## parameters in order are:
  ## bigmatrix, directory with bigmatrix, snp.info object, gc.dat object, file name, regenerate?
  dir <- validate.dir.for(dir,c("gc","big"))
  # get genome wide human GC data at megabase resolution 
  # load average data for 10^6 windows prepared using 'extractGCwindows.R'
  gc.sav <- cat.path(dir$gc,"GCAv6.RData")
  uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  if(!file.exists(gc.sav))
  {
    gc.dat <- get.gc.human(10^6,ucsc=ucsc,n.cores=n.cores)
    save(gc.dat,file=gc.sav)
  } else {
    gc.dat <- get(paste(load(gc.sav)))
  }
  #####print(is(bigMat)); print(dim(bigMat))
  # get list of snps for each valid megabase windows that has at least 10 snps
  med.list <- get.ranges.for.median(snp.info=snp.info,gc.dat=gc.dat,bigMat=bigMat)
  rng <- med.list$gc.dat
  if(!(med.chunk.fn=="")) { 
    med.chunk.fn <- cat.path(dir$gc,med.chunk.fn)
    if(!file.exists(med.chunk.fn) | !restore.mode) { 
      gotMedian <- F 
    } else {
      ## if median matrix has been calculated previously load it and skip this step
      med.store <- get(paste(load(paste(med.chunk.fn))))
      if(!((ncol(med.store)==ncol(bigMat)) & (nrow(med.store)==nrow(rng))))
      { warning("can't used loaded object 'med.store', dim: ",paste(dim(med.store),collapse=","),
                ". Object has the wrong dimensions versus other parameters")
        gotMedian <- F
      } else {
        gotMedian <- T
      }
    }
  } else {
    cat(" Warning: missing 'med.chunk.fn'; may result in unnecessary recalculation of gc%\n")
    gotMedian <- F
  }
  if(!gotMedian | !exists("med.store")) {
    ## calculate the medians in each range for each sample, about 1hr per 100K-snps*5K-samples
    med.store <- do.median.for.ranges(ranges.list=med.list[[1]],bigMat=bigMat,dir=dir,cont=!med.list[[3]],
                                      use.big.list=med.list[[2]], med.sav.nm=med.chunk.fn,n.cores=n.cores)
  } 
  r_GC <- S_WF <- S_GCWF <- numeric(ncol(bigMat))
  cat(" calculating Diskin constants...\n")
  for (j in 1:ncol(bigMat)) {
    r_GC[j] <- cor(rng$gc,med.store[,j],use="pairwise.complete.obs")
    S_WF[j] <- (1-2*(as.numeric(r_GC[j]<0))) *
      median(abs(med.store[,j]-median(med.store[,j],na.rm=T)),na.rm=T)
  }
  S_GCWF <- S_WF * abs(r_GC)
  dat.out <- cbind(r_GC,S_WF,S_GCWF)
  rownames(dat.out) <- colnames(bigMat)
  colnames(dat.out) <- c("r_GC","S_WF","S_GCWF")
  return(dat.out)
}


make.sample.info <- function(dir="",id.list=NULL,phenotype=NULL,plate=NULL,grp=NULL,sex=NULL,
                             other=list(),overwrite=T,verbose=F,...) {
  # create sample info dataframe, combining an id list with phenotype, plate, group,
  # or any other data, inputted as a named list of objects or file names.
  dir <- validate.dir.for(dir,"ano")
  pheno.default.loc <- "pheno.lookup.txt"; phenot <- NULL
  if(is.null(id.list) & !is.null(dir)) {
    id.list <- get.subIDs(dir,"combine",verbose=verbose)
  } else {
    if(is.file(id.list,dir$ano,dir)) {
      id.list <- find.file(id.list,dir$ano,dir)
    } 
  }
  id.list <- force.vec(id.list)
  id.list <- unique(id.list)
  sample.info <- data.frame(row.names=paste(id.list),grp=rep(1,times=length(id.list)))
  if(!is.null(grp)) {
    sample.info <- add.to.sample.info(sample.info,grp,"grp",overwrite=overwrite,verbose=verbose)
  }
  if(!is.null(dir)) {
    if(is.character(plate) & length(plate)==1) { plate <- find.file(plate,dir$ano,dir) }
    if(is.character(phenotype) & length(phenotype)==1) { phenotype <- find.file(phenotype,dir$ano,dir) }
    if(is.character(sex) & length(sex)==1) { sex <- find.file(sex,dir$ano,dir) }
    plate.look <- get.plate.info(dir,...) 
    if(!is.null(plate.look)) {
      sample.info <- add.to.sample.info(sample.info,plate.look[[1]],c("plate","well"),overwrite=overwrite,verbose=verbose)
      if(any(colnames(sample.info) %in% "plate")) { plate <- NULL }
    }
    if(pheno.default.loc %in% list.files(dir$ano)) {
      phenot <- read.table(file=cat.path(dir$ano,pheno.default.loc))
      sample.info <- add.to.sample.info(sample.info,phenot,"phenotype",overwrite=overwrite,verbose=verbose)
    } else {
      if(is.null(phenotype)) {
        warning(paste("phenotype not specified, nor was the file",pheno.default.loc,"found in dir$ano"))
      }
    }
    sample.info <- samp.update.qc.fail(sample.info,dir,verbose=verbose)
  } 
  if(!is.null(plate)) {
    sk <- F;  if(is.character(plate)) { if(plate=="") { warning("no plate.lookup.txt file found"); sk <- T } }
    if(!sk) { sample.info <- add.to.sample.info(sample.info,plate,c("plate","well"),overwrite=overwrite,verbose=verbose) }
  }
  if(!is.null(phenotype) & is.null(phenot)) {    
    sk <- F;  if(is.character(phenotype)) { if(phenotype=="") { warning("no phenotype.lookup.txt file found"); sk <- T } }
    sample.info <- add.to.sample.info(sample.info,phenotype,"phenotype",overwrite=overwrite,verbose=verbose)
  }
  if(!is.null(sex)) {
    sk <- F;  if(is.character(sex)) { if(sex=="") { warning("no sex.lookup.txt file found"); sk <- T } }
    if(!sk) { 
      sample.info <- add.to.sample.info(sample.info,sex,"sex",overwrite=overwrite,verbose=verbose) 
      sample.info$sex <- force.sex.codes(sample.info$sex)
    }
  }
  if(is.list(other)) {
    if(length(other)>0) {
      for(cc in 1:length(other)) {
        sample.info <- add.to.sample.info(sample.info,other[[cc]],names(other)[cc],overwrite=overwrite,verbose=verbose)
      }
    }
  } else {
    if(!is.null(other)) {
      sample.info <- add.to.sample.info(sample.info,other,overwrite=overwrite,verbose=verbose)
    }
  }
  sample.info <- validate.samp.info(sample.info,dir,verbose=verbose)
  return(sample.info)
}


add.to.sample.info <- function(sample.info,more.info,to.add=NULL,overwrite=F,verbose=F)
{
  # adds columns 'to.add' or else all valid columns from table 'more.info' to
  # the object sample.info, for matching rownames. autodetects id column in 
  # more.info based on the rownames of sample.info. can optionally overwrite existing columns
  if(verbose) { cat("\nAdding batch information to 'sample.info' table\n") }
  if(is.null(colnames(more.info))) 
  { 
    fake.names <- T # ie, indicating that no colnames present so using fake ones
    #colnames(more.info) <- paste("col",1:ncol(more.info)) 
  } else {
    fake.names <- F
  }
  more.info <- force.frame(more.info)
  if(verbose) { print(head(more.info)) }
  mat.fil.list <- find.id.col(more.info,ids=rownames(sample.info))
  if(verbose) { cat("",round(mat.fil.list$maxpc*100,1),
                    "% of ids in batch lookup table matched sample.info\n") }
  indx <- mat.fil.list$index
  
  if(length(narm(indx))>0) {
    whichoz <- c(1:ncol(more.info))
    # remove ID column and if overwrite=F, any columns already in sample.info
    if(!overwrite) {
      nn <- which(colnames(sample.info) %in% colnames(more.info))
    } else { nn <- NULL }
    mm <- (mat.fil.list$col); mm <- c(mm[mm!=0],nn[nn!=mm])
    if(length(mm)!=0) { whichoz <- whichoz[-mm] }
    # left with list of all valid columns in more.info
    to.try <- paste(colnames(more.info)[whichoz])
    to.use <- to.add[toupper(to.add) %in% toupper(to.try)] # select only requested ones
    if(length(to.use)==0) {
      to.use <- to.try # if no request colnames found, use all possible
    }
    ## decide what column names will be in sample.info
    if(fake.names & length(to.add)==1 & length(to.use)==1) {
      #only 1 var to be added, only 1 in dataframe but unnamed, so use to.add name
      to.name <- to.add
    } else {
      to.name <- to.use
    }
    for(jj in 1:length(to.use)) {
      # allow matching to version with any upper/lower-case configuration
      wHERE <- which(toupper(colnames(more.info)) %in% toupper(to.use[jj]))
      if(length(wHERE)>0)
      {  
        if(length(wHERE)>1) {
          # in case same variable exists in different case, if possible choose the exact match
          wHERE2 <- which(colnames(more.info) %in% to.use[jj])
          if(length(wHERE2)>0) { wHERE <- wHERE2[1] }
        }
        sample.info[[to.name[jj]]] <- more.info[[wHERE[1]]][indx]
      }
    }
    if(verbose) { cat(" added columns",to.name,"[from columns",
                      to.use,"in 'more.info'] to sample.info\n",sep=", ") }
  } else {
    stop("no subject ids matched ids in annotation lookup file\n")
  }
  return(sample.info)
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



samp.cr.summary <- function(bigMat,histo=F,print=F) {
  # calculate call rate using percentage of missing by column, in a bigmatrix
  must.use.package(c("bigmemory","biganalytics"))
  col.ms <- colna(bigMat)
  crs <- 1-(col.ms/nrow(bigMat))
  cat("mean sample call rate is now:",mean(crs,na.rm=T))
  if(histo) { textogram(crs) }
  call.rate.summary(cr.vec=crs,print=print)
  names(crs) <- colnames(bigMat)
  return(crs)
}


snp.cr.summary <- function(bigMat,histo=F,print=F,n.cores=1) {
  # calculate call rate using percentage of missing by column, in a bigmatrix
  must.use.package(c("multicore","bigmemory"))
  countNA <- function(x) { length(which(is.na(x))) }
  row.ms <- multi.fn.on.big(bigMat,1,FUN=countNA,dir=dir$big,by=200,n.cores=n.cores)
  row.ms <- unlist(row.ms)
  #row.ms <- apply(bigMat,1,countNA)
  crs <- 1-(row.ms/nrow(bigMat))
  cat("mean snp call rate is now:",mean(crs,na.rm=T))
  if(histo) { textogram(crs) }
  call.rate.summary(cr.vec=crs,print=print)
  names(crs) <- rownames(bigMat)
  return(crs)
}



list.rowsummary <- function(snpMatLst,mode="row",dir="",warn=T,n.cores=1)
{
  # performs 'row.summary' or 'col.summary' snpStats functions
  # on a list of snpMatrix objects
  fail <- F
  typz <- sapply(snpMatLst,is)[1,]
  if(all(typz=="SnpMatrix"))
  { HD <- F } else {
    if (all(typz=="character")) { HD <- T } else {
      warning("snpMatLst doesn't appear to be a list of SnpMatrix objects",
              "or a list of file locations, row.summary impossible")
      fail <- T
    } 
  }
  # this line allows this function to be either row.summary or col.summary (mode)
  if(mode!="col" & mode!="row") { mode <- "row" }
  my.summary <- switch(mode,row=row.summary,col=col.summary)
  if (!fail) {
    must.use.package("snpStats",T)
    rowsum.list <- vector("list",length(snpMatLst))
    if(!warn) { 
      #avoid alarming user given known issue with genotypes that are uniformly zero (empty) in column mode
      options(warn=-1) 
    }
    if(n.cores>1 & !HD) {
      must.use.package("multicore")
      cat(" processing elements in",min(n.cores,length(snpMatLst)),"parallel cores ...")
      rowsum.list <- mclapply(snpMatLst,my.summary)
    } else {
      cat(" processing element ")
      for (dd in 1:length(snpMatLst))
      {
        cat(dd,"..",sep="")
        if(HD) {
          fl.nm <- cat.path(dir$cr,snpMatLst[[dd]])
          obj.nm <- paste(load(fl.nm))
          #print(dim(get(obj.nm))) ; badness <- which(is.na(get(obj.nm)),arr.ind=T); print(head(badness)); print(dim(badness))
          rowsum.list[[dd]] <- my.summary(get(obj.nm))
        } else {
          rowsum.list[[dd]] <- my.summary(snpMatLst[[dd]])
          #print(dim(snpMatLst[[dd]]))
        }
      }
    }
    cat("done\n")
    if(!warn) { options(warn=0) }
    return(do.call("rbind",rowsum.list))
  } else {
    return(NULL)
  }
}


list.colsummary <- function(snpChrLst,mode="col",dir=dir,warn=F,n.cores=1)
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
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
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
  rownames(list.spec) <- c("Samples","SNPs")
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
    must.use.package("multicore"); snpChrLst <- mclapply(1:n.chr,one.chr,mc.cores=min(n.cores,4))
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
  chr.set <- chr.nums.of.info(snp.info); n.chr <- length(chr.set)
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



LRR.gc.correct <- function(dir,snp.info,bigLRR,pref="GC",write=F,add.means=T,n.cores=1)
{
  ## using GC%, run correction for GC-wave on a dataset
  load.all.libs()
  dir <- validate.dir.for(dir,c("big","gc"))
  origMat <- getBigMat(bigLRR,dir)
  cat("\nRunning GC correction, using LRR-dataset:\n")
  print.big.matrix(origMat,name="preCorrectedMat")
  # get filenames now to add to result later
  rN <- rownames(origMat); cN <- colnames(origMat)
  # run pca.correction using ectors (PCs) and alues from LRR.PCA
  if(is(snp.info)[1]=="RangedData") {
    gc.dat <- get.gc.markers(dir=dir, snp.info=snp.info, n.cores=n.cores)
    mySnpsGC <- gc.dat$gc[match(rownames(snp.info),rownames(gc.dat))]
    gc <- mySnpsGC[match(rN,rownames(snp.info))]
    cat(length(which(is.na(gc))),"/",length(rN)," markers missing GC data replaced with mean\n",sep="")
    gc[is.na(gc)] <- mean(gc,na.rm=T)
    gc_log <- log(gc+1)
  } else {
    stop("Error: snp.info invalid")
  }
  
  # create new matrix same size, ready for corrected values
  nR <- nrow(origMat); nC <- ncol(origMat)
  cat(" creating new file backed big.matrix to store gc-corrected data...")
  gcCorMat <- filebacked.big.matrix(nR,nC, backingfile=paste(pref,"Bck",sep=""),
                                    backingpath=dir$big, descriptorfile=paste(pref,"Descr",sep=""))
  cat("done\n")
  if(!is.filebacked(gcCorMat) | !is.filebacked(origMat)) {
    warning("at least one of the big.matrices is not filebacked, memory problems may be encountered")
  }
  # rows are subjects / samples
  col.sel <- 1:ncol(origMat)
  cat(" correcting by gcwave percentage, taking the LRR lm-residual for each SNP\n")
  jj <- proc.time()
  num.samps <- ncol(origMat); snpz <- 1:nrow(origMat)
  stepz <- round(seq(from=1,to=num.samps+1,by=100))
  if((tail(stepz,1)) != num.samps+1) { stepz <- c(stepz,num.samps+1) }
  split.to <- length(stepz)-1
  big.extras <- T # flush memory every 'n' iterations.
  flush.freq <- 5
  for (dd in 1:split.to)
  {
    x1 <- stepz[dd]; x2 <- stepz[dd+1]-1 #subset row selection
    # use of this 'sub.big.matrix' structure, stops the memory leak behaviour which spirals
    # the memory relating to 'origMat' out of control.
    next.cols <- sub.big.matrix(origMat, firstCol=x1, lastCol=x2, backingpath=dir$big )
    # next.rows is now a pointer to a matrix subset, must use 'as.matrix' to coerce to a regular R object 
    #    if(add.means) { sample.means <- rep(colmean(next.rows,na.rm=T),each=nrow(next.rows)) } else { sample.means <- 0 }
    gcCorMat[,x1:x2] <- gcCorrect(as.matrix(next.cols),gc) #+ sample.means    
    loop.tracker(dd,split.to)
    ## Every 'flush.freq' iterations, clean up the memory, remove the 
    ##  big.matrix object 'gcCorMat' and re-attach it 
    if(dd %% flush.freq == 0) {    
      fl.suc <- flush(gcCorMat) & flush(next.cols)
      if(!fl.suc) { cat("flush failed\n") } 
      gc()  # garbage collection
      if(big.extras) {
        RR <- describe(gcCorMat)
        rm(gcCorMat)
        gcCorMat <- attach.big.matrix(RR,path=dir$big)
      }
    }
    rm(next.cols) # remove the sub-matrix pointer each iteration or this memory builds up 
  }
  options(bigmemory.allow.dimnames=TRUE)
  rownames(gcCorMat) <- rN;  colnames(gcCorMat) <- cN 
  ll <- proc.time()
  cat(paste(" LRR GC-Correction took",round((ll-jj)[3]/3600,3),"hours\n"))
  flush(gcCorMat) # should allow names to take  
  cat("\nGC-corrected dataset produced:\n")
  print.big.matrix(gcCorMat,name="gcCorMat")
  
  mat.ref <- describe(gcCorMat)
  if(write) {
    big.fn <- "describeGCcorrect.RData"
    ofn <- cat.path(dir$big,big.fn)
    save(mat.ref,file=ofn)
    cat(paste("~wrote GC-corrected data description file to file:",ofn,"\n"))
    return(big.fn)
  } else {
    return(mat.ref)
  }
}


get.gc.markers <- function(dir, snp.info=NULL, snp.fn="snpdata.map", anot="map3", wndw=(10^6/2), ucsc="hg18", 
                           ret=c("gc","bio")[2], reset=F, n.cores=1,...)
{
  # for a list of SNPs, get the average GC% in the 1MB surrounding window
  # load annotation package from bioconductor
  dir <- validate.dir.for(dir,c("ano","gc"),warn=F)
  store.fn <- cat.path(dir$gc,"marker.gc.RData")
  must.use.package(c("BiocGenerics","genoset"),T)
  cat("\nGenerate GC percentage for SNPs\n")
  if(is.null(snp.info)) {
    snp.info <- make.snp.info(dir, snp.fn=snp.fn, anot=anot, verbose=F, ucsc=ucsc)
    cat(" made new snp.info object\n")
  }
  if(is(snp.info)[1]=="RangedData") {
    uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
    snp.info <- select.autosomes(snp.info)
    chr.set <- chr.nums.of.info(snp.info)
    cat(" using UCSC build",ucsc,"and chromosomes",paste(chr.set,collapse=","),"\n")
  } else { stop("Error: GC calculation failed due to invalid snp.info object") }
  hgFn <- paste("BSgenome.Hsapiens.UCSC.",ucsc,sep="")
  suppressWarnings(suppressMessages(must.use.package(c(hgFn),bioC=T)))
  # recalculate: chrLens,chrStarts for bioconductor annotation (could be different to other annotation)
  # these are just the lengths and genome starting positions for each autosome
  chrLs <- as.numeric(seqlengths(Hsapiens)[chr.set])
  # Need to create this info exactly according to the bioconductor annotation
  snp.info.B <- remove.boundary.snps.window(snp.info,window=wndw,chrLens=chrLs)
  ## find replace all chr.info with snp.info/equiv
  
  gotGC <- F
  if(!reset & file.exists(store.fn)) {
    snp.info.new <- snp.info.B
    obj.names <- paste(load(store.fn))
    if(("gc.dat" %in% obj.names) & ("snp.info.B" %in% obj.names)) {
      if(all(dim(snp.info.new)==dim(snp.info.B))) {
        if(all(rownames(snp.info.new)==rownames(snp.info.B))) {
          cat(" loaded marker-GC from file (current snp list seemed to exactly matched stored result)\n")
          gotGC <- T
        }
      }
    }
    if(!gotGC) { warning("existing GC file will be replaced as was not valid for this dataset") }
    snp.info.B <- snp.info.new 
    # put back our newer snp.info object into the b, else will always keep the old one
  }
  if(!gotGC) {
    cat("\n loading GC data for SNPs..."); st.time <- proc.time()[3]
    gc.dat <- suppressWarnings(calcGC(snp.info.B, expand = wndw, bsgenome = Hsapiens, return.bio=T, n.cores=n.cores))
    cat("took",round((proc.time()[3]-st.time)/60,1),"minutes\n")
  }
  misn <- which(gc.dat$gc==0); cat("",length(misn),"zero GC values set to missing (NA)\n")
  gc.dat$gc[misn] <- NA # set any zero values to NA (as zero = missing)
  # save now to save time later
  save(gc.dat,snp.info.B,file=store.fn)
  # choose from 3 possible formats to return GC data
  return(switch(ret,gc=round(gc.dat$gc,2),bio=gc.dat))
}


get.gc.human <- function(windowsize=10^6,ucsc="hg18",ret=c("bio","gc")[1], n.cores=1) {
  # get the mean GC% for every 1MB window of the human genome
  # can set UCSC build ('ucsc') to hg18 or hg19
  hgFn <- paste("BSgenome.Hsapiens.UCSC.",ucsc,sep="")
  suppressWarnings(suppressMessages(must.use.package(c("genoset","BiocGenerics",hgFn),bioC=T)))
  # get chromosome lengths from annotation (possibly slightly different to other sources)
  chrLens <- as.numeric(seqlengths(Hsapiens)[1:22])
  chrsz <- chrLens %/% windowsize
  windowsts <- list()
  for (cc in 1:22) {	
    windowsts[[cc]] <- (((1:(chrsz[cc]-1))-1)*windowsize)+1 
  }
  leno <- sapply(windowsts,length)
  coco <- rep(c(1:length(leno)),leno)
  stz <- as.integer(as.numeric(unlist(windowsts)))
  st.time <- proc.time()[3]
  cat("\nExtracting GC% for each",windowsize,"window of the human genome...")
  lData <- RangedData(ranges=IRanges(start=(stz+(windowsize/2)),width=1),
                      space=coco,universe=ucsc)
  lData <- toGenomeOrder(lData,strict=T)
  suppressWarnings(gc.dat <- calcGC(lData, expand = windowsize/2, bsgenome = Hsapiens,return.bio=T, n.cores=n.cores))
  #print(is(gc.dat)); print(length(gc.dat)); print(dim(gc.dat)); print(gc.dat)
  #if(ret=="bio") {
  #  gc.dat <- lData; gc.dat[["gc"]] <- gc.vec
  #  gc.dat$ranges <- (flank(gc.dat$ranges,width=windowsize/2,both=T)) # accurately reflect GC ranges (expand)
  #}
  cat("took",round((proc.time()[3]-st.time)/60,1),"minutes\n")
  #cat("NB: please ignore warnings about previous imports and 'number of items to replace' - this is normal\n")
  det.nm <- paste("package:",hgFn,sep="")
  detach(det.nm,character.only = TRUE) # un-load large annotation package.
  return(switch(ret,bio=gc.dat,gc=gc.dat$gc))
}


calcGC <- function(object, bsgenome, expand=1e6, return.bio=T, missing.as=NA, n.cores=1) {
  # taken from genoset package
  if (!requireNamespace("BSgenome",quietly=TRUE)) {
    stop("Failed to require BSgenome package.\n")
  }
  if (!requireNamespace("Biostrings",quietly=TRUE)) {
    stop("Failed to require Biostrings package.\n")
  }
  if(n.cores>1) { must.use.package("multicore") }
  calc.one.chr <- function(chr.name) {
    range = seq.int(chr.ind[chr.name, 1], chr.ind[chr.name, 2])
    seq = bsgenome[[chr.name]]
    v = suppressWarnings( Views(seq,  start=start[range], end=end[range]) )
    alf = alphabetFrequency(v, as.prob = TRUE)
    gc = rowSums(alf[,  c("G",  "C"), drop=FALSE])
  }
  chr.ind = chrIndices(object)
  rownames(chr.ind) = paste0("chr", rownames(chr.ind))
  start = start(object) - expand
  end = end(object) + expand
  gc.list = mclapply(rownames(chr.ind),calc.one.chr,mc.cores=n.cores)
  gc = do.call(c, gc.list)
  # decide what to do with missing values if they occur
  gc[is.na(gc)] <- missing.as
  if(return.bio) {
    gc.dat <- object; gc.dat[["gc"]] <- gc
    gc.dat$ranges <- (flank(gc.dat$ranges,width=expand,both=T)) # accurately reflect GC ranges (expand)
    return(gc.dat)
  } else {
    return(gc)
  }
}


get.all.snp.fails <- function(dir,verb=F)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl2"),warn=F)
  snps.to.cut <- NULL
  all.fls <- list.files(dir$excl2)
  if(length(all.fls)==0) { return(character(0)) } 
  for (cc in 1:length(all.fls))
  {
    pr.len <- length(snps.to.cut)
    next.fl <- readLines(cat.path(dir$excl2,all.fls[cc]))
    snps.to.cut <- unique(c(snps.to.cut,next.fl))
    if(verb) {  cat("",(length(snps.to.cut)-pr.len),
                    "Snps to remove, found in",all.fls[cc],"\n") }
  }
  return(snps.to.cut)
}


get.all.samp.fails <- function(dir,verb=F)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl"),warn=F)
  samps.to.cut <- NULL
  all.fls <- list.files(dir$excl)
  if(length(all.fls)==0) { return(character(0)) }
  for (cc in 1:length(all.fls))
  {
    pr.len <- length(samps.to.cut)
    next.fl <- readLines(cat.path(dir$excl,all.fls[cc]))
    samps.to.cut <- unique(c(samps.to.cut,next.fl))
    if(verb) {  cat("",(length(samps.to.cut)-pr.len),
                    "samples to remove (across all groups), found in",all.fls[cc],"\n") }
  }
  return(samps.to.cut)
}


plot.all.ranges <- function(cnv.ranges,DT=NULL,file="all.ranges.pdf",dir="",pc.flank=5,snp.info=NULL,
                            col1="black",scheme="mono",LRR=T,BAF=F,bafOverlay=F,hzOverlay=F,
                            PREPOSTPC=F,baf.file="big.baf",lrr.file="big.pcc",n.cores=1,n.pcs=NA,...) {
  must.use.package(c("multicore","bigmemory","biganalytics")); must.use.package("genoset",T)
  dir <- validate.dir.for(dir,c("big","res"))
  if(is(cnv.ranges)[1]!="RangedData") { warning("not a RangedData object"); return(NULL) }
  if(any(!c("id") %in% colnames(cnv.ranges))) { warning("cnv.ranges must contain id"); return(NULL) }
  idz <- cnv.ranges$id; stz <- start(cnv.ranges); enz <- end(cnv.ranges); chrz <- chr(cnv.ranges); wz <- width(cnv.ranges)
  if(!is.data.tracker(DT) & all(dir!="")) { DT <- read.data.tracker(dir,warn.only=T) }
  if((!is.character(lrr.file) | !is.character(baf.file))) { lrr.file <- baf.file <- "";
                                                            warning("lrr.file/baf.file must be character locations, defaults will be used from tracker") }
  XX <- get.tracker.for.plot(DT=DT,dir=dir,PREPOSTPC=PREPOSTPC,baf.file=baf.file,lrr.file=lrr.file,
                             LRR=LRR,BAF=BAF,samples=idz,n.pcs=n.pcs)
  if(bafOverlay) { pfb.rng <- get.pfb.ranges(dir) } else { pfb.rng <- NULL }
  if(hzOverlay) { hz.rng <- get.snpqc.ranges(dir,col="het") } else { hz.rng <- NULL }
  if(all(c("lrr.file","baf.file","lrr.file.raw") %in% names(XX))) {
    lrr.file <- XX$lrr.file; baf.file <- XX$baf.file; lrr.file.raw <- XX$lrr.file.raw
    if(!PREPOSTPC) {
      # load big objects to save loading each time for many ranges
      if(LRR) { lrr.file <- getBigMat(lrr.file,dir$big) }
      if(BAF) { baf.file <- getBigMat(baf.file,dir$big) }
    } else {
      warning("PREPOSTPC==T is slow for multiple ranges as data must be reloaded each iteration")
    }
  } else {
    warning("lookup of baf.file and/or lrr.file file(s) failed - will try to proceed but failure is likely")
  }
  ofn <- cat.path(dir$res,file)
  if(LRR & !PREPOSTPC) {
    excl.miss1 <- which(!cnv.ranges$id %in% colnames(lrr.file))
    #print(head(cnv.ranges$id)); print(head(colnames(lrr.file)))
    if(length(excl.miss1)>0) { warning(length(excl.miss1)," samples excluded as not in LRR matrix"); cnv.ranges <- cnv.ranges[-excl.miss1,] }
  }
  if(BAF & !PREPOSTPC) {
    excl.miss2 <- which(!cnv.ranges$id %in% colnames(baf.file))
    #print(head(cnv.ranges$id)); print(head(colnames(baf.file)))
    if(length(excl.miss2)>0) { warning(length(excl.miss2)," samples excluded as not in BAF matrix"); cnv.ranges <- cnv.ranges[-excl.miss2,] }
  }
  n.to.plot <- nrow(cnv.ranges)
  if(!is(snp.info)[1]=="RangedData") { snp.info <- read.snp.info(dir) }
  if(!is.na(scheme) & ("color" %in% colnames(snp.info))) { snp.info <- snp.info[,-match("color",colnames(snp.info))] } 
  if(!all(c("QCfail","gindx") %in% colnames(snp.info))) {
    snp.info <- add.gindx.to.Ranges(snp.info)
    snp.info <- snp.update.qc.fail(snp.info,dir=dir)
  }
  cL <- get.chr.lens(dir)
  if(n.to.plot<1) { return(NULL) } else { cat(n.to.plot,"CNVs to plot:\n")}
  #pdf(ofn)
  if(n.cores>1 & n.to.plot>1) {
    cc <- 1; nc <- 0; nutin <- nufin <- list()
    while(cc <=n.to.plot) {
      nc <- nc + 1
      poz <- c( max(0,(stz[cc]-(pc.flank*wz[cc]))), min((enz[cc]+(pc.flank*wz[cc])),cL[as.numeric(chrz[cc])]) ) #; print(poz)
      nufin[[nc]] <- multicore::parallel(suppressWarnings({
        cnv.plot(DT=DT,cnvPlotFileName=paste(file,cc,sep="."),samples=idz[cc],Chr=chrz[cc],Pos=poz,Cnv=c(stz[cc],enz[cc]),dir=dir,snp.info=snp.info,
               col1=col1,scheme=scheme,lrr.file=lrr.file,baf.file=baf.file,PREPOSTPC=PREPOSTPC,n.pcs=n.pcs,
               LRR=LRR,BAF=BAF,bafOverlay=bafOverlay,hzOverlay=hzOverlay,pfb.rng=pfb.rng,hz.rng=hz.rng,...) }))
      cc <- cc + 1
      if(nc==n.cores | cc==n.to.plot) { nutin <- collect(nufin); nc <- 0 }
    }
  } else {
    pdf(ofn)
    for(cc in 1:n.to.plot) {
      poz <- c( max(0,(stz[cc]-(pc.flank*wz[cc]))), min((enz[cc]+(pc.flank*wz[cc])),cL[as.numeric(chrz[cc])]) ) #; print(poz)
      suppressWarnings({
        cnv.plot(DT=DT,samples=idz[cc],Chr=chrz[cc],Pos=poz,Cnv=c(stz[cc],enz[cc]),dir=dir,snp.info=snp.info,
               col1=col1,scheme=scheme,lrr.file=lrr.file,baf.file=baf.file,PREPOSTPC=PREPOSTPC,n.pcs=n.pcs,
               LRR=LRR,BAF=BAF,bafOverlay=bafOverlay,hzOverlay=hzOverlay,pfb.rng=pfb.rng,hz.rng=hz.rng,...) })
    }
    dev.off()
  }
  #dev.off()
  cat("\n~wrote file:",ofn,"\n")
}

#cnv.plot(DT,"5299865121_R01C01",Chr=1,Pos=c(242,243)*10^4)

force.chr.pos <- function(Pos,Chr,snp.info) {
  if(length(Pos)==2) {
    maxln <- end(tail(snp.info[paste(Chr)],1)) # force start and end to be within 1:chr.len
    mbs <- min(max(1,Pos[1]),(maxln-1)); mbe <- min(max(2,Pos[2]),maxln)
    return(c(mbs,mbe))
  } else {
    Pos <- NA; warning("Pos needs to be numeric length 2, min, max")
  }
  return(Pos)
}


unique.in.range <- function(ranged,chr,pos,full.overlap=F, unit=c("b","kb","mb","gb"), rmv.dup=T) {
  # select everything in a ranges format that is in a window of a chromosome; 
  if(length(pos)>2 | !is.numeric(pos)) { warning("pos should be a start and end numeric range"); return(NULL) }
  if(length(pos)==1) { pos <- rep(pos,2) }
  if(length(chr)>1 | is.na(as.numeric(chr))) { warning("chr should be a single number"); return(NULL) }
  if(!any(is(ranged)[1] %in% c("RangedData","IRanges","GRanges","RangesList"))) { warning("'ranged' should be a RangedData type or similar"); return(NULL) }
  unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  # get set of genes in a position range for a chromosome
  chr.genez <- ranged[paste(chr)]; 
  if(full.overlap) {
    ranged <- chr.genez[which(start(chr.genez)>min(pos) & end(chr.genez)<max(pos)),]
  } else {
    # any overlap
    ranged <- chr.genez[which(end(chr.genez)>min(pos) & start(chr.genez)<max(pos)),]
  }
  if(rmv.dup) {
    # remove duplicate genes/exons
    ranged <- ranged[!(duplicated(start(ranged)) & duplicated(end(ranged)) & duplicated(width(ranged))),]
  }
  return(ranged)
}


plot.gene.annot <- function(gs=NULL, chr=1, pos=NA, x.scl=10^6, y.ofs=0, width=1, txt=T, chr.pos.offset=0,
                            ucsc="hg18", dir="", box.col="green", txt.col="black", join.col="red", ...)
{
  # add gene annotation to existing plot with y range 'ylim',
  #  in range 'sect' using 'parse.annot.file' formmated
  #  annotation data read in from dataframe 'DF'. step is a scaling
  #  factor that will adjust the gene position units based on the 
  #  the step size of the rest of the plot.
  # more args to 'rect' ...
  dir <- validate.dir.for(dir,"ano")
  if(is(gs)[1]!="RangedData") { gs <- get.gene.annot(dir=dir,ucsc=ucsc) }
  if(!"gene" %in% colnames(gs)) { warning("didn't find 'gene' column in annotation") ; return(NULL) }
  Col <- c("green", "darkgreen")
  if(all(is.na(pos))) { pos <- c(1,Inf) } 
  # get set of genes in range of the graph section + remove duplicate genes/exons
  rng.genez <- unique.in.range(gs,chr,pos,full.overlap=F, rmv.dup=T)
  if(nrow(rng.genez)<1) { warning("no genes found in range") ; return(NULL) }
  old.cc <- 1 ; tp <- 2
  # set vertical alignments for annotation
  old.auto <- F
  y.cent <- y.ofs
  y.bot <- y.ofs-(width/2)
  y.top <- y.ofs+(width/2)
  # text alignment
  tps <- y.bot + c(.18,.29,.46,.64,.82)[c(1,3,5)]*width
  # x position with scaling (e.g, Mb units = 10^6)
  cnrlo <- start(rng.genez)/x.scl+chr.pos.offset
  cnrhi <- end(rng.genez)/x.scl+chr.pos.offset
  gnnm <- (rng.genez$gene)
  n.genes <- length(cnrlo)
  txt.cex <- .75; if(n.genes>10) { txt.cex <- .5 } ; if(n.genes>100) { txt.cex <- .35 } # more genes = smaller labels
  for (cc in 1:n.genes) {
    # draw rectangle and label for each gene
    rect(cnrlo[cc],y.top,cnrhi[cc], y.bot,border=box.col,...)
  }
  for (cc in 1:n.genes) {
    if (gnnm[old.cc]==gnnm[cc] & cc!=1)
    {
      link <- c(cnrhi[old.cc],cnrlo[cc])
      if(link[1]<link[2]) {
        lines(link,y=rep(y.cent,2),lwd=1,col=join.col,lty="dotted")
      }
    } else {
      if(txt) {
        if(cnrlo[cc] < min(pos/x.scl)) { 
          txt.x <- mean(c(min(pos/x.scl),min(cnrhi[cc],max(pos/x.scl))),na.rm=T)  
        } else { txt.x <- cnrlo[cc] }
        text(txt.x,tps[tp],gnnm[cc],col=txt.col,cex=txt.cex,pos=4,las=2,offset=0)
      }
    }
    old.cc <- cc  
    tp <- tp+1; if(tp==4) {tp <- 1}; #if(n.genes <=5) { tp <- 2 }
  }
  return(rng.genez)
}


draw.cnv.bounds <- function(cnv,chr.offset=NA,pos=NA,cnv.lty="dashed",cnv.col="orange",plot.scl=10^6) {
  # do abline at start and end of all CNV positions
  done <- F
  if(is(cnv)[1]=="RangedData") {
    # highlight CNVs in rangedData object
    if(all(is.na(chr)) & all(is.na(pos))) {
      cnv <- unique.in.range(cnv,chr,pos,full.overlap=F, rmv.dup=T)
    }
    if(nrow(cnv)<1) { warning("no CNVs in range"); return(NULL) }
    if("cn" %in% colnames(cnv)) {
      # plot zones
      e.lrr <- c(-1.65,-1,0,.586,1)
      for (jj in 1:nrow(cnv)) { 
        mid <- e.lrr[(cnv[jj,]$cn+1)]
        x1 <- (start(cnv)[jj]+chr.offset)/plot.scl
        x2 <- (end(cnv)[jj]+chr.offset)/plot.scl
        y1 <- (mid-.15); y2 <- (mid+.15)
        #     print(paste(x1,y1,x2,y2,sep="/"))
        rect(x1,y1,x2,y2,lty=cnv.lty,border=cnv.col) 
      }
      done <- T
    } else { 
      cnv <- cbind(start(cnv),end(cnv))
    }
  }
  if(!is.null(dim(cnv)) & !done) {
    if(ncol(cnv)==2) {
      # read cnvs s/es from columns
      if(nrow(cnv)<=5) {
        for (jj in 1:nrow(cnv)) { abline(v=(cnv[jj,]+chr.offset)/plot.scl,lty=cnv.lty,col=cnv.col) }
      } else {
        for (jj in 1:nrow(cnv)) { lines(x=(cnv[jj,]+chr.offset)/plot.scl,y=-1.5,lty="solid",col=cnv.col,lwd=5) }
      }
    } else {
      warning("'cnv' should have 2 columns - start and end position"); return(NULL)
    }
  } else {
    if(length(cnv)==2 & is.numeric(cnv)) {
      # plot a single CNV
      #print(cnv/plot.scl)
      abline(v=(cnv+chr.offset)/plot.scl,lty=cnv.lty,col=cnv.col)
    } else {
      warning("invalid 'Cnv' specification")
    }
  }
  return(cnv)
}


get.tracker.for.plot <- function(DT=NULL,dir="",PREPOSTPC=F,baf.file="",lrr.file="",
                                 LRR=T,BAF=F,samples=NULL,n.pcs=NA) {
  ## hidden function just for cnv.plot and plot.all.ranges
  lrr.file.raw <- ""
  if(is.data.tracker(DT)) {
    # if running with a DT, datatracker object then many datasources can be sourced automatically
    if(is.ch(getSlot(DT,"dir")) & !is.ch(dir)) { dir <- getSlot(DT,"dir") }
    DTs.required <- c("snp.info","big.lrr","big.baf")
    is.valid <- req.valid(DT,DTs.required)
    if(LRR) {
      if(PREPOSTPC) {
        # include comparison of pre vs post QC LRR at the locus, so need to import pre-pc too
        DTs.required <- c("big.qc","big.pcc")
        if(req.valid(DT,DTs.required)) {
          lrr.file <- "big.pcc"
          lrr.file.raw <- "big.qc"
          lrr.file.raw <- getSlot(DT,lrr.file.raw)
        } else {
          warning("couldn't make pre/post comparison as required data not in 'DT' tracker"); PREPOSTPC <- F
        }
      }
    } else { lrr.file <- NULL }
    if(is.valid) {
      if(is.null(dir)) { dir <- getSlot(DT,"dir") }; bdf <- F
      if(BAF) {
        if(is.null(baf.file)) { bdf <- T }; if(all(baf.file=="")) { bdf <- T }
        if(!is.file(baf.file,dir$big,dir)) { bdf <- T }
        if(bdf) { baf.file <- getSlot(DT,"big.baf") } # basically unless it's a real file, ignore 'baf.file' entered
      } else { baf.file <- "" }
      if(LRR) {
        if(lrr.file %in% c("big.lrr","big.filt","big.qc","big.pcc")) {
          if((lrr.file %in% c("big.lrr","big.filt")) & (getSlot(DT,"ngrps")>1)) {
            ## if it's one of these two and there are multiple groups, need to work out which file(s) data will be in for sample(s)
            si <- getSlot(DT,"sample.info",ret.obj=T); if(is.list(si)) { si <- si[[1]] }
            fail.this <- T
            if(is.data.frame(si)) {
              if(all(samples %in% rownames(si))) {
                grpz <- unique(si[paste(samples),"grp"]) # lookup which grp sample(s) is/are in
                if(length(grpz)==1) {
                  fail.this <- F
                  lrr.file <- getSlot(DT,lrr.file,grps=grpz,n.pcs=n.pcs); if(is.list(lrr.file)) { lrr.file <- lrr.file[[1]] } 
                } else { wrn <- "more than 1 grp in sample list, but data for grps is split across files" }
              } else { wrn <- "'grp' column not in sample.info" }
            } else { wrn <- "sample.info not valid" }
            if(fail.this) {
              warning("tried to use big matrix when ",wrn)
              return()
            }
          } else {
            lrr.file <- getSlot(DT,lrr.file,n.pcs=n.pcs) # one of the ones that should always work
          }
        } else {
          # non of the 'DT' abbreviations were used so assume it refers to an actual filename, but if not existing, try to use big.pcc
          if(!file.exists(lrr.file) | all(lrr.file=="")) { 
            if(req.valid(DT,"big.pcc")) {
              if(!all(lrr.file=="")) {
                warning("file 'lrr.file' did not exist, defaulting to PC-corrected LRR matrix") 
              }
              lrr.file <- getSlot(DT,"big.pcc",n.pcs=n.pcs)
            } else {
              stop("file 'lrr.file' did not exist and didn't find a PC corrected matrix to default to")
            }
          }
        }
      } else {
        lrr.file <- ""
      }
    } else {
      # probably gonna fail as it doesn't seem like the required files are there
      warning("datatracker object DT didn't contain files needed for plotting, will attempt to continue using remaining parameters")
    }
  }
  return(list(lrr.file=lrr.file[[1]],baf.file=baf.file[[1]],lrr.file.raw=lrr.file.raw[[1]]))
}


get.pfb.ranges <- function(dir) {
  must.use.package("genoset",T)
  dir <- validate.dir.for(dir,"cnv")
  pfb.fn <- cat.path(dir$cnv,"BAF.pfb") # bafOverlay
  if(file.exists(pfb.fn)) { pfb <- reader(pfb.fn,quiet=T) } else { warning("BAF.pfb file not found"); return(NULL) }
  if("PFB" %in% colnames(pfb)) {
    pfbData <- RangedData(ranges=IRanges(start=pfb$Position,end=pfb$Position,names=rownames(pfb)),
                          space=pfb$Chr,PFB=pfb$PFB)
    pfbData <- toGenomeOrder(pfbData,strict=T)
    pfbData <- add.gindx.to.Ranges(pfbData,ucsc=ucsc,absolute=T,label="gindx")     
    pfbData$PFB[is.na(pfbData$PFB)] <- 0
  }
  return(pfbData)
}


get.snpqc.ranges <- function(dir,col="het",def=NA) {
  must.use.package("genoset",T)
  dir <- validate.dir.for(dir,"cr")
  hz.fn <- cat.path(dir$cr,"snpqc.txt") # hzOverlay [heterozygosity]
  if(file.exists(hz.fn)) { hz <- reader(hz.fn,quiet=T) } else { warning("snpqc.txt file not found"); return(NULL) }
  if(col %in% colnames(hz)) {
    hzData <- RangedData(ranges=IRanges(start=hz$start,end=hz$end,names=hz$names),
                         space=hz$space)
    hzData <- toGenomeOrder(hzData,strict=T)
    hzData[[col]] <- hz[[col]]
    hzData <- add.gindx.to.Ranges(hzData,ucsc=ucsc,absolute=T,label="gindx")   
    hzData[[col]][is.na(hzData[[col]])] <- def
  }
  return(hzData)
}


cnv.plot <- function(dir="",samples="",LRR=T,BAF=F,PREPOSTPC=F,n.pcs=NA,
                     chr.expand=0,Chr=1:22,Pos=NA,Cnv=Pos,scl=10^6,
                     tag.cnvs=F,medSmooth=F,med.rat=10,gcOverlay=F,geneOverlay=F,exons=F,
                     bafOverlay=F,pfb.rng=NULL,hzOverlay=F,hz.rng=NULL,
                     cnvPlotFileName="",show.fail=F,
                     lrr.file="LRRFiltSortdescrFile",baf.file="BAFFiltSortdescrFile",
                     snp.info=NULL,scheme=NA,col1="lightblue",col2="lightgreen", 
                     cnv.col="orange", cnv.lty="dashed", gc.col="blue",ov.lty="dotted",
                     hz.col="navy",pfb.col="brown",gene.col="green",
                     c.xlim=NULL,c.ylim=NULL,cust.sub="",ucsc="hg18",DT=NULL) {
  customLims <- rngOn <- zoom <- F # change this later if valid limits received
  must.use.package(c("bigmemory","biganalytics")); must.use.package("genoset",T)
  #lrr.file.raw <- ""
  if(!LRR & !BAF) { warning("selected no LRR and no BAF so no plot(s) produced"); return(NULL) }
  if(!is.character(samples)) { warning("samples should be a list (character) of sample ids"); return(NULL) }
  if(!is.data.tracker(DT)) { DT <- read.data.tracker(dir,warn.only=T) }
  if(PREPOSTPC & show.fail) { warning("setting 'PREPOSTPC'=T, will prevent 'show.fail'=T from showing QC failing points in LRR") }
  dir <- validate.dir.for(dir,c("big","res","gc"))
  # gc overlay file
  gc.fn <- cat.path(dir$gc,"GCAv6.RData") # gcOverlay
  if((!is.character(lrr.file) | !is.character(baf.file))) {
    if(PREPOSTPC) {
      warning("PREPOSTPC==TRUE overrides custom object entry of lrr.file/baf.file; defaults for these files will be used")
      skip.locs <- F
    } else {
      if(is.big.matrix(lrr.file) | is.matrix(lrr.file) | is.data.frame(lrr.file)) { lrr.obj <- T } else { lrr.obj <- F }
      if(is.big.matrix(baf.file) | is.matrix(baf.file) | is.data.frame(baf.file)) { baf.obj <- T } else { baf.obj <- F }
      if(!(lrr.obj & baf.obj) & (lrr.obj | baf.obj) & (LRR & BAF)) {
        warning("only 1 of lrr.file/baf.file entered as object rather than location - object ignored and defaults will be used")
        skip.locs <- F
      } else {
        skip.locs <- baf.obj | lrr.obj
      }
    }
  } else { skip.locs <- F }
  if(!skip.locs) {
    XX <- get.tracker.for.plot(DT=DT,dir=dir,PREPOSTPC=PREPOSTPC,baf.file=baf.file,lrr.file=lrr.file,
                               LRR=LRR, BAF=BAF,samples=samples,n.pcs=n.pcs)
    if(all(c("lrr.file","baf.file","lrr.file.raw") %in% names(XX))) {
      lrr.file <- XX$lrr.file; baf.file <- XX$baf.file; lrr.file.raw <- XX$lrr.file.raw
    } else {
      warning("lookup of baf.file and/or lrr.file file(s) failed - will try to proceed but failure is likely")
    }
  } 
  ## over-ride settings incompatible with viewing a partial chromosome if ZOOM is on:
  if(is(snp.info)[1]!="RangedData") { snp.info <- read.snp.info(dir) }
  snp.info <- toGenomeOrder(snp.info,strict=T)
  chr.set <- chr.nums.of.info(snp.info); 
  if(any(!paste(Chr) %in% paste(chr.set))) { Chr <- Chr[paste(Chr) %in% paste(chr.set)] }
  if(length(Chr)<length(chr.set)) { rngOn <- T }
  if(length(Chr)==1 | length(chr.set)==1){
    if(all(!is.na(Pos))){
      Pos <- force.chr.pos(Pos,Chr,snp.info)
      if(length(Pos)==2) { zoom <- T; mbs <- Pos[1]; mbe <- Pos[2] }
    }
  } else { 
    if(all(!is.na(Pos))){ warning("'Pos' ignored as more than one chromosome selected") }
  }
  if(zoom) 
  { 
    rngOn <- F #### JUST DID THIS HERE
    #rngOn <- T ; Chr <- Chr[1] 
    #med.rat <- max(10,med.rat)
    # make black defaul when only plotting 1 chromosome
    if(all(is.na(scheme)) & col1=="lightblue" & col2=="lightgreen") { scheme <- "mono"; col1 <- "black" }
    if(medSmooth) { medSmooth <- F; warning("medSmooth is set to false when Pos is in use") }
    if(chr.expand!=0) { 
      warning("when using a window 'Pos' on 'Chr', chr.expand will be set to zero")
      chr.expand <- 0 
    }
    rng.mb <- range(c(mbs,mbe))/10^6
    if(is.null(c.xlim)) { c.xlim <- rng.mb }
    #cat("using zoomed plot of Chr",Chr,"from",mbs,"mB to",mbe,"mB\n")
  } else { 
    # if not zooming into a range within a single chromosome
    rng.mb <- NA
  }  
  if(LRR) { 
    bigMat2 <- getBigMat(lrr.file,dir$big) ;
    if(PREPOSTPC) { 
      bigMat1 <- getBigMat(lrr.file.raw,dir$big)
      #    print.big.matrix(bigMat1,dir=dir$big,name="LRR_Matrix")
      not.in.lrr2 <- !(samples %in% colnames(bigMat1))
      if(any(not.in.lrr2))
      {
        samples <- samples[!not.in.lrr2]
        lll <- length(which(not.in.lrr2))
        if(length(samples)<1) {
          warning("sample(s) not found in LRR big matrix datafile (pre/post)")
          return(NULL)
        } else {
          warning(paste("Excluding",lll,"samples not in the pre-PCA LRR matrix"))
        }
      }
    }
    #    print.big.matrix(bigMat2,dir=dir$big,name="LRR_Matrix")
    not.in.lrr <- !(samples %in% colnames(bigMat2))
    if(any(not.in.lrr))
    {
      samples <- samples[!not.in.lrr]
      if(length(samples)<1) {
        warning("sample(s) not found in LRR big matrix datafile")
        return(NULL)
      } else {
        warning(paste("Excluding",length(which(not.in.lrr)),"samples not in the LRR matrix"))
      }
    }
  }
  if(BAF) { 
    bigMat3 <- getBigMat(baf.file,dir$big) 
    #    print.big.matrix(bigMat3,dir=dir$big,name="BAF_Matrix")
    not.in.baf <- !(samples %in% colnames(bigMat3))
    if(any(not.in.baf))
    {
      samples <- samples[!not.in.baf]
      if(length(samples)<1) {
        warning("sample(s) not found in big matrix BAF datafile - setting BAF=F")
        BAF <- F
      } else {
        warning(paste("Excluding",length(which(not.in.baf)),"samples not in the BAF matrix"))
      }
    }
  }
  if(!is.null(c.ylim)) {
    L1 <- as.numeric(c.ylim[1]); 	L2 <- as.numeric(c.ylim[2])
    if(!is.na(L1) & !is.na(L2)) { if (L1>-20 & L2<20) { customLims <- T }  }
  }
  auto.lab <- LRR | BAF  #LRR & !BAF
  numplotsin1 <- sum(as.numeric(c(LRR,BAF)))
  chrN <- vector("list",length(samples)); chrLab <- character(length(samples))
  names(chrN) <- names(chrLab) <- samples
  for (jj in 1:length(samples)) { chrN[[samples[jj]]] <- Chr  } #[1] 
  chrLab <- rep(paste("Chromosome",Chr[1]),times=length(chrLab))
  if(rngOn){
    chrLab[samples] <- ""
    # titl <- c(paste("Sample",samples[cc]),paste(chrLab[cc]))
  }
  c.post.n <- c.pre.n <- chr.expand 
  
  ## plot all samples in samples
  if(all(cnvPlotFileName!="")) { cnvPlotFileName <- paste(dir$res,cnvPlotFileName,sep="") }
  # set default to select all chromosomes for plotting
  chr.select <- 1:length(chr.set) # should be overwritten after first plot
  # if smoothing option in use then adjust y axis limits accordingly
  if (medSmooth) { limz <- c(-1,1) } else { limz <- c(-2,2) }
  blimz <- c(0,1)
  if (gcOverlay & medSmooth) { limz <- limz + c(0, 0.2) }
  if (customLims) { limz <- c(L1,L2) }
  xl <- paste(c("Genome","Chromosome")[1+as.numeric((length(Chr)+chr.expand)==1)],"Position (Megabases)")
  if(all(cnvPlotFileName=="") & length(samples)>1) { cnvPlotFileName <- "CNV_BAF_LRR.pdf" }
  if(all(cnvPlotFileName!="")) { 
    ofn <- cat.path(dir$res,cnvPlotFileName)
    pdf(ofn); plot.to.file <- T
  } else { plot.to.file <- F }
  ## loop through each sample plotting the desired plots
  if(length(samples)>1) { cat(paste("plotting",length(samples),"samples...\n")) }
  for (cc in 1:length(samples))
  {
    titl <- c(paste("Sample",samples[cc]))
    if(cust.sub!="") { titl <- c(titl,cust.sub) }
    par(mfrow=c(numplotsin1,1))
    # start plotting 
    if(LRR) {
      if (BAF) { xl1 <- "" ; xt <- "n" } else { xl1 <- xl; xt <- "s" }
      par(mar=list(c(5, 4, 10, 2),c(0.2, 4, 5, 2))[[numplotsin1]])
      # full LRR data for chromosomes (full colour)
      #print(c.xlim)
      #print(samples[cc])
      
      #if(PREPOSTPC) { col1a <- col2a <- "white" ; schemea = "mono" } else
       # points.only back to T
      if(!PREPOSTPC)  { 
        col1a <- col1; col2a <- col2; schemea = scheme 
        #plot(c(1,1),c(1,1),ylim=limz, xlim=c.xlim, main = titl,  xlab = xl1, ylab="Log-R Ratio", xaxt=xt,col="white")
        out.list <- col.plot.lrr(samples[cc], bigMat2, dir=dir, centre.chr=chrN[[samples[cc]]], 
                               scheme=schemea, col1=col1a, col2=col2a, snp.info=snp.info, plotAdj=rngOn, points.only=F, 
                               show.fail=show.fail, add.leg=show.fail,scl=scl,
                               c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=medSmooth, ucsc=ucsc, trunc=limz,
                               ratio=med.rat, ylim=limz, xlim=c.xlim, rng.mb=rng.mb, main = titl,  xlab = xl1, ylab="Log-R Ratio", xaxt=xt) 
      }
      if(PREPOSTPC) {
        ## put down pre-pc LRR data first for comparison at locus
        out.list <- col.plot.lrr(samples[cc], bigMat1, dir=dir, centre.chr=chrN[[samples[cc]]], 
                                  scheme=scheme, col1="grey", col2="darkgrey", snp.info=snp.info, plotAdj=rngOn,
                                  points.only=F, show.fail=show.fail, add.leg=F,scl=scl,
                                  c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=medSmooth, ucsc=ucsc,  trunc=limz,
                                  ratio=med.rat, ylim=limz, xlim=c.xlim, rng.mb=rng.mb, main = titl,  xlab = xl1, ylab="Log-R Ratio", xaxt=xt)
        null.list <- col.plot.lrr(samples[cc], bigMat2, dir=dir, centre.chr=chrN[[samples[cc]]], 
                                  scheme=scheme, col1=col1, col2=col2, snp.info=snp.info, plotAdj=rngOn, points.only=T,
                                  show.fail=show.fail, add.leg=F,scl=scl,
                                  c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=medSmooth, ucsc=ucsc, trunc=limz,
                                  ratio=med.rat, rng.mb=rng.mb)
        legend("bottomleft",legend="Pre-PC LRR",pch=19,col="grey",cex=0.8,bty="n")
        legend("bottomright",legend="Post-PC LRR",pch=19,col=col1,cex=0.8,bty="n")
      }
      x.coords <- out.list$x.coords; chr.select <- out.list$chr.select  
      # legend("topleft",legend=c("Raw LRR data [coloured by chromosome]"),cex=.7,
      #         pch=c(21),col=c("grey"),bty="n")
      chrn <- chrN[[samples[cc]]]; 
      if(length(chrn)>1) { 
        chrn <- "number"
        # add chr labels with target coloured red
        null.result <- add.chr.top.axis(chr.select, badChrN=NULL, x.coords, nC=length(chr.set), sub=F)
      }
      if(auto.lab) { mtext(paste("Chromosome",chrn),side=3,line=0.5,cex=.9) }
      if(geneOverlay) {
        if(exons) { dat <- get.exon.annot(dir) } else { dat <- get.gene.annot(dir) }
        if(zoom) {
          plot.gene.annot(gs=dat, chr=Chr[1], pos=Pos, x.scl=scl, y.ofs=(as.numeric(medSmooth)*.5)-1, width=.5, txt=T,
                          ucsc=ucsc, box.col=gene.col, txt.col="black", join.col="red", dir=dir)
        } else {
          warning("no gene overlay option when plotting more than 1 chromosome")
        }
      }
      if(tag.cnvs ) {
        chr.offset <- 0
        if(is(Cnv)[1]=="RangedData") { if("id" %in% colnames(Cnv)) { Cnv <- Cnv[Cnv$id %in% samples,] } }
        if((length(Chr)+chr.expand)>1) {
          if(Chr[1]>1 & (length(Chr)>21 | chr.expand==0)) { 
            chrLens <- get.chr.lens(dir,ucsc=ucsc)[chr.set]
            chrStarts <- c(0,cumsum((chrLens)))
            chr.offset <- chrStarts[which(chr.set %in% Chr[1])]
            if(length(Chr)>1) { warning("assume Cnv coordinates refer to ",Chr[1]," given 'Chr' has length>1") }
          }
        }
        uu <- draw.cnv.bounds(cnv=Cnv,chr.offset=chr.offset,pos=Pos,cnv.lty=cnv.lty,cnv.col=cnv.col,plot.scl=scl)
        legend("bottom",legend="CNV",lty=cnv.lty,col=cnv.col,bty="n",cex=.8)
      }
    }
    if(zoom & gcOverlay) { if((diff(Pos)/scl)<3) {  gcOverlay <- F ; warning("gcOverlay ignored as resolution in window is too narrow") } }
    if (!gcOverlay & LRR) {
      # load average data for 10^6 windows prepared using 'extractGCwindowsFromSangerTxt.R'
      abline(h=0,col=gc.col,lty=ov.lty)
      legend("topleft",legend="normal copy baseline (0)",lwd=1,bty="n",col=gc.col,lty=ov.lty,cex=0.8)
    }
    if (gcOverlay & file.exists(gc.fn) & LRR) {
      # load average data for 10^6 windows prepared using 'extractGCwindowsFromSangerTxt.R'
      gc.avs6 <- reader(gc.fn,quiet=T); if(median(gc.avs6$gc)>2) { gc.avs6$gc <- gc.avs6$gc/100 }
      gc.avs6$gc[gc.avs6$gc<.1] <- NA  # erase v-low/zero values as likely just due to gaps in annotation
      gcer <- add.gindx.to.Ranges(gc.avs6)
      ofs <- out.list$offset; # print(ofs)
      xes <- (gcer$gindx/scl)-ofs
      lines(xes,gcer$gc,col=gc.col,lty=ov.lty)
      legend("topleft",legend="GC Percentage",lwd=1,bty="n",col=gc.col,lty=ov.lty,cex=0.8)
    }
    if(BAF) {
      if (LRR) { titl <- "" }
      par(mar=list(c(5, 4, 10, 2),c(5, 4, 0.2, 2))[[numplotsin1]])
      if (F) { titl <- "" ; #if(LRR) instead of if(F)
                 # drawing this as a blank slate helps to perfectly align the BAF and LRR plotting areas
                 out.list <- col.plot.lrr(samples[cc], bigMat2, dir=dir, centre.chr=chrN[[samples[cc]]],
                                          scheme="mono", col1="white", col2=col2, snp.info=snp.info, plotAdj=rngOn, 
                                          points.only=F, show.fail=show.fail, add.leg=F,scl=scl,
                                          c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=F, ucsc=ucsc, trunc=limz,
                                          ylim=blimz, xlim=c.xlim, rng.mb=rng.mb, main=titl, xlab=xl, ylab="B Allele Frequency")
                 pO <- T
      } else { pO <- F }
      # full LRR data for chromosomes (full colour)
      out.list <- col.plot.lrr(samples[cc], bigMat3, dir=dir, centre.chr=chrN[[samples[cc]]],
                               scheme=scheme, col1=col1, col2=col2, snp.info=snp.info, plotAdj=rngOn, 
                               points.only=pO, show.fail=show.fail, add.leg=F,scl=scl,
                               c.pre.n=c.pre.n, c.post.n=c.post.n, m.smooth=F, ucsc=ucsc, trunc=limz,
                               ylim=blimz, xlim=c.xlim, rng.mb=rng.mb, main=titl, xlab=xl, ylab="B Allele Frequency")
      
      x.coords <- out.list$x.coords; chr.select <- out.list$chr.select
      # add chromosome labels with bad ones coloured red
      if (!LRR & length(Chr)>1) { null.result <- add.chr.top.axis(chr.select, badChrN= NULL, x.coords, nC=length(chr.set), sub=F) }
      #legend("topleft",legend=c("Raw BAF data [coloured by chromosome]"),cex=.7,
      #           pch=c(21),col=c("grey"),bty="n") 
      if(tag.cnvs & is(Cnv)[1]!="RangedData") {
        uu <- draw.cnv.bounds(cnv=Cnv,chr.offset=chr.offset,pos=Pos,cnv.lty=cnv.lty,cnv.col=cnv.col,plot.scl=scl)
      }
    }
    if (bafOverlay & is(pfb.rng)[1]=="RangedData" & BAF) {
      # use pfb file information to plot BAF group averages across markers
      chromo <- chrN[[samples[cc]]]
      if(zoom) {
        ss <- pfb.rng[chromo]; ss <- ss[start(ss)>=min(Pos) & end(ss)<=max(Pos),]
        lines(start(ss)/scl,ss$PFB,col=pfb.col,lty=ov.lty)
      } else {
        ofs <- out.list$offset
        lines((pfb.rng$gindx-ofs)/scl,pfb.rng$PFB,col=pfb.col,lty=ov.lty)
      }
      legend("right",legend="Av. BAF",lwd=1,bty="n",col=pfb.col,lty=ov.lty,cex=0.8)
    }
    if (hzOverlay & is(hz.rng)[1]=="RangedData" & BAF) {
      # use hz file information to plot BAF group averages across markers
      chromo <- chrN[[samples[cc]]]
      if(zoom) {
        ss <- hz.rng[chromo]; ss <- ss[start(ss)>=min(Pos) & end(ss)<=max(Pos),]
      #  print(ss)
        lines(start(ss)/scl,ss$het,col=hz.col,lty=ov.lty)
      } else {
        ofs <- out.list$offset
        lines((hz.rng$gindx-ofs)/scl,hz.rng$het,col=hz.col,lty=ov.lty)
      }
      legend("left",legend="Av. Heterozygosity",lwd=1,bty="n",col=hz.col,lty=ov.lty,cex=0.8)
    }
    cat(".")  # progress dots
  }
  if(plot.to.file) { dev.off() ; cat("~wrote file:",ofn) }
  return(NULL)
}



#' Smooth a dataseries by the median or mean.
#' 
#' Smooth a dataseries towards the median or the mean
#' Vary the window-size and the ratio of median vs original signal
#'
#' @param X numeric vector of data to be smoothed
#' @param n the window size for calculating the median/mean
#' @param method default is median, otherwise 'mean'
#' @param ratio the ratio of how much median/mean to mix in with 1 part
#'  of raw signal. E.g, ratio=1 would give half and half.
#' @return returns vector X rescaled towards the median/mean
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@cimr.cam.ac.uk}
#' @examples
#' tab <- rbind(1:16,MedianScale(1:16,n=4,ratio=1),
#'     MedianScale(1:16,n=4,ratio=10^8)) # large 'ratio' value gives ~median
#' rownames(tab) <- c("raw.X","scaled.result","window.median"); tab
#' sd10 <- rnorm(10000,0,10)
#' summary(sd10) # raw data
#' summary(MedianScale(sd10,method="mean")) #scaled data
MedianScale <- function(X, n=100, method="median", ratio=10)
{
  ## smooth data series by median / mean
  rem <- length(X) %% n
  YorN <- (rem != 0) 
  n.chunks <- (length(X) %/% n) + as.integer(YorN)
  if (YorN) {  X <- c(X,rep(NA,times=(n-rem))) }
  Xmat <- matrix(X,nrow=n)
  if(method=="median") {
    meds.out <- apply(Xmat,2,median,na.rm=T) 
  } else {
    meds.out <- apply(Xmat,2,mean,na.rm=T)
  }
  meds.out <- rep(meds.out,each=n)[1:length(X)]
  meds.out <- ((ratio*meds.out)+X)/(ratio+1)
  return(meds.out)
}



get.dat.for.ranges <- function(ranges,dir,pcCor=T,BAF=F,snp.names=T,snp.pos=F,genome=F) {
  ## based on ranges object, provide a plumbCNV directory and choose whether
  # LRR or BAF, and if LRR, whether pre/post pc-correction
  # also choose chr/genome coords, whether use snp.names as labels, and whether
  # to return in x,y format with snp positions.
  dir <- validate.dir.for(dir,c("big","ano"))
  DT <- read.data.tracker(dir)
  typ <- if(BAF) { "big.baf" } else { if(pcCor) { "big.pcc" } else { "big.qc" } }
  fn <- getSlot(DT, typ)
  if(!file.exists(fn)) { warning("bigMatrix datafile ",typ,"not found"); return(NULL) }
  bigMat <- getBigMat(fn,dir$big)
  snp.info <- read.snp.info(dir)
  snp.info <- snp.info[snp.info$QCfail==0,]
  first.last.snps <- range.snp(snp.info,ranges)
  inf <- x.y.for.snp.range(first.last.snps,bigMat,snp.info,genome=genome,unord=BAF)
  besr <- big.extract.snp.ranges(first.last.snps,ranges$id,bigMat,snp.info=snp.info)
  if(snp.pos) { besr.x <- besr }
  if(snp.names) { 
    for(cc in 1:length(besr)) {
      names(besr[[cc]]) <- inf$SNP[[cc]]
      if(snp.pos) { besr.x[[cc]] <- inf$Pos[[cc]] }
    }
  }
  if(snp.pos){
    return(list(x=besr.x,y=besr))
  } else {
    return(besr)
  }
}


x.y.for.snp.range <- function(fls,bigMat,snp.info,genome=F,unord=F) {
  ## for a frame of first and last snps (fls), extract these
  # ranges from bigMat, either in genome coords, or just in chr-relative coords
  if(is.null(dim(fls))) { dim(fls) <- c(1,length(fls)) }
  if(unord) { snp.info <- toGenomeOrder(snp.info) }
  rns <- rownames(snp.info)
  rnb <- rownames(bigMat)
  locs <- if(genome) { genoPos(snp.info) } else { start(snp.info) }
  stz <- match(fls[,1],rnb)
  enz <- match(fls[,2],rnb)
  if(unord) {
    stz <- match(fls[,1],rns)
    enz <- match(fls[,2],rns)
    lst <- vector("list",length(stz))
    for (cc in 1:length(stz)) {
      lst[[cc]] <- rns[stz[cc]:enz[cc]]
    }
  }
  ind <- pos <- nms <- vector("list",nrow(fls))
  for (cc in 1:nrow(fls)) {
    if(unord) { 
      ind[[cc]] <- narm(match(lst[[cc]],rnb))
    } else {
      ind[[cc]] <- stz[cc]:enz[cc]
    }
    nms[[cc]] <- rnb[ind[[cc]]]
    sind <- match(nms[[cc]],rns)
    nms[[cc]] <- nms[[cc]][!is.na(sind)]
    pos[[cc]] <- locs[narm(sind)]
  }
  return(list(SNP=nms,Pos=pos,Index=ind))
}


col.plot.lrr <- function (ID, bigMat, snp.info=NULL, centre.chr=1:22, rng.mb=NA, show.fail=F, add.leg=F,
                          plotAdj=F, samp.chr.means=NULL, whole.genome=F, points.only=F, trunc=NA,
                          c.pre.n=2, c.post.n=2, m.smooth=F, scl=10^6, col1="black", col2="grey",
                          ratio=10, ucsc="hg18", dir="",scheme=c("mono","alt","unique","hilite"),...) 
{
  # Color Type Plot
  # colour type plot of LRR data for a sample over a given range, optional smoothing, overlay
  # of means, etc.
 # print(head(snp.info))
  #targ.chr.ref sets the x axis scale start point (e.g., relative to start of what chromosome)
  if(is.null(dim(bigMat))) { stop("bigMat is dimensionless") }
  if(is(snp.info)[1]=="RangedData" & (all(c("QCfail","gindx") %in% colnames(snp.info)))) { 
    # reduce info to only snps in the data
    if(any(!rownames(bigMat)%in% rownames(snp.info))) { warning("some SNP names in bigMat were not in snp.info"); return(NULL) }
    snp.info <- snp.info[(rownames(snp.info) %in% rownames(bigMat)),]
    if(!isGenomeOrder(snp.info)) { warning("bigMat rows (snps) are not in genome order"); return(NULL) }
    uv <- tolower(universe(snp.info)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  #  if(force.autosomes) { snp.info <- select.autosomes(snp.info) }
    chr.set <- chr.nums.of.info(snp.info); nC <- length(chr.set)
  } else { warning("snp.info invalid, plot skipped"); return(NULL) }
  chrLens <- get.chr.lens(dir,ucsc=ucsc)[chr.set]
  chrStarts <- c(0,cumsum((chrLens)))[1:nC] # genome position of chromosome starts
  # extract samples' LRR data from big matrix
  loc.col <- match(ID,colnames(bigMat))
  LRR.mat <- bigMat[match(rownames(snp.info),rownames(bigMat)),loc.col]
  targ.chr.ref <- which(chr.set %in% paste(centre.chr[1]))  
  if(!whole.genome & all(!is.na(rng.mb))) {
    XXv <- as.logical(snp.info[targ.chr.ref]$QCfail)
    XX <- start(snp.info[targ.chr.ref])/scl
    chr.x.labs <- NULL
    #print(rng.mb)
    LRR.mat <- LRR.mat[paste(chr(snp.info))==paste(targ.chr.ref)]
    ## PLAYING AROUND HERE - ALL A MESS - AT LEAST CHR22
  } else {
    if(whole.genome) { whole.genome <- F; warning("whole.genome set to FALSE as rng.mb was entered") }
    if(!plotAdj) { targ.chr.ref <- NULL }
    XXv <- as.logical(snp.info$QCfail)
    XX <- (snp.info$gindx)/scl
    chr.x.labs <- ((chrStarts[1:nC]+(chrLens[1:nC]/2))/scl)
  }
  #print(length(XX))
  #XX <- (CHR.INFO$long.gnmIndx/scl)
  # extract samples' LRR data from big matrix
#  loc.col <- match(ID,colnames(bigMat))
#  LRR.mat <- bigMat[1:nrow(bigMat),loc.col] #CHR.INFO$long.gnmIndx
  #LRR.dat <- sub.big.matrix(bigMat, firstCol = loc.col, lastCol = loc.col, backingpath=b.dir)
  #LRR.mat <- LRR.dat[1:nrow(LRR.dat)]
  # set indices if only plotting a limited range of genome, eg only 2 adjacent chrs, etc
  if (plotAdj) {
    ii <- range(as.integer(centre.chr))
    chr.select <- max(1,(ii[1]-c.pre.n)):min((ii[2]+c.post.n+1),nC); #print(ii); print(chr.select)
    rng <- chrStarts[range(chr.select)]/scl
    rng[1] <- floor(rng[1]); rng[2] <- ceiling(rng[2])
    #print(rng)
    select <- c( head(which(XX>rng[1]),1): tail(which(XX<rng[2]),1) ) 
  } else {
    select <- c(1:length(LRR.mat)) 
    if(whole.genome) { chr.select <- 1:nC } else { chr.select <- targ.chr.ref }
  }
  if(!"color" %in% colnames(snp.info)) {
    snp.info <- add.color.to.snp.info(snp.info=snp.info,scheme=scheme,col1=col1,col2=col2,hilite=targ.chr.ref)
  }
  # apply smoothing (mean/median) to LRR data
  if(length(LRR.mat)!=length(XX)) { stop("mismatch in length between x axis and bigMat file") }
  if (m.smooth) { 
    yy <- MedianScale(LRR.mat[select],ratio) 
    xx <- MedianScale(XX[select],ratio=0)
    ccc <- (snp.info$color[select])[seq(1,length(LRR.mat[select]),length.out=length(xx))]
  } else {
    xx <- XX[select]; yy <- LRR.mat[select] ; 
   # print("XX"); print(xx)
    ccc <- snp.info$color[select] 
  }
  xxv <- XXv[select]
  # Do the large colour plot
  offset <- (chrStarts[targ.chr.ref]/scl)
  if (plotAdj) 
  {    
    # ie adjust scale if zooming into a subset of chromosome
    #print("yes")
    xx <- xx-offset
    chr.x.labs <- chr.x.labs-offset
  } #else { offset <- 0 }
  if(all(!is.na(rng.mb))) {
    #nb: targ.chr.ref is the chr number in snp.info, centre.chr[1] is the actual number (usually but not always the same)
    rng.mb <- suppressWarnings(force.chr.pos(rng.mb*scl,centre.chr[1],snp.info))/scl
    if(all(!is.na(rng.mb))) {
      # filter just to range within the chromosome
      #print("rng.mb"); print(rng.mb)
      select.sub <- (xx>= rng.mb[1] & xx<=rng.mb[2])
      yy <- yy[select.sub]
      ccc <- ccc[select.sub]
      xxv <- xxv[select.sub]
      xx <- xx[select.sub]
    } else {
      warning("invalid range entered"); return(NULL)
    }
  }
  yyp <- yy[!xxv]; cccp <- ccc[!xxv]; xxp <- xx[!xxv] # passing snps
  yyf <- yy[xxv]; cccf <- ccc[xxv]; xxf <- xx[xxv] # failing snps
  
  if(!show.fail | (length(yyf)<1)) { yy <- yyp; xx <- xxp; ccc <- cccp } else { yy <- yyf; xx <- xxf; ccc <- "red" }
  mypch <- 19; mycex <- 1; if(length(yy) > 10) { mycex <- .75 }; 
  if(length(yy) > 100) { mycex <- .5 }; if(length(yy) > 1000) { mycex <- .2 }; if(length(yy) > 2000) { mycex <- 1; mypch="." }
#  print(length(xx)); print(summary(xx)); print(table(xxv))
  if(length(xx)>0) {
    if(!all(is.na(trunc))) {
      # truncate y values that go off the plot so they'll fit
      if(length(trunc)==2 & is.numeric(trunc)) {
        yy[yy>max(trunc)] <- max(trunc); yy[yy<min(trunc)] <- min(trunc)
      }
    }
    if(points.only) {
      points(xx,yy,col=ccc,pch=mypch,cex=mycex,...)
    } else {
      plot(xx,yy,col=ccc,pch=mypch,type="p",cex=mycex,...) 
    }
    #print("xx"); print(xx); print("yy"); print(yy)
  } else { warning("no points left to graph") }
  if(add.leg & show.fail) { legend("topright",legend="QC failing SNPs",pch=19,col="red",cex=0.8,bty="n") }
  if(show.fail & (length(yyf)>0)) { 
    yy <- yyp; xx <- xxp; ccc <- cccp ; mycex <- 1.5; mypch <- "."
#    print(length(xx)); print(summary(xx));
    mypch <- 19; mycex <- 1; if(length(yy) > 10) { mycex <- .75 }; if(length(yy) > 100) { mycex <- .5 }
    if(length(yy)>0) {  points(xx,yy,col=ccc,pch=mypch,cex=mycex,...)  }
  }
  # add lines between the means of each chromosome
  if (!is.null(samp.chr.means) & length(chr.select)>1 & whole.genome) {
    lines(chr.x.labs[chr.select],samp.chr.means[chr.select]) }
  out <- list(chr.x.labs,chr.select,select,offset)
  names(out) <- c("x.coords","chr.select","select","offset")
  return(out)
}


dlrs <- function(X,na.rm=T)
{
  # calculate dlrs
  out <- (sd(diff(X),na.rm=na.rm))
  return(out/sqrt(2))
}


get.exon.annot <- function(dir=NULL,ucsc="hg18",range.out=T,list.out=F) {
  ## load exon annotation (store locally if not there already)
  from.scr <- T
  if(!is.null(dir)) {
    ex.fn <- cat.path(dir$ano,"exonAnnot.RData")
    if(file.exists(ex.fn)) {
      ts <- get(paste(load(ex.fn)))
      if((range.out|!list.out) & is.data.frame(ts)) { from.scr <- F }
      if((list.out & !range.out) & is.list(ts)) { from.scr <- F }
    }
  } 
  must.use.package("GenomicFeatures",T)
  if(from.scr) {
    must.use.package("gage",T)
    # get transcripts from UCSC table 'knownGene'
    success <- tryCatch(txdb <- makeTranscriptDbFromUCSC(genome=ucsc,
                   tablename="knownGene")  ,error=function(e) { F } )
    if(is.logical(success)) { 
      if(!success) {
        warning("Couldn't reach ucsc website! try again later or, \n",
                "if in europe/uk, there may still be a bug in rtracklayer; \n",
                "Installing the latest version of R and bioconductor\n",
                "and running biocLite('rtracklayer'), should fix this")
        return(NULL) }
    }
    ts = transcriptsBy(txdb, by="gene")
    data(egSymb)
    select <- match(names(ts),egSymb[,1])
    names(ts)[!is.na(select)] <- egSymb[,2][select[!is.na(select)]]
  }
  if(!list.out | range.out) {
    if(from.scr) {
      ts <- as.data.frame(ts)
      chrs <- paste(ts$seqnames); chrs <- gsub("chr","",chrs)
      ts$seqnames <- chrs
      if(all(colnames(ts)==c("element","seqnames","start","end","width","strand","tx_id","tx_name"))) {
        colnames(ts) <- c("gene","chr","start","end","width","strand","txid","txname")
      } else {
        cat(" unexpected colnames found using makeTranscriptDbFromUCSC()\n")
        if(range.out) { cat(" therefore returning data.frame instead of RangedData object\n") ;
                        range.out <- F }
      }
      if(exists("ex.fn")) { save(ts,file=ex.fn) }
    }
    if(range.out) {
      ts <- RangedData(ranges=IRanges(start=ts$start,end=ts$end),
                       space=ts$chr,gene=ts$gene, strand=ts$strand,
                       txid=ts$txid, txname=ts$txname,universe=ucsc)
      ts <- toGenomeOrder(ts,strict=T)
    }
  } else {
    if(from.scr) { if(exists("ex.fn")) { save(ts,file=ex.fn) } }
  }
  return(ts)
}


get.gene.annot <- function(dir=NULL,ucsc="hg18",range.out=T,duplicate.report=F) {
  # faster than exon, but only contains whole gene ranges, not transcripts
  # allows report on duplicates as some might be confused as to why some genes
  # have more than one row in the listing (split across ranges usually)
  must.use.package(c("biomaRt","genoset","gage"),T)
  from.scr <- T
  if(!is.null(dir)) {
    gn.fn <- cat.path(dir$ano,"geneAnnot.RData")
    if(file.exists(gn.fn)) {
      dat <- get(paste(load(gn.fn)))
      from.scr <- F
    }
  }
  if(from.scr) {
    if(ucsc=="hg18") {
      ens <- useMart("ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     host="may2009.archive.ensembl.org",
                     path="/biomart/martservice",
                     archive=FALSE)
    } else {
      ens <- useMart("ensembl")
    }
    ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
    data(egSymb)
    dat <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
                                "start_position", "end_position", "band"), filters = "hgnc_symbol",
                 values = egSymb[,2], mart = ens)
    if(exists("gn.fn")) { save(dat,file=gn.fn) }
  } 
  if(range.out) {
    outData <- RangedData(ranges=IRanges(start=dat$start_position,end=dat$end_position),
                          space=dat$chromosome_name,gene=dat$hgnc_symbol, band=dat$band, universe=ucsc)
    outData <- toGenomeOrder(outData,strict=T)
    if(duplicate.report) {
      dup.genes <- gene.duplicate.report(outData)
    }
  } else {
    outData <- dat; colnames(outData) <- c("gene","chr","start","end","band")
  }
  return(outData)
}


gene.duplicate.report <- function(ga,full.listing=F) {
  # for a RangedData object, report on any multiple listings for the same gene
  if(is(ga)[1]!="RangedData") { warning("not a RangedData object") ; return(NULL) }
  if("gene" %in% c)
    gene.col <- (which(tolower(colnames(ga)) %in% c("gene","genes","geneid")))
  if(length(gene.col)>0) { gene.col <- gene.col[1] } else { warning("no 'gene' column"); return(NULL) }
  colnames(ga)[gene.col] <- "gene" #force this colname
  if(duplicate.report) {
    culprits <- unique(ga$gene[which(duplicated(ga$gene))])
    n.gene.multi.row <- length(culprits)
    culprit.ranges <- ga[ga$gene %in% culprits,]
    total.culprit.rows <- nrow(culprit.ranges)
    start.same.ct <- end.same.ct <- 0; which.ss <- NULL
    for (cc in 1:length(culprits)) { 
      mini <- (ga[ga$gene %in% culprits[cc],]) 
      if(full.listing) {
        cat("Gene:",culprits[cc],"# same start:",anyDuplicated(start(mini)),
            "# same end:",anyDuplicated(end(mini)),"\n") }
      start.same.ct <- start.same.ct+anyDuplicated(start(mini))
      end.same.ct <- end.same.ct+anyDuplicated(end(mini))
      if(anyDuplicated(start(mini)) | anyDuplicated(end(mini))) { which.ss <- c(which.ss,cc) }
    }
    cat(" genes with split ranges:\n"); print(culprits,quote=F); cat("\n")
    cat(" genes with same start or end:\n"); print(culprits[which.ss],quote=F); cat("\n")
    cat(" total gene-segments with same start",start.same.ct,"; total with same end:",end.same.ct,"\n")
  }
  return(culprits)
}


dgv.subset.cov.by.n.snps <- function(n.snps=10,dir=NULL,snp.set=NULL,dgv=NULL,max.freq=NA,
                                     snp.range=1:20,do.plot=T,add.col=F,rmv.fail=T,...) {
  # return the subset of the DGV covered on the current chip by 'n.snps', or if 'add.col' then just adds a snp count column to dgv dataset
  if(is(dgv)[1]!="RangedData") { if(!is.null(dir)) { dgv <- get.dgv.ranges(dir=dir) ; cat("read dgv from file\n")
                        } else { warnings("no valid dgv or dir parameter"); return(NULL) } }
  if(!is.na(max.freq) & ("frequency" %in% colnames(dgv))) {
    # if specified then filter dgv for CNVs below a maximum frequency in the sample
    dgv <- dgv[!is.na(dgv$frequency),][narm(dgv$frequency<max.freq),]
  }
  if(is(snp.set)[1]!="RangedData") { if(!is.null(dir)) { snp.set <- read.snp.info(dir=dir) ; cat("read snp.info from file\n")
                        } else { warnings("no valid dgv or dir parameter"); return(NULL) } }
  if(!is.numeric(n.snps)) { warning("n.snps must be an integer scalar value") }
  if("QCfail" %in% colnames(snp.set)) {
    if(rmv.fail){
      snp.set$QCfail[is.na(snp.set$QCfail)] <- 1
      cat(" removing",nrow(snp.set)-length(which(snp.set$QCfail==0)),"QC failing SNPs from snp.set\n")
      snp.set <- snp.set[snp.set$QCfail==0,]
    }
  }
  dg <- find.overlaps(dgv[,-c(1:ncol(dgv))],ref=snp.set[,-c(1:ncol(snp.set))],...) # e.g, n.cores & -c() speeds up
  dgL <- sapply(dg,clean.fn,fail=0,fn=length)
  mydg <- function(n) { sapply(n,function(n) { length(which(dgL>n)) } ) }
  if(do.plot & length(snp.range)>1 & is.numeric(snp.range)) {  
    plot(x=snp.range,mydg(snp.range)/nrow(dgv),type="l",xlab="number of snps",ylab="% of DGV") }
  if(add.col) {
    dgv[["chip.coverage"]] <- 0
    dgv$chip.coverage[as.numeric(names(dg))] <- dgL
    return(dgv)
  } else {
    return(dgv[as.numeric(names(dg)[dgL>round(n.snps[1])]),])
  }
}


count.cnvs <- function(cnv.hit.list,cnv.sample.map=NULL,by.phenotype=F,tot.cnvs=NA,
                       content.txt="reference region",print.result=F) {
  ## uses result from 'find.overlaps()'
  ## a function to count the number of CNVs, genes, exon, etc. sapply?
  ## number of unique gene/ex/dgv per CNV; number of gene/ex/dgv hits across all CNVs; 
  ## number of unique gene/ex/dgv across all CNVs
  ## all these could have percentage filter applied relative to reference or source
  if(!is.list(cnv.hit.list)) { stop("Error: cnv.hit.list should be a list, 1 element per CNV, each listing ids/regions hit") }
  n.cnvs <- length(cnv.hit.list); ph.mode <- F
  ## Test whether sample object has a valid column for id, and derive a unique CNV ID
  if(!is.null(dim(cnv.sample.map))) {
    samp.mode <- T
    if(!is.null(rownames(cnv.sample.map))) { cids <- rownames(cnv.sample.map) } else { cids <- paste(1:nrow(cnv.sample.map)) }
    id.col <- (which(tolower(colnames(cnv.sample.map)) %in% c("id","ids","sample","samples","subjects")))
    if(length(id.col)>0) {
      sids <- cnv.sample.map[,id.col[1]]; if(is(sids)[1]=="RangedData") { sids <- sids[[1]] }
    } else { 
      warning("no id column found in cnv.sample.map, no sample-based summaries will be generated")
      samp.mode <- F
    } 
    if(by.phenotype) {
      # look for phenotype data in 'cnv.sample.map'
      ph.col <- (which(tolower(colnames(cnv.sample.map)) %in% c("pheno","phenotype")))
      if(length(ph.col)>0) {
        ph.dat <- cnv.sample.map[,ph.col]; ph.mode <- T
        if(is(ph.dat)[1]=="RangedData") { ph.dat <- ph.dat[[1]] } # do this to remove 'ranges' and just get a vector
      } else { 
        warning("no phenotype column found in cnv.sample.map, no phenotype-based summaries will be generated")
        ph.dat <- NULL; ph.mode <- F
      }
    }
  } else {
    if(is.character(cnv.sample.map) | is.numeric(cnv.sample.map)) { 
      sids <- cnv.sample.map ; cids <- paste(1:nrow(cnv.sample.map)); samp.mode <- T
    } else {
      samp.mode <- F
    }
  }
  if(samp.mode) { 
    if(length(sids)!=length(cids)) { stop(" sample id length must match cnv id length\n") }
    cnv.num.per.cnv.hit <- match(names(cnv.hit.list),cids)
    if(any(is.na(cnv.num.per.cnv.hit))) { 
      if(length(cids)!=n.cnvs) {
        stop("Error: unequal length between cnv.sample.map and cnv.hit.list, and cnv.ids don't match")
      } else {
        nmm <- length(which(is.na(cnv.num.per.cnv.hit)))
        warning(paste("",nmm,"/",n.cnvs," mismatchs in cnv.sample.map and cnv.hit.list names, assuming exact order match\n",sep=""))
        cnv.num.per.cnv.hit <- 1:n.cnvs
      }
    }
    sids <- sids[cnv.num.per.cnv.hit]; cids <- cids[cnv.num.per.cnv.hit]
    if(!is.null(ph.dat)) { ph.dat <- ph.dat[cnv.num.per.cnv.hit] }
  }
  ## do overall calculations
  overall.result <- basic.cnv.counts(cnv.hit.list, tot.cnvs=tot.cnvs, content.txt=content.txt)
  ## do phenotype calculations (essentially call function again filtering for subset)
  pheno.result <- NULL
  if(ph.mode) {
    phenos <- sort(unique(ph.dat))
    n.phenos <- length(phenos); pheno.result <- vector("list",n.phenos)
    if(length(phenos)>1) {
      for(dd in 1:length(phenos)) {
        cnv.hit.sublist <- cnv.hit.list[which(ph.dat==phenos[dd])]
        # print(length(cnv.hit.list)); print(length(cnv.hit.sublist)); print(length(ph.dat)); print(dim(ph.dat))
        #pheno.result[[dd]] <- count.cnvs(cnv.hit.sublist,cnv.sample.map=NULL,by.phenotype=F,tot.cnvs=NA,content.txt=content.txt)
        pheno.result[[dd]] <- basic.cnv.counts(cnv.hit.sublist, tot.cnvs=NA, content.txt=content.txt)
        if(samp.mode) {
          ## if also have sample information, then add this to make sample*phenotype summary
          s.result <- sample.cnv.counts(cnv.hit.sublist, sids[which(ph.dat==phenos[dd])], content.txt)[c(3,4)]
          for (j in 1:length(names(s.result))) {
            pheno.result[[dd]][[names(s.result)[j]]] <- s.result[[j]]
          }
        }
      }
      names(pheno.result) <- paste(phenos)
    } else { warning("only 1 phenotype found in CNV object\n no phenotype summary generated") }    
  }
  ## do sample calculations
  if(samp.mode) {
    sample.result <- sample.cnv.counts(cnv.hit.list, sids, content.txt)
  } else {
    sample.result <- NULL
  }
  result.out <- list(overall.result); nm.out <- c("overall.counts")
  if(!is.null(pheno.result)) { result.out <- list(overall.result,pheno.result); nm.out <- c(nm.out,"phenotype.counts") }
  if(!is.null(sample.result)) { result.out <- c(result.out,list(sample.result)); nm.out <- c(nm.out,"sample.counts") }
  names(result.out) <- nm.out
  if(print.result) {
    tl <- paste("Overall",content.txt,"vs CNV count summary") 
    cat("\n",tl,"\n",paste(rep("-",nchar(tl)),collapse=""),"\n",sep="")
    print(overall.result[4:7]) # i.e, skip the long lists in output
    if(!is.null(pheno.result)) {
      tl <- paste(toheader(content.txt),"count by phenotype summary:")
      cat("\n",tl,"\n",paste(rep("-",nchar(tl)),collapse=""),"\n",sep="")
      for (dd in 1:length(pheno.result)) { 
        cat("\nSummary for phenotype =",phenos[dd],":\n")
        print(pheno.result[[dd]][c(c(4:7),c(-1,0)+length(pheno.result[[dd]]))]) 
      } }
    if(!is.null(sample.result)) {
      tl <- paste(toheader(content.txt),"count by sample summary:")
      cat("\n",tl,"\n",paste(rep("-",nchar(tl)),collapse=""),"\n",sep="")
      print(sample.result[3:4]); cat("\n")
    }
  }
  return(result.out)
}


do.CNV.all.overlaps.summary <- function(cnvResults,dir,comps=c(1:5),dbs=c(1:3),...) {
  # pass the results list from plumbCNV() and include any args of 'find.overlaps()' for filtering
  fail <- F
  if(!is.list(cnvResults)) { fail <- T }
  if(length(cnvResults)!=5) { fail <- T }
  if(!all(sapply(cnvResults,is)[1,]=="RangedData")) { fail <- T }
  if(fail) { 
    warning("cnvResults must be a list produced by plumbCNV() containing:\n",
            " 5 RangedData objects: allCNV, allDel, allDup, rareDel, rareDup")
    return(NULL) 
  }
  # initialise values for each run of the loop
  if(!all(comps %in% 1:5)) { comps <- 1:5 }; n.c <- length(comps)
  if(!all(dbs %in% 1:3)) { dbs <- 1:3 }; n.d <- length(dbs)
  DUPs <- c(T,F,T,F,T)[comps]
  DELs <- c(T,T,F,T,F)[comps]
  compz <- c("AllCNV","AllDel","AllDup","RareDEL","RareDUP")[comps]
  set.n <- c(1:5)[comps]
  db <- c("gene","exon","dgv")[dbs]
  lab <- c("Genes","Exons","DGV regions")[dbs]
  dbz <- vector("list",n.d); names(dbz) <- db
  # loop for finding overlaps and doing counts
  for (dd in 1:n.d) {
    # setup list substructures and names
    dbz[[dd]] <- vector("list",2); names(dbz[[dd]]) <- c("overlaps","counts")
    dbz[[dd]]$overlaps <- vector("list",n.c); names(dbz[[dd]]$overlaps) <- compz
    dbz[[dd]]$counts <- dbz[[dd]]$overlaps
    for (cc in 1:n.c) {
      tl <- paste("Running overlap analysis for:",compz[cc],"with",lab[dd])
      cat(tl,"\n")  #,paste(rep("=",nchar(tl)),collapse=""),"\n\n",sep="")
      dbz[[dd]]$overlaps[[cc]] <- oo <- find.overlaps(cnvResults[[set.n[cc]]],DEL=DELs[cc],DUP=DUPs[cc],
                                                      db=db[dd],dir=dir,quiet=T,...)
     # cat(".. done,")
      #print(headl(oo)); print(dim(oo)); print(length(oo)); print(is(oo))
      dbz[[dd]]$counts[[cc]] <- count.cnvs(oo,content.txt=lab[dd],cnv.sample.map=cnvResults[[set.n[cc]]],by.phenotype=T)
     # cat(" counts ... done",cc,"\n")
    }
  }
  return(dbz)
}


## not visible
summary.counts.table <- function(CNV.overlaps.summary, print.result=T, verbose=T) {
  # do a table for overlap with each database at a time
  # columns are type of CNV, rows are type of count
  ############### start local function ###########
  summary.counts.table.1 <- function(dbz, print.result=T, verbose=T) {
    # do a table for overlap for a database
    # this is a function only used locally
    #if(verbose) {  cat("flattening dB...") }
    flat.list <- reduce.list.to.scalars.2(dbz)
    #if(verbose) {  cat("done\n") }
    # derive table size and column names
    n.rows <- length(unlist(flat.list$counts[[1]]))
    n.cols <- length(flat.list$counts)
    r.data <- matrix(integer(0),ncol=n.cols,nrow=n.rows)
    colnames(r.data) <- names(flat.list$counts)
    Category <- Subgroup <- Count <- character(0)
    #if(verbose) { cat(" extracting summary data:\n") }
    #st.t <- proc.time()
    ## loop through extracting data for each column
    for (jj in 1:n.cols) {
      rdat <- integer(0)
      # Add overall counts (must be present)
      overall <- flat.list$counts[[jj]]$overall.counts
      if(length(unlist(overall))>1) {
       # (headl(overall))
        #print(is(overall)); print(length(overall))
        rdat <- c(rdat,as.integer(unlist(overall)))
        if(jj==1) {
          nz <- (names(overall))
          Count <- c(Count,nz)
          Category <- c(Category,rep("OVERALL",length(nz)))
          Subgroup <- c(Subgroup,rep("",length(nz)))
        }
      }
      # Add sample counts if present
      samp <- flat.list$counts[[jj]]$sample.counts
      if(length(unlist(samp))>1) {
        rdat <- c(rdat,as.integer(samp))
        if(jj==1){
          nz <- (names(samp))
          Count <- c(Count,nz)
          Category <- c(Category,rep("SAMPLE",length(nz)))
          Subgroup <- c(Subgroup,rep("",length(nz)))
        }
      }
      # Add phenotype counts if present
      p.cnt <- flat.list$counts[[jj]]$phenotype.counts
      phenos <- names(p.cnt); if(is.null(phenos)) { phenos <- 1:length(p.cnt) }
      if(length(unlist(p.cnt))>1) {
        for (cc in 1:length(p.cnt)) {
          rdat <- c(rdat,as.integer(p.cnt[[cc]]))
          if(jj==1){
            nz <- (names(p.cnt[[cc]]))
            Count <- c(Count,nz)
            Category <- c(Category,rep("PHENOTYPE",length(nz)))
            Subgroup <- c(Subgroup,rep(paste("group",phenos[cc],sep="="),length(nz)))
          }
        }
      }
      r.data[,jj] <- rdat
      #if(verbose) { loop.tracker(jj,n.cols,st.t) }
    }
    Count <- gsub("number of ","",Count)
    result <- cbind(Category,Subgroup,Count,r.data)
    if(print.result) { print(result,quote=F) ; cat("\n") }
    return(result)
  }
  #debug(summary.counts.table.1)
  
  ################ end local function ############
  dbz <- CNV.overlaps.summary
  if(!is.list(dbz)) { stop("Error: CNV.overlaps.summary must be a list") }
  if(!all(names(dbz) %in% c("gene","exon","dgv"))) { warning("CNV.overlaps.summary was designed for sublists called gene, exon, dgv; proceed with caution") }
  n.db <- length(dbz)
  sct <- vector("list",n.db)
  for (cc in 1:n.db) {
    if(print.result) { cat("\nCNV Count Summary for",toheader(names(dbz)[cc]),"\n") }
    sct[[cc]] <- summary.counts.table.1(dbz[[cc]],print.result=print.result, verbose=verbose)
  }
  return(sct)
}


pheno.ratios.table <- function(dir,sum.table)
{
  sample.info <- read.sample.info(dir)
  if(!is.list(sum.table)) { stop("sum.table was not a list") }
  if(all(sapply(sum.table,is)[1,]!="matrix")) { stop("sum.table was not a list of matrices") }
  # NB: for consistency will for now not put CNV-QC failures in the denominator
  p.cnts <- table(sample.info$phenotype[sample.info$QCfail==0 | sample.info$QCfail==6])
  n.tabs <- length(sum.table);
  out <- vector("list",n.tabs); names(out) <- names(sum.table)
  for (dd in 1:n.tabs)
  {
    t1 <- (sum.table[[dd]])
    p1 <- t1[t1[,1]=="PHENOTYPE",-1]
    pnum <- as.numeric(sapply(strsplit(p1[,1],"="),"[",2))
    grps <- unique(pnum)
    vals <- pcs <- vector("list",length(grps)); j <- 1; names(vals) <- names(pcs) <- paste(grps)
    for (cc in grps) {
      vals[[j]] <- matrix(as.numeric(p1[pnum==cc,3:ncol(p1)]),nrow=length(which(pnum==cc)))
      colnames(vals[[j]]) <- colnames(p1)[3:ncol(p1)]; rownames(vals[[j]]) <- p1[1:nrow(vals[[j]]),2]
      pcs[[j]] <- vals[[j]]/p.cnts[paste(cc)]
      j <- j+1
    }
    if(length(grps)>1) {
      ratios <- vector("list",length(grps)-1); j <- 1; names(ratios) <- paste(grps[-1],"/",grps[1],sep="")
      for (cc in grps[-1]) {
        ratios[[j]] <- pcs[[paste(cc)]]/pcs[[paste(grps[1])]]
        j <- j+1
      }
      out[[dd]] <- list(raw.counts=vals,per.sample=pcs,pheno.ratios=ratios)
    } else {
      out <- "no phenotype ratios printed as only 1 group was found"
    }
    
  }
  return(out)
}


find.overlaps <- function(cnv.ranges, thresh=0, geq=T, rel.ref=T, pc=T, ranges.out=F, vals.out=F, vec.out=F, delim=",",
                          ref=NULL, none.val=0, fill.blanks=T, DEL=T, DUP=T, min.sites=0, len.lo=NA, len.hi=NA,
                          db=c("gene","exon","dgv"), txid=F, ucsc="hg18", n.cores=1, dir="", quiet=F) {
  # like 'findOverlaps' but according to a specific percentage
  # specify own reference comparison or use standard 'gene', 'exon' or 'dgv'
  # can filter on several criteria prior to matching (length/snps/dels/dups), or after matching (overlap %)
  # can get a list of overlap hits, or just a filtered RangedData object
  if(is(cnv.ranges)[1]!="RangedData") {
    warning("'cnv.ranges' must be 'RangedData' type; returning null")
    return(NULL)
  } else {
    if(nrow(cnv.ranges)<1) { warning("cnv.ranges contains zero rows, returning NULL"); return(NULL) }
  }
  uv <- tolower(universe(cnv.ranges)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { ucsc <- uv } }
  dir <- validate.dir.for(dir,"ano")
  if(ranges.out) { fill.blanks <- T; none.val <- 0 }
  ## filter the cnv dataset for only DELs or only DUPs if either set false
  if(!DEL & !DUP) { DEL <- DUP <- T ; warning("reset DEL=T, DUP=T as both were F") }
  if("cn" %in% colnames(cnv.ranges)) {
    bef <- nrow(cnv.ranges)
    if(!DUP) { cnv.ranges <- cnv.ranges[cnv.ranges$cn<2,] }
    if(!DEL) { cnv.ranges <- cnv.ranges[cnv.ranges$cn>2,] }
    aft <- nrow(cnv.ranges); if(bef>0 & aft==0) {
      warning("DUP/DEL filter has removed all CNVs from comparison")
      return(NULL)
    }
  } else { if(!DEL | !DUP) { warning("'cn' column not found in cnv.ranges") } }
  ## filter the cnv dataset for a max and min CNV length
  selectlo <- T; selecthi <- T; do.l.filt <- F 
  if(all(!is.na(len.lo))){
    ## make sure lo is lower than hi
    if(all(!is.na(len.hi))) { both <- range(len.lo,len.hi,na.rm=T); len.lo <- both[1]; len.hi <- both[2] }
    selectlo <- width(cnv.ranges)>=len.lo; do.l.filt <- T
  }
  if(all(!is.na(len.hi))){
    selecthi <- width(cnv.ranges)<=len.hi; do.l.filt <- T
  }
  select <- selectlo & selecthi ; 
  if(length(which(select))==0) { warning("lo/hi length selection has removed all CNVs from set"); return(NULL) }
  if(do.l.filt) { cnv.ranges <- cnv.ranges[which(select),] }
  ## filter the cnv dataset for a min number of sites/snps
  if(is.numeric(min.sites) & ("numSnps" %in% colnames(cnv.ranges))) {
    bef <- nrow(cnv.ranges); cnv.ranges <- cnv.ranges[cnv.ranges$numSnps>=min.sites,]; aft <- nrow(cnv.ranges)
    if(aft<bef) { cat(" filtered",bef-aft,"CNVs with fewer than",min.sites,"sites\n") }
  }
  ## check compatibility of settings if vals.out selecteed
  if(vals.out & (thresh>0 | !vec.out)) { 
    warning("vals.out=T only valid when vec.out=T and thresh=0. Setting to FALSE")
    vals.out <- F
  }
  ## overlap with genes or exons or DGV
  ## percentages with respect to the CNVs or with respect to the reference
  #e.g, look in db to find exon, gene or dgv best match ; default is gene
  if(is(ref)[1]=="RangedData" & all(db==c("gene","exon","dgv"))) {
    db <- "custom ranges"
  }
  db <- tolower(paste(db)[1])
  if(is(ref)[1]!="RangedData") {
    ; ano <- "gene" 
    if(length(grep("exon",db))>0) { ano <- "exon" }
    if(length(grep("dgv",db))>0) { ano <- "dgv" }
    if(!quiet) { cat(" loading",ano,"annotation...\n") }
    ref <- switch(ano,exon=get.exon.annot(dir=dir,ucsc=ucsc),gene=get.gene.annot(dir=dir,ucsc=ucsc),
                  dgv=get.dgv.ranges(dir=dir,ucsc=ucsc,compact=T)) # choose which reference to use
    
  } else {
    # ref is already a ranges object for comparison passed in by parameter
    if(any(tolower(db) %in% c("exon","gene","dgv")))  { ano <- db } else { ano <- "custom ranges" ; txid <- F }
  }
  if(ano=="gene" | ano=="exon") { nbg <- T } else { nbg <- F } # use gene names if gene or exon
  if(ano=="dgv") { if(!(DUP & DEL)) {
    if(!DUP) { if("Loss" %in% colnames(ref)) { ref <- ref[!is.na(ref$Loss),] } else { cat("Loss column not found\n")} }
    if(!DEL) { if("Gain" %in% colnames(ref)) { ref <- ref[!is.na(ref$Gain),] } else { cat("Gain column not found\n")} }
  } }
  # rel.query=T gives %'s relative to the query (e.g gene), false gives %'s relative to subject (ie, cnv)
  if(!quiet) { cat(" testing CNV set for overlaps with",ano,"...\n") }
  mm <- overlap.pc(query=ref,subject=cnv.ranges,name.by.gene=nbg,
                   rel.query=rel.ref,fill.blanks=fill.blanks, txid=txid,
                   n.cores=n.cores,none.val=none.val)
  if(vec.out & vals.out & !ranges.out) {
    # summarise results into 1 column separated by commas(or 'delim')
    if(pc) {
      out.vec <- jlapply(mm[[1]],mm[[3]],pc=T,collapse=delim) # print gene=%
    } else {
      out.vec <- jlapply(mm[[1]],mm[[2]],pc=T,collapse=delim) # print gene=width
    }
    return(out.vec)
  }
  if(pc) { if(thresh>1) { thresh <- thresh/100 } ; thresh <- max(0,thresh,na.rm=T) }
  if(geq) { op <- ">=" } else { op <- "<=" }
  if(thresh>0) {
    if(thresh>0) { fun <- op } else { fun <- NULL }
    pass.thresh <- jlapply(list1=mm[[3]],list2=thresh,pc=pc,op=fun) # get list of passing threshold
    out <- jlapply(list1=mm[[1]],list2=pass.thresh,select=T)
  } else {
    out <- mm[[1]]
  }
  if(ranges.out) {
    ## if ranges.out return the set of CNVs implied as a RangedData object
    sel.r <- (sapply(out,length)>0)
    if(!quiet) { cat(" selected",length(which(sel.r)),"overlapping ranges from total set of",length(sel.r),"\n") }
    return(cnv.ranges[which(sel.r),])
  }
  if(vec.out) {
    # concatenate result into a single column of delimited text
    out <- as.character(unlist(sapply(out,function(x) { paste(x,collapse=delim) })))
  }
  return(out)
}


#' Combine the output of two identically shaped lists (Joint-List apply)
#'
#' A wrapper for mapply in the case of two equally dimensioned lists.
#' Also has more advanced functionality for combining output, merging,
#' and displaying values as percentages.
#' For instance, might be lists of overlapping genes with percentages.
#' Allows inputting of an optional function of both lists used to produce a result.
#
#' @examples
#' nms <- c("TED","MARK","KIM","TAI","LEE","TINA","PAM","MAY","FRED","UMA")
#' wds <- list(0,c(1:5),c(10:4)) # list of width vectors
#' lns <- list(2,5,10) # list of length vectors
#' pcs <- mapply("/",wds,lns) # list of percentage vectors
#' ids <- list(nms[1],nms[2:6],nms[10:4]) # list of id/name vectors
#' sel <- jlapply(pcs,0.5,"<"); sel  # tests whether each sublist item is <0.5
#' jlapply(ids,pcs,pc=F)  # returns: c(id1=dc1, id2=dc2, ... )
#' # convert to a percentage, and collapse
#' jlapply(ids,pcs,pc=T,collapse=",")  # prints: "id1=xx%, id2=xx%, ... "
#' # basic usage is same as mapply()
#' all.equal(jlapply(wds,pcs,"/"),mapply("/",wds,pcs))
#' # use custom function 'FUN'
#' jlapply(wds,pcs,FUN=function(x1,x2) { x1[x2>.5] <- -x1[x2>.5]; x1 } )
#' jlapply(wds,sel,select=T,collapse=",") # use selection
jlapply <- function(list1, list2, FUN=NULL, select=F, pc=F, collapse=NULL) {
  use.op <- F; if(!is.null(FUN)) { 
    if(exists("FUN",mode="function")) { use.op <- T } else {
      if(exists(FUN,mode="function")) { use.op <- T } }
  }
  if(!(use.op & length(list2)==1) & !(use.op & (length(list2)==length(list1))) )
  {
    if(length(list1)!=length(list2)) { warning("lists of unequal length, returning NULL"); return(NULL) }
    if(any(sapply(list1,length)!=sapply(list2,length))) { warning("sublists of unequal length(s), returning NULL"); return(NULL) }
  }
  if(!use.op) {
    if(select) {
      if(is.numeric(unlist(list2)) | is.logical(unlist(list2))) {
        FUN <- function(l1,l2) { return(l1[l2]) } 
      } else {
        warning("second list must be logical or numeric when select=T"); return(NULL)
      }
    } else {
      if(pc) {
        FUN <- function(l1,l2) { if(any(l1!="")) { paste(l1,"=",round(l2*100),"%",sep="") } else { "" } }
      } else {
        FUN <- function(l1,l2) { if(any(l1!="")) { paste(l1[l1!=""],"=",l2[l1!=""],sep="") } else { "" } }
      }
    }
  }
  new.list <- mapply(FUN,list1,list2)
  if(is.character(collapse)) {
    add.delim <- function(x) { paste(x,collapse=collapse) }
    new.list <- as.character(unlist(sapply(new.list,add.delim)))
  }
  return(new.list)
}


plot.ranges <- function(rangedat,labels=NULL,do.labs=T,skip.plot.new=F,...) {
  # plot a set of ranges from a RangedData object
  #if(is(ranges)[1]!="RangedData") { warning("need RangedData object") ; return(NULL) }
  chk <- chr.nums.of.info(rangedat)
  if(length(chk)>1) { warning(length(chk)," chromosomes in 'rangedat', only using the first, chr",chk[1]) ; rangedat <- rangedat[1] }
  xl <- range(c(start(rangedat),end(rangedat)))
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  nr <- nrow(rangedat); if(is.null(nr)) { nr <- length(rangedat) }
  yl <- c(0,(nr+2))
  colz <- get.distinct.cols(nr); if(nr>22) { colz <- rep("black",nr) }
  if(is.null(labels)) { lab <- rownames(rangedat) } else { lab <- rangedat[[labels]] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) }
  if(!skip.plot.new) {
    plot(x=c(start(rangedat[1,]),end(rangedat[1,])),y=c(1,1),xlim=xl,ylim=yl,type="l",col=colz[1],...)
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      lines(x=c(start(rangedat[cc,]),end(rangedat[cc,])),y=c(cc,cc),col=colz[cc])
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      text(x=start(rangedat[cc,]),y=cc+0.5,labels=lab[cc],cex=0.6,pos=4,offset=0)
    }
  }
}


overlap.pc <- function(query,subject,name.by.gene=F,rel.query=T,fill.blanks=T,
                       text.out=F,delim=",",none.val=0,n.cores=1, txid=F, prog=F) {
  # uses 'findOverlaps' and reports the percentage overlap for each region/etc
  # most sense when: x is reference (e.g, genes); y is regions of interest, eg.  CNVs
  must.use.package(c("genoset","IRanges"),T)
  if(is(query)[1]!="RangedData" | is(subject)[1]!="RangedData") {
    warning("'query' and 'subject' must both be 'RangedData' type; returning null")
    return(NULL)
  }  
  blnk.out <- list(ids=NULL,width=NULL,pc=NULL)
  if(nrow(query)<1 | nrow(subject)<1) { return(blnk.out) }
  # force this so order can be constructed with respect to the initial order
  subj.nrow <- nrow(subject)
  rownames(subject) <- paste(1:subj.nrow)
  # query names only matter for result text, so only add names if no text.
  # genes/exons have a separate column treated specially as otherwise there are duplicate entries not suited to rownames
#print(length(which(is.na(rownames(query)))))
  if(is.null(rownames(query))) { rownames(query) <- paste(1:nrow(query)) }
  ## reduce dimension of subject; 
  to.cut <- which(tolower(colnames(subject)) %in% c("cn","numsnps","fid","score"))
  if(length(to.cut)>0) { subject <- subject[,-to.cut] } # removing unnecessary cols speeds up
  to.cut2 <- which(!tolower(colnames(query)) %in% c("gene","txid","txids","txname","txnames"))
  if(length(to.cut2)>0) { query <- query[,-to.cut2] } # removing unnecessary cols speeds up
  xx=toGenomeOrder(select.autosomes(query),strict=T)
  yy=toGenomeOrder(select.autosomes(subject),strict=T)
  chr.set.x <- chr.nums.of.info(xx)
  chr.set.y <- chr.nums.of.info(yy)
  common <- which(sort(chr.set.x) %in% sort(chr.set.y))
  if(length(common)<1 | nrow(xx)<1 | nrow(yy)<1) { return(blnk.out) }
  chr.set <- paste(sort(chr.nums.of.info(xx[paste(sort(chr.set.x)[common])])))
  x <- ranges(xx[chr.set]); y <- ranges(yy[chr.set]) # keep only common + convert to IRangeslist
  xx <- xx[chr.set]; yy <- yy[chr.set] # keep only common in originals (these store rownames, unfiltered positions)
  n.chr <- length(chr.set)
  by.chr.list <- vector("list",n.chr); names(by.chr.list) <- chr.set
  if(prog) { cat("|") }
  if(n.cores>1) {
    must.use.package("multicore")
    by.chr.list <- multicore::mclapply(X=1:n.chr, FUN=get.overlap.stats.chr, x=x, y=y, xx=xx, yy=yy, 
                                       name.by.gene=name.by.gene, rel.query=rel.query, txid=txid, mc.cores=n.cores,prog=prog)
  } else {
    by.chr.list <- lapply(X=1:n.chr, FUN=get.overlap.stats.chr, x=x, y=y, xx=xx, yy=yy, 
                          name.by.gene=name.by.gene, rel.query=rel.query, txid=txid, prog=prog)
  }
  if(prog) { cat("|\n") }
#  return(by.chr.list) ; # save(by.chr.list,file="by.chr.list.RData")
  namez <- as.character(unlist(sapply(lapply(by.chr.list,"[[",1),names)))
#  namez <- as.character(unlist(sapply(by.chr.list,names)))
  idz <- unlist(sapply(by.chr.list,"[[",1),recursive=F)
  if(is.null(idz)) { if(name.by.gene) { warning("blank names: nb: name.by.gene should be false if gene names not present") } }
  widz <- unlist(sapply(by.chr.list,"[[",2),recursive=F)
  pcz <- unlist(sapply(by.chr.list,"[[",3),recursive=F)
  if(is.null(widz) | is.null(pcz)) { nomatches <- T } else { nomatches <- F }
#  print(idz); print(widz); print(pcz)
#  print(head(namez))
  if(!nomatches) { names(widz) <- names(idz) <- names(pcz) <- namez }
  by.chr.list <- list(ids=idz,width=widz,pc=pcz)
#  print(headl(by.chr.list))
  if(name.by.gene) { names(by.chr.list)[1] <- "genes" }
  if(fill.blanks) {
    replc <- list("",none.val,none.val)
    if(length(by.chr.list)==3) {
      for (dd in 1:3) { by.chr.list[[dd]] <- fill.blank.matches(by.chr.list[[dd]],replc[[dd]],full.len=subj.nrow) }
    }
  }
  #if(text.out) {
  #  add.delim <- function(x) { paste(x,collapse=delim) }
  #  out <- as.character(unlist(sapply(out,add.delim)))
  #}
  return(by.chr.list)
}


annot.cnv <- function(cnvResults, gs=NULL, vec.out=T, delim=";", txid=F){
  ## which genes overlap with each CNV - don't forget to account for the empty ones!
  must.use.package("genoset")
  if(is(gs)[1]!="RangedData" | is(cnvResults)[1]!="RangedData") {
    warning("'gs' and 'cnvResults' must both be 'RangedData' type; returning null")
    return(NULL)
  }
  gene.vec <- (find.overlaps(cnvResults,db="gene",ref=gs,vec.out=T,delim=";"))
  if(nrow(cnvResults)==length(gene.vec)) { cnvResults[["gene"]] <- gene.vec ; return(cnvResults) } else { 
     stop("gene vector returned doesn't match nrow of cnvResult") }
  ## aLLL the rest is useless?
#   rownames(cnvResults) <- paste(1:nrow(cnvResults))
#   olp <- findOverlaps(gs,cnvResults)
#   gene.nums.per.cnv <- tapply(queryHits(olp),subjectHits(olp),c)
#   match.num.to.names <- function(num) { return(gs$gene[num]) }
#   for.exon <- function(num) { return(gs$txname[num]) }
#   if(txid & ("txname" %in% colnames(gs))) { match.num.to.names <- for.exon }
#   out <- sapply(gene.nums.per.cnv,match.num.to.names)
#   # ensure a blank entry for cnvs with no overlaps
#   no.gene.list <- which(!c(1:nrow(cnvResults)) %in% names(out))
#   if(length(no.gene.list)>0) {
#     none.list <- vector("list",length(no.gene.list))
#     names(none.list) <- rownames(cnvResults)[no.gene.list]
#     out <- c(out,none.list)
#   }
#   # unlist into a vector separated by eg, semicolons
#   if(vec.out) {
#     add.delim <- function(x) { paste(x,collapse=delim) }
#     out <- as.character(unlist(sapply(out,add.delim)))
#   } 
#   return(out)
}


Ranges.to.txt <- function(chr.list) {
  ## takes a RangedData object from some annotation lookup functions and converts to standard text positions
  if(is(chr.list)[1]!="RangedData") { stop("Not a RangedData object")}
  text.out <- paste("chr",chr(chr.list),":",format(start(chr.list),scientific=F,trim=T),"-",
                    format(end(chr.list),scientific=F,trim=T),sep="")
  return(text.out)
}


data.frame.to.ranged <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,ucsc="hg18") 
{
  ## convert any data frame with chr,start,end data into a RangedData object
  must.use.package(c("genoset","IRanges"),T)
  if(!is.matrix(dat)) { dat <- as.data.frame(dat) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,start,end,chr,width)
  if(!all(key.nms %in% colnames(dat))) { warning("columns not found") }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { id <- rownames(dat) } else { id <- paste(1:nrow(dat)) } }
  if(!is.null(chr)) { ch <- gsub("chr","",dat[[chr]],ignore.case=T) } else { ch <- NULL }
  if(!is.null(start)) { st <- dat[[start]] } else { st <- NULL }
  if(!is.null(end)) { en <- dat[[end]] } else { en <- NULL }
  if(!is.null(width)) { wd <- dat[[width]] } else { wd <- NULL }
  outData <- RangedData(ranges=IRanges(start=st,end=en,names=id),space=ch,universe=ucsc[1])
  outData <- toGenomeOrder(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to resort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      outData[[more.cols[cc]]] <- dat[[more.cols[cc]]][reorder]
    }
  }
  return(outData)
}


get.immunog.locs <- function(ucsc=c("hg18","hg19"),bioC=F,text=F) {
  nchr <- 22
  if(ucsc[1]=="hg19") {
    #hg19
    chr <- c(22,14,2,14)
    stz <- c(22385572,105994256,89156874,22090057)
    enz <- c(23265082,107281230,89630187,23021097)
  } else {
    # hg18
    chr <- c(22,14,2,14)
    stz <- c(20715572,105065301,88937989,21159897)
    enz <- c(21595082,106352275,89411302,22090937)
  }
  nmz <- c("ig_c22","ig_c14_a","ig_c2","ig_c14_b")
  reg.dat <- rep("immunoglobin",length(chr))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chr,
                          reg=reg.dat,universe=ucsc[1])
    outData <- toGenomeOrder(outData)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- vector("list",nchr); names(outData) <- paste("chr",1:nchr,sep="")
    for (cc in 1:nchr) {
      if(cc %in% chr) {
        outData[[cc]] <- list(start=stz[chr==cc],end=enz[chr==cc])
      }
    }
  }
  return(outData)
}



get.dgv.ranges <- function(dir=NULL,ucsc="hg18",bioC=T,text=F,shortenID=T, compact=F, alt.url=NULL)
{
  ## download or use local version of DGV (database for Genomic Variants)
  from.scr <- T
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  old.colnm.core <- c("VariationID","Start","End","Chr","Gain","Loss") # order: ids,st,en,chr,...
  colnm.core <- c("variantaccession","start","end","chr","observedgains","observedlosses") # order: ids,st,en,chr,...
  if(!is.null(dir)) {
    dg.fn <- cat.path(dir$ano,"dgvAnnot.RData")
    if(file.exists(dg.fn)) {
      tt <- get(paste(load(dg.fn)))
      from.scr <- F 
    }
  }
  if(from.scr) {
    #dgv.path <- paste("http://projects.tcag.ca/variation/downloads/variation.",ucsc[1],".v10.nov.2010.txt",sep="")
    dgv.path <- switch(ucsc,hg19="http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2013-05-31.txt",
    hg17="http://dgv.tcag.ca/dgv/docs/NCBI35_hg17_variants_2013-05-31.txt",
    hg18="http://dgv.tcag.ca/dgv/docs/NCBI36_hg18_variants_2013-05-31.txt")
    if(!is.null(alt.url)) { dgv.path <- alt.url }
    downloc <- cat.path(dir$ano,"DGV.variants.txt")
    download.file(url=dgv.path,downloc,quiet=T)
    tt <- read.delim(downloc,header=T)
    if(!any(colnames(tt) %in% c(colnm.core,old.colnm.core))) {
      warning("the URL or file format for the DGV seems to have changed. Go to http://dgv.tcag.ca/\n",
              "to find the latest url, and then run get.dgv.ranges() using alt.url to create a\n",
              "local copy of the DGV in your working dir\n")
      return(NULL)
    } else {
      # test whether the old or new names are in use
      oldL <- length(which(colnames(tt) %in% old.colnm.core))
      newL <- length(which(colnames(tt) %in% colnm.core))
      if(oldL>newL) { colnm.core <- old.colnm.core }
    }
    if(exists("dg.fn")) { save(tt,file=dg.fn) }
  }
  # shorten the default ID and make the source obvious!
  if(shortenID) {
    if(old.colnm.core[1] %in% colnames(tt)) { tt[,old.colnm.core[1]] <- gsub("Variation_","DGV_",tt[,old.colnm.core[1]]) } 
    if(colnm.core[1] %in% colnames(tt)) { 
      tt[,colnm.core[1]] <- paste("DGV",tt[,colnm.core[1]],sep="_") }   }
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    if(compact) { to.cut <- colnames(tt)[!colnames(tt) %in% colnm.core] } else { to.cut <- NULL }
    outData <- data.frame.to.ranged(tt,ids=colnm.core[1],start=colnm.core[2],end=colnm.core[3],
                                    width=NULL,chr=colnm.core[4],ucsc=ucsc, exclude=to.cut)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- tt 
  }
  matches <- which(colnames(outData) %in% colnm.core)
  indx <- match(colnames(outData)[matches],colnm.core)
  if(length(matches)>1) { colnames(outData)[matches] <- old.colnm.core[indx] }
  if(all(c("Gain","Loss","samplesize","variantsubtype") %in% colnames(outData))) {
    ff <- suppressWarnings(apply(((cbind(outData[["Gain"]],outData[["Loss"]]))),1,max,na.rm=T)/outData[["samplesize"]])
    ff[outData[["variantsubtype"]]=="Gain+Loss"] <- ff[outData[["variantsubtype"]]=="Gain+Loss"]/2
    ff[ff>1] <- NA
    outData[["frequency"]] <- ff
  }
  return(outData)
}


get.centromere.locs <- function(dir="",ucsc=c("hg18","hg19"),bioC=F,text=F)
{
  dir <- validate.dir.for(dir,c("ano"),warn=F); success <- T
  local.file <- paste(dir$ano,"cytoBand.txt.gz",sep="")
  if(!file.exists(local.file)) {
    golden.path <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",ucsc[1],"/database/cytoBand.txt.gz",sep="")
    success <- tryCatch( download.file(url=golden.path,local.file,quiet=T),error=function(e) { F } )
  }
  if(is.logical(success)) {
    if(!success) { warning("couldn't reach ucsc website! try sourcing cytoband data elsewhere"); return(NULL) } }
  tt <- read.table(local.file)
  nchr <- 22 #default
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",1:nchr,sep="")
  for (cc in 1:nchr) {
    just.centros <- tt[paste(tt[,5])=="acen",]
    just.chr <- just.centros[which(paste(just.centros[,1])==names(my.chr.range)[cc]),]
    my.chr.range[[cc]] <- list(start=min(just.chr[,2]), end=max(just.chr[,3]))
  }
  reg.dat <- rep("centromere",nchr)
  nmz <- paste(reg.dat,1:nchr,sep="_")
  stz <- sapply(my.chr.range,"[[",1)
  enz <- sapply(my.chr.range,"[[",2)
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=1:nchr,
                          reg=reg.dat,universe=ucsc[1])
    outData <- toGenomeOrder(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


get.telomere.locs <- function(dir="",kb=500,ucsc=c("hg18","hg19"),bioC=F,text=F)
{
  chr.lens <- get.chr.lens(dir=dir,ucsc=ucsc[1])
  nchr <- 22 #default
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",1:nchr,sep="")
  for (cc in 1:nchr) {
    one <- c(1,kb*1000)
    two <- chr.lens[[cc]]+c(-kb*1000,0)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  reg.dat <- rep("telomere",nchr*2)
  chrz <- rep(1:nchr,each=2)
  nmz <- paste(reg.dat,chrz,rep(c("a","b"),times=nchr),sep="_")
  stz <- as.vector(sapply(my.chr.range,"[[",1))
  enz <- as.vector(sapply(my.chr.range,"[[",2))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chrz,
                          reg=reg.dat,universe=ucsc[1])
    outData <- toGenomeOrder(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


get.chr.lens <- function(dir="",len.fn="humanChrLens.txt",ucsc=c("hg18","hg19")[1],n=c(1:22),delete.after=F)
{
  # retrieve chromosome lengths from local annotation file, else download from UCSC
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  chrlens.f <- cat.path(dir$ano,len.fn) # existing or future lengths file
  # backups for offline use
  if(ucsc=="hg18") {
    offline.backup <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
      146274826,140273252,135374737,134452384,132349534,114142980,106368585,
      100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432)
  } else {
    offline.backup <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
      146364022,141213431,135534747,135006516,133851895,115169878,107349540,
      102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  }
  if(file.exists(chrlens.f))
  {
    # file seems to be in annotation directory already
    chrLens <- readLines(chrlens.f)
    if (length(as.numeric(chrLens))!=length(n))
    {
      cat(paste("Length of chromosome file doesn't match expected",length(n),"\n"))
      notGot <- T
    } else { notGot <- F}
  } else { notGot <- T }
  if (notGot) {
    #download from UCSC
    cat("attempting to download chromosome lengths from UCSC\n")
    urL <- switch(ucsc,
                  hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
                  hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",)
    success <- T
    success <- tryCatch(download.file(urL, chrlens.f,quiet=T),error=function(e) { F } )
    if(success) {
      print("download succesfful")
      chrL.f <- readLines(chrlens.f)
      len.lst <- strsplit(chrL.f,"\t")
      nmz <- sapply(len.lst,"[",1)
      lnz <- sapply(len.lst,"[",2)
      want.chr.names <- match(paste("chr",n,sep=""),nmz)
      want.chr.names <- want.chr.names[!is.na(want.chr.names)]
      chrLens <- lnz[want.chr.names]
      names(chrLens) <- want.chr.names
    } else {
      warning("couldn't reach ucsc website, so have used offline versions of chr lengths")
      chrLens <- paste(offline.backup); names(chrLens) <- paste(1:22)
      delete.after <- F
    }
    if(length(dir)==1 & dir[1]=="" & delete.after) {
      unlink(chrlens.f)
    } else {
      writeLines(chrLens,con=chrlens.f) # save file for future use
    }
  }
  return(as.numeric(chrLens))
}



import.marker.data <- function(dir, markerinfo.fn="snpdata.map",snp.fn="snpNames.txt", anot=c("bim","vcf","map","map3")[4],
                               snp.col=NA, pos.col=NA, chr.col=NA, verbose=F )
{
  # calculate chromosome-wise position and ID for each SNP in list()
  # dir should be contain the annotation directory assumed to contain the snp list (dir.ano)
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  # read in gene annotation (vcf file/bim file) if not passed as fn. arg.
  # sort info by chromosome and position
  ######## READ FILE WITH CHR,POS,LABEL for SNPS ###########
  if(verbose) { cat("\nRetrieving SNP information (chr,pos,id) from",anot,"file\n") }
  markerinfo.fn <- cat.path(dir$ano,markerinfo.fn)
  if(file.exists(markerinfo.fn)) {
    del <- get.delim(markerinfo.fn)
    vcf.file <- switch(anot,
                       bim=read.delim(markerinfo.fn,comment.char="#",header=F,stringsAsFactors=F,sep=del),
                       vcf=read.delim(markerinfo.fn,comment.char="#",header=T,stringsAsFactors=F,sep=del),
                       map=read.delim(markerinfo.fn,header=F,stringsAsFactors=F,sep=del),
                       map3=read.delim(markerinfo.fn,header=F,stringsAsFactors=F,sep=del) )
    if(verbose | T) { cat(" file preview:\n"); print(head(vcf.file,4)); cat("\n") }
  } else {
    linez <- character()
    linez[1] <- (paste("Error: expecting file: '",markerinfo.fn,"' from function parameter 'markerinfo.fn'.",sep=""))
    linez[2] <- ("This file should be a vcf, bim, map, or map3 file (see plink website")
    linez[3] <- ("for description of map and map3 formats). 'Map3' is the simplest with")
    linez[4] <- ("chromosome number in column 1, snp id in column 2, snp position in column 3.")
    linez[5] <- ("Can be space, tab or comma delimited.")
    cat(linez,"\n"); stop()
  }
  # if no user value entered, set column with SNP labels to default for given mode
  if(is.na(snp.col)) { snp.col <- switch(anot,bim=2,vcf=2,map=2,map3=2) }
  # if no user value entered, set column with SNP positions to default for given mode
  if(is.na(pos.col)) { pos.col <- switch(anot,bim=4,vcf=4,map=4,map3=3) }
  if(is.na(chr.col)) { chr.col <- 1 } # (these main file types all have chr in col 1)
  if(verbose) {
    cat(" assuming columns are:  SNP-id:",snp.col,"; Chr:",chr.col,"; Pos:",pos.col,"\n")
    cat(" if this does not match file preview above, please stop and change file 'anot', or set")
    cat(" the values of:\n snp.col, pos.col, chr.col, explicitly in functions passing args to 'calc.chr.ind'.\n")
  }
  ######## READ FILE WITH SUBSET OF SNPS ###########
  if(is.null(snp.fn) | length(snp.fn)>10)
  {
    if(is.null(snp.fn))
    {
      if(verbose) { cat(" no subset of snps selected, default to use all in map/bim/vcf file\n") }
      snp.list <- vcf.file[,snp.col]
      match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
    } else {
      if(verbose) { cat(" using vector of snp IDs '",basename(snp.fn),"' as subset\n",sep="") }
      snp.list <- snp.fn
      match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
      pc.missing <- (length(which(is.na(match.list.to.vcf)))/length(snp.list))
      if(pc.missing>0) {
        cat("",paste(round(100*pc.missing,1),"% snps missing from annotation file\n"))
      }
    }
  } else {
    if(is.file(snp.fn,dir$ano,dir)) { snp.fn <- find.file(snp.fn,dir$ano,dir) }
    if(!file.exists(snp.fn)) 
    {
      linez <- character()
      linez[1] <- (paste("Error: expecting file: '",snp.fn,"' from parameter 'snp.fn' in ",dir$ano,sep=""))
      linez[2] <- ("This file should be a list of snp ids (one per line) to include in the current process")
      linez[3] <- ("This could be a list of all snps, or any subset of snps in the map/bim/vcf file.")
      linez[4] <- ("Alternatively pass in a character() list of ids, or NULL to include all in the map/bim/vcf file.")
      cat(linez,"\n") ; stop()
    } else {
      snp.list <- readLines(snp.fn)
    }
    # match snp list to vcf file
    match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
    pc.missing <- (length(which(is.na(match.list.to.vcf)))/length(snp.list))
    if(verbose) { cat("",paste(round(100*pc.missing,1),"% snps missing from annotation file\n")) }
    if(pc.missing>.5) { stop("too many missing. comment out [#] line in function 'calc.chr.ind' if intentional")}
  }
  snpData <- RangedData(ranges=IRanges(start=vcf.file[,pos.col],width=1,
                                       names=paste(vcf.file[,snp.col])),space=vcf.file[,chr.col])
  snpData <- toGenomeOrder(snpData,strict=T)
  return(snpData)
}


stats.on.CNV.file <- function(fn,use.dat=F) 
{
  ## return detailed summary of CNVs in a plink .cnv file
  dat <- read.table(fn,header=T)  
  tt <- as.numeric(table(dat$IID))[order(as.numeric(table(dat$IID)))]
  # headers in file
  #FID  IID  CHR  BP1	BP2	TYPE	SCORE	SITES
  lenz <- (dat$BP2 - dat$BP1)
  cat("Summary of CNVs in file:",fn,"\n")
  cat("\nSummary of CNV lengths (ALL)","\n\n")
  print(summary(lenz),digits=8) ; cat("\n")
  if(use.dat) {	do.lens.summary(lenz,dat) } else { do.lens.summary(lenz) }
  cat("\nlog histogram\n") 
  textogram(log(1+lenz))
  cat("\nSummary of CNV lengths (by chromosome)","\n")
  cat("Mean:\n");	print(tapply(lenz,factor(dat$CHR),mean),digits=2)
  cat("\nSummary of CNV lengths (by cnv type)","\n")
  cat("Mean:\n");	print(tapply(lenz,factor(dat$TYPE),mean),digits=2)
  cat("Range:\n");  print(tapply(lenz,factor(dat$TYPE),range),digits=2)
  denz <- lenz/dat$SITES
  #print(cbind(lenz,dat$SITES)[lenz<20,],"\n")
  cat("\nSummary of CNV densitys (ALL, length per SNP)","\n")
  print(summary(denz))
  cat("\nlog histogram\n") 
  textogram(log(1+denz))
  cat("\nSummary of CNV densitys (by chromosome)","\n")
  cat("Mean:\n");	print(tapply(denz,factor(dat$CHR),mean),digits=2)
  cat("\nSummary of CNV densitys (by cnv type)","\n")
  cat("Mean:\n");	print(tapply(denz,factor(dat$TYPE),mean),digits=2)
  cat("Range:\n"); print(tapply(denz,factor(dat$TYPE),range),digits=2)
  cat("\nSummary of number of SNPS per CNV (ALL)","\n")
  print(summary(dat$SITES)) 
  cat("\nlog histogram\n") 
  textogram(log(1+dat$SITES))
  cat("\nSummary of number of SNPS per CNV (by chromosome)","\n")
  cat("Mean:\n");	print(tapply(dat$SITES,factor(dat$CHR),mean),digits=2)
  cat("\nSummary of number of SNPS per CNV (by cnv type)","\n")
  cat("Mean:\n");	print(tapply(dat$SITES,factor(dat$TYPE),mean),digits=2)
  cat("Range:\n"); print(tapply(dat$SITES,factor(dat$TYPE),range),digits=2)
  #stats on cnvs
  cat("\nSummary of number of CNVs per Sample","\n")
  print(summary(tt),digits=2)
}


getSlot <- function(DT,what="description",grps=NA,n.pcs=NA,dir=NULL,
                    result.name=NULL,ret.obj=F,force.valid=T,warns=F)
{
  # retrieve a field from a data tracker object (using short labels), can also do for grp subset
  # list valid values for 'what' by type of file they refer to
  mdt <- methods.data.tracker(T)
  what <- what[1] # only first value used
  stdz <- function(x) { tolower(gsub(".","",x,fixed=T)) } # avoid mismatches due to simple error
  if(!stdz(what) %in% stdz(mdt$valid))
  {
    warning(paste("illegal attribute,",what,"requested, returning NULL"), 
            paste("\nlegal attributes are:",paste(mdt$valid,collapse=",")))
    return(NULL)
  }
  if(!is.data.tracker(DT)) {
    stop("Error: could not update because 'DT' was not a data.tracker object")
  }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,c("base"))
    if(dir$base!=DT$Info$BASE) {
      warning("directory object did not match existing. No change to data tracker object")
      return(DT)
    }
  } else {
    if(is.ch(DT$Info$DIR)) { dir <- DT$Info$DIR } 
  }
  # those with grps = # of files [note that multi in 1 group stored in a sublist]
  if(stdz(what)=="rawlrr")     { out <- DT$RawDataGrps$LRR }
  if(stdz(what)=="shapelrr")   { out <- DT$ShapedDataGrps$LRR }
  if(stdz(what) %in% c("geno","genofile"))  { out <- DT$SnpDataGrps$GENO }
  if(stdz(what)=="subsamples") { out <- DT$BigRawDataGrps$SAMPLES }
  if(stdz(what)=="rawbaf")     { out <- DT$RawDataGrps$BAF }
  if(stdz(what)=="shapebaf")   { out <- DT$ShapedDataGrps$BAF }
  # those with ugrps = # of files
  if(stdz(what)=="biglrr")     { out <- DT$BigRawDataGrps$LRR }
  if(stdz(what)=="bigfilt")    { out <- DT$BigFiltDataGrps$LRR }
  if(stdz(what)=="pcas")        { out <- DT$BigPCADataGrps$LRR }
  if(stdz(what)=="eigens")      { out <- DT$BigPCADataGrps$EIGEN }
  if(stdz(what)=="snpfields")  { out <- DT$SnpDataGrps$FIELDS }
  if(stdz(what)=="medians") { out <- DT$BigRawDataGrps$MEDIAN }
  if(stdz(what)=="stats") { out <- DT$BigRawDataGrps$STATS }
  if(all(!is.na(grps))) {
    #check grp based input
    exgrps <- DT$Info$N.GRPS
    nmgrps <- DT$Info$GRP.NAMES
    ugrps <- unique(grps)
    nugrps <- length(ugrps)
    if(all(ugrps %in% nmgrps)) { 
      grps <- match(grps,nmgrps); ugrps <- unique(grps)
    }
    if(max(ugrps)>exgrps) {
      stop("Error: illegal group number specified in 'grps', exceeds number of grps")
    }
    if(nugrps<exgrps) {
      if(exists("out") & warns) {
        warning("fewer grps specified than 'DT$Info$N.GRPS' so assuming retrieval applies to subset")
      }
    } else {
      if(nugrps!=exgrps) {
        stop(paste("Error: number of grps specified (",nugrps,
                   ") exceeds n.grps in 'DT' (",exgrps,")"),sep="")
      }
    }
    lluu <- length(ugrps)
    if(lluu>0 & exists("out") & stdz(what)!="snpfields") {
      out <- out[ugrps] 
    }
    #if(length(out)==lluu) {
    #  if(stdz(what)=="rawbaf")   { out <- DT$RawDataGrps$BAF[ugrps] }
    #  if(stdz(what)=="shapebaf") { out <-DT$ShapedDataGrps$BAF[ugrps] }
    #  try.grped.baf <- F
    #}
  } 
  # update combined (unique) inputs
  if(stdz(what)=="dir") { out <- DT$Info$DIR }
  if(stdz(what) %in% c("numgrps","ngrps","ngrp")) { out <- DT$Info$N.GRPS }
  if(stdz(what) %in% c("procdone","processdone","processesdone")) { out <- DT$Info$PROC.DONE }
  if(stdz(what) %in% c("numfiles","nfiles","nfile")) { out <- DT$Info$N.FILES }
  if(stdz(what) %in% c("groupnames","grpnames")) { out <- DT$Info$GRP.NAMES }
  if(stdz(what) %in% c("mapfile")) { out <- DT$Info$MAP.FILE }
  if(stdz(what) %in% c("maptype")) { out <- DT$Info$MAP.TYPE }
  if(stdz(what) %in% c("snpstats")) { out <- DT$SnpDataGrps$SNPMAT }
  if(stdz(what)=="description") { out <- DT$Info$DESCRIPTION }
  if(stdz(what)=="bigqc") { out <- DT$BigQCData$LRR }
  if(stdz(what)=="bigbaf") { out <- DT$BigQCData$BAF }
  if(stdz(what)=="gcgenome") { out <- DT$BigQCData$INFO$gc.genome }
  if(stdz(what)=="gcmarker") { out <- DT$BigQCData$INFO$gc.marker }
  if(stdz(what) %in% c("sampleinfo","sampleinfofn")) { out <- DT$BigQCData$INFO$sample.info }
  if(stdz(what) %in% c("snpinfo","snpinfofn")) { out <- DT$BigQCData$INFO$snp.info }
  if(stdz(what)=="snps") { out <- DT$ShapedDataGrps$IDs$snps }
  if(stdz(what)=="samples") { out <- DT$ShapedDataGrps$IDs$samples }
  # settings has 8 slots, return the last valid
  if(stdz(what)=="settings") { out <- DT$Info$SETTINGS } # return all settings
# if(stdz(what)=="settings") {
#     ttt <- sapply(DT$Info$SETTINGS,length)
#     last.non.empty <- tail(which[ttt>1],1)
#     if(length(ttt)>0) {
#       out <- DT$Info$SETTINGS[[ttt]]
#     } else {
#       out <- NULL
#     }
#   }
  # update named inputs
  n.pcs <- n.pcs[1]
  if(stdz(what) %in% c("bigpcc","eigen","median","stat","cnvresults","cnvresult","cnvcalls")) {
    if(!is.na(n.pcs)) {
      if(stdz(what)=="bigpcc") { out <- DT$BigCorrectedData$LRR[[paste("pc",n.pcs,sep="")]] }
      if(stdz(what)=="eigen") { out <- DT$BigCorrectedData$EIGEN[[paste("pc",n.pcs,sep="")]] }
      if(stdz(what)=="median") { out <- DT$BigCorrectedData$MEDIAN[[paste("pc",n.pcs,sep="")]] }
      if(stdz(what)=="stat") { out <- DT$BigCorrectedData$STATS[[paste("pc",n.pcs,sep="")]] }
      if(stdz(what) %in% c("cnvresults","cnvresult","cnvcalls")) { out <-DT$CnvResults$RESULT[[paste("pc",n.pcs,sep="")]] }
    } else {    
      if(warns) {  
        warning(paste("'n.pcs' not passed to getSlot(",what,"); default is to select the first entry",sep=""))
      }
      if(stdz(what)=="bigpcc") { out <- DT$BigCorrectedData$LRR[[1]] }
      if(stdz(what)=="eigen") { out <- DT$BigCorrectedData$EIGEN[[1]] }
      if(stdz(what)=="median") { out <- DT$BigCorrectedData$MEDIAN[[1]] }
      if(stdz(what)=="stat") { out <- DT$BigCorrectedData$STATS[[1]] }
      if(stdz(what) %in% c("cnvresults","cnvresult","cnvcalls")) { out <-DT$CnvResults$RESULT[[1]] }
    }
  }
  #result.name <- result.name[1]
  #if(is.character(result.name)) {
  #  if(stdz(what)=="cnvresult") { out <-DT$CnvResultData$RESULT[[result.name]] }
  #}
  if(!exists("out")) { out <- NULL }
  if(length(out)>0 & stdz(what)!="snpfields") {
    if(force.valid & (stdz(what) %in% stdz(mdt$filenames))) {
      out <- find.file(out,dir) #add directory if needed
    }
    if(stdz(what) %in% stdz(mdt$numbers)) {
      if(!is.numeric(out)) {
        if(length(narm(as.numeric(out)))==length(out)) {
          out <- as.numeric(out) # in case stored as text
        }
      }
    }
    if(ret.obj) {
      if(is.file(out,dir)) {
        # if using this option, the object rather than just the filename is requested
        if(stdz(what) %in% stdz(mdt$dats)) {
          out <- lapply(out,force.frame)
        }
        if(stdz(what) %in% stdz(mdt$ranged)) {
          out <- lapply(unlist(out),function(fn) { get(paste(load(fn))) })
        }
        if(stdz(what) %in% stdz(mdt$bigs) & length(out[[1]])==1) {
          out <- getBigMat(out[[1]][1],dir)  
        }
        if(stdz(what) %in% stdz(mdt$massives)) {
          warning("requested object is too large to retrieve, returning filename only")
        }
        if(stdz(what) %in% stdz(mdt$vecs)) {
          #out <- lapply(out,force.vec)
          #print(headl(out))
          if(is.list(out)) {
            for (jj in 1:length(out)){
              if(length(unlist(out[[jj]]))<1000) {
                ## if out is big can waste lots of time checking whether each is a file
                if(is.file(out[[jj]],dir)) {
                  out[[jj]] <- sapply(find.file(out[[jj]],dir),function(x) { lapply(x,readLines) } )
                }
              }
            }
          } else {
            out <- lapply(out,readLines)
          }
        }
      } else {
        warning(paste("filename",out,"does not exist"))
        return(NULL)
      }
    }
    return(out)
  } else {
    if(length(out)==0) { warning("output of data tracker lookup is empty") }
    return(out)
  }
}


is.data.tracker <- function(DT) 
{
  # check whether object is in correct 'data.tracker' format  
  is.it <- F
  if(is.null(DT)) { return(F) }
  if(is.list(DT)) {
    if(!is.null(names(DT))) {
      if(all(names(DT)==c("Info","RawDataGrps", "SnpDataGrps", "ShapedDataGrps", "BigRawDataGrps",
                          "BigFiltDataGrps","BigPCADataGrps","BigQCData", "BigCorrectedData",
                          "CnvResults") )) { is.it <- T } } }
  return(is.it)
}


make.snp.info <- function(dir=NULL, map=NULL, snpIDs=NULL, anot="map3", ucsc="hg18",
                          snp.info.fn="snpdata.map",sav=F,verbose=F,filt.by=NULL,make.sort=F,
                          extended=T,non.autosomes.excl.file=F,autosomes.only=F,absolute.index=T,
                          scheme=NA,col1="black",col2="grey",...)
{
  # make a snp.info object from scratch by combining a list of snps (in order of datafile)
  # with an annotation file (map3/bim/vcf/etc), optionally filtering for only 'filt.by' snps
  if(is.null(dir)) { if(is.null(map)) { warning("dir not specified, using current") } ; dir <- getwd() }
  dir <- validate.dir.for(dir,c("ano","excl2","sort2"))
  if(!is.null(map)) { 
    id.names <- c("Id","snp.name","snp","snpid","snp.id",
    "label","rsid","rs.id","marker","marker.name","probe","probeid","probe.id")
    sobj <- map
    # coerce chr,pos,id columns to those names if possible
    sobj <- column.salvage(sobj,"chr",c("Chr","Chromosome","space"),ignore.case=T)
    sobj <- column.salvage(sobj,"pos",c("Pos","position","location"),ignore.case=T)
    if(!any(id.names %in% colnames(sobj))) { sobj[["id"]] <- rownames(sobj) }
    sobj <- column.salvage(sobj,"id",id.names,ignore.case=T)
    want <- c("chr","id","pos"); if(!all(want %in% colnames(sobj))) { 
      warning("invalid map object, need chr, id, pos columns"); return(NULL) }
    write.table(sobj[,want],cat.path(dir$ano,snp.info.fn),sep="\t",col.names=F,row.names=F)
  }
  snp.info.fn <- find.file(snp.info.fn,dir$ano,dir) 
  #CHR.INFO <- calc.chr.ind(dir=dir, snp.fn=snpIDs, anot=anot, 
  #                         markerinfo.fn=snp.info.fn, sav=F, verbose=verbose,...)
  if(!is.null(snpIDs)) { snpnames <- force.vec(cat.path(dir$ano,snpIDs)) } else { snpnames <- NULL }
  
  snp.info <- get.chr.info.filt(dir=dir, extended=T, filt.names=filt.by,verbose=verbose,
                                snp.info.fn=snp.info.fn,ucsc=ucsc,autosomes.only=autosomes.only,
                                absolute=absolute.index,snp.fn=snpnames,anot=anot,scheme=scheme,col1=col1,col2=col2,...)
  if(make.sort) {
    # if not present will make the snpsort.txt file; needed to initialise the start of the pipeline
    #print("tere")
    null.thing <- choose.best.sort.file(dir$sort2,ref.list=rownames(snp.info))
  }  
  if(non.autosomes.excl.file) {
    if(is.null(snpnames)) { snpnames <- rownames(snp.info) }
    chr23.removed <- paste(snpnames[!snpnames %in% rownames(select.autosomes(snp.info))])
    if(length(chr23.removed)>0) {
      anot.file <- cat.path(dir$excl2,"nonAutosomes.txt")
      writeLines(chr23.removed,con=anot.file)
      cat("~wrote excluded non-autosome snp list to:",dir$excl2,"\n")
    }
  }
  return(snp.info)
}


check.data.size <- function(nsamp,nsnp) {
  maxsize <- ((2^31)-2)
  if(!is.numeric(nsamp) | !is.numeric(nsnp)) { stop("Error: check.data.size takes numeric args")}
  size <- nsamp*nsnp; pc.max <- round(((size/2^31)*100),1)
  valid <- size<maxsize
  if(any(!valid)) 
  {
    warning(paste("*warning:",(nsamp*(nsnp/nsnp))[!valid],"x",(nsnp*(nsamp/nsamp))[!valid],"is",
              pc.max[!valid],"% of the allowed object size for SnpMatrix\n"),
    "Please run the SNP-QC in plink or split the import files into smaller chunks")
  } 
  return(valid)
}


read.penn.cnv.file <- function(filename,readtable=T) {
  cN <- c("coords","n.snps","length","copy.number","file","first.snp","last.snp")
  if(readtable) {
    file.out <- read.table(filename,header=F)
    colnames(file.out) <- cN
  } else {
   raw.dat <- readLines(filename)
   file.out <- as.data.frame(matrix(nrow=length(raw.dat),ncol=length(cN)))
   colnames(file.out) <- cN
   splitto <- strsplit(raw.dat,"\t| +")
   file.out[,] <- t(sapply(splitto,c))
  }
  rownames(file.out) <- NULL
  return(file.out)
}


read.plink.file <- function(filename,readtable=T) {
  cN <- c("FID","IID","CHR","BP1","BP2","TYPE","SCORE","SITES")
  if(readtable) {
    file.out <- read.table(filename,header=T)
    if(any(!colnames(file.out) %in% cN)) { 
      warning("column names [",colnames(fileout),"]did not match expected")  
      colnames(file.out) <- cN
    }
  } else {
    raw.dat <- readLines(filename)
    file.out <- as.data.frame(matrix(nrow=length(raw.dat)-1,ncol=length(cN)))
    colnames(file.out) <- cN
    splitto <- strsplit(raw.dat,"\t| +")
    nlines <- length(splitto); sN <- 100
    smpl <- unique(c(1:sN,sample(nlines,sN,T),((nlines-sN):nlines)))
    smpl <- smpl[smpl %in% 1:nlines] # select first, random and last hundred
    nC <- length(cN)
    if(all(as.numeric(names(table(sapply(splitto[smpl],length))))==c(nC,nC+1))) {
      # everything seems to be as expected
      if(!all(cN %in% splitto[[1]])) { 
        warning("didn't find expected headers [",cN,"] in plink file, parse may fail") }
      splitto <- splitto[-1]
      for (cc in 1:nC) {
        file.out[,cc] <- (sapply(splitto,"[",cc))
      }
      rownames(file.out) <- NULL
    } else {
      # unfortunately file type seems to have changed
      warning("plink file format has not imported as expected, using slower method to parse")
      file.out <- read.table(filename)
    }
  }
  return(file.out)
}


rmv.dir.plink.file <- function(...)
{
  rmv.dir.penn.cnv.file(...,plink=T)
}


rmv.dir.penn.cnv.file <- function(filename,append=".nodir",ext=F,verbose=F,
                                  plink=F,readtable=T) {
  # whether to read files using 'read.table()' or manual method.. speed varies
  raw.file <- readLines(filename)
  if(plink) {
    p.file <- read.plink.file(filename,readtable=readtable); cL <- "IID"
  } else {
    p.file <- read.penn.cnv.file(filename,readtable=readtable); cL <- "file"
  }
  all.dirs <- dirname(p.file[,cL])
  un.all <- unique(all.dirs)
  un.dirs <- dir.force.slash(un.all)
  if(length(un.all)>0) {
    for (cc in 1:length(un.dirs)) {
      if(verbose) { cat(" removing path text:",un.dirs[cc],"\n") }
      raw.file <- gsub(un.dirs[cc],"",raw.file)
    }
  }
  if(ext) {
    all.ext <- get.ext(p.file[,cL])
    un.all <- unique(all.ext)
    if(length(un.all)>0) {
      un.ext <- paste(".",un.all,sep="")
      for (cc in 1:length(un.ext)) {
        if(verbose) { cat(" removing extension text:",un.ext[cc],"\n") }
        raw.file <- gsub(un.ext[cc]," ",raw.file)
      }
    }
  }
  ofn <- paste(filename,append,sep="")
  writeLines(raw.file,con=ofn)
  if(verbose | T) { cat("~wrote file",ofn,"with unwanted path text removed\n") }
  return(ofn)
}


q.cmd <- function(file.name,output.dir="",id="") {
  return(paste("qsub -q ",id," -o ",output.dir," -j y ",file.name,sep=""))
}


convert.textpos.to.data <- function(text) {
  # convert text locations like chr5:12344-5253553 into a matrix with cols: chr, start, end
  do.one <- function(X) {
    chr.pos <- strsplit(X,":",fixed=T)[[1]]  
    chr <- as.integer(gsub("chr","",chr.pos[[1]],ignore.case=T))
    pos12 <- as.integer(strsplit(chr.pos[[2]],"-",fixed=T)[[1]])
    out <- c(chr,pos12[1],pos12[2])
    names(out) <- c("chr","start","end")
    return(out)
  }
  return(t(sapply(text,do.one)))
}


convert.penn.to.plink <- function(penn.in,plink.out=NULL,fixed.width=F) {
  ## convert from penn cnv format to plink format
  ## plink has tab delimited and fixed width formats. can do either
  penn <- read.penn.cnv.file(penn.in)
  if(ncol(penn)<1 | nrow(penn)<1) { stop("file ",penn.in," was empty, returning NULL") }
  colnames(penn) <- c("location","numsnp","length","cn","id","startsnp","endsnp")
  loc <- as.data.frame(convert.textpos.to.data(penn$location))
  plink <- data.frame(FID=1,IID=penn$id,CHR=loc$chr,BP1=loc$start,BP2=loc$end,
                      TYPE=extract.val.penn(penn$cn),SCORE=0,
                      SITES=extract.val.penn(penn$numsnp))
  if(!is.null(plink.out)) {
    if(fixed.width) {
      plink <- rbind(colnames(plink),plink)
      plink <- conv.fixed.width(plink)
      write.table(plink,file=plink.out,row.names=F,col.names=F,quote=F)
    } else {
      write.table(plink,file=plink.out,row.names=F,quote=F,sep="\t")
    }
  } else {
    return(plink)
  }
}


make.bad.region.file <- function(dir,ucsc="hg18",fn="badRegions.txt") {
  dir <- validate.dir.for(dir,"cnv.plk")
  excl.reg <- c(get.telomere.locs(text=T,ucsc=ucsc,kb=500),
                get.centromere.locs(text=T,ucsc=ucsc),
                get.immunog.locs(text=T,ucsc=ucsc))  #IG must be last for ordered chr
  writeLines(paste(sort(excl.reg)),con=cat.path(dir$cnv.qc,fn))
}


add.info.to.ranges <- function(dir,cnv.ranges,target="phenotype",info.file="pheno.lookup.txt") {
  ## add information to a RangedData object (like phenotype)
  min.chars <- 4 # column name in 'file' must have at least the first 'min.chars' letters of 'target' 
  if(is(cnv.ranges)[1]!="RangedData") { warning("must be a RangedData object"); return(cnv.ranges) }
  if(!"id" %in% colnames(cnv.ranges)) { warning("cnv.ranges must have column 'id'"); return(cnv.ranges) }
  phenot <- as.data.frame(reader(fn=find.file(info.file,dir$ano,dir)))
  id.col <- find.id.col(phenot,ids=cnv.ranges$id,ret="col")
  if(id.col==0) { idz <- rownames(phenot) } else { idz <- phenot[,id.col] }
  if(ncol(phenot)==1) { phenoz <- phenot[,1] } else {
    colz <- which(tolower(colnames(phenot)) %in% substr(rep(tolower(target),
                                                            ((nchar(target)-min.chars)+1)), 1,(min.chars:nchar(target))) )
    if(length(colz)>0) { phenoz <- phenot[,colz[1]] } else { phenoz <- NULL }
  }
  if(!is.null(phenoz)) {
    row.select <- match(cnv.ranges$id,idz)
    cnv.ranges[[target]] <- phenoz[row.select]
  } else {
    warning(paste("found ",target," file but couldn't find column named '",target,"' so skipped",sep=""))
  }
  return(cnv.ranges)
}


plot.pheno.cnvs <- function(fn,type="DEL",pref="",dir)
{
  ## make plot of CNV regions across chromosomes in competing phenotypes for rare dels/dups
  if(!file.exists(fn)) { warning("cnv.summary file not found"); return(NULL) }
  if(length(grep("cnv.summary",fn))<1) { warning("doesn't look like a plink cnv.summary file: may fail") }
  tt <- reader(fn,treatas="txt")
  chr.list <- unique(tt[,1])
  if(length(chr.list)>0) {
    if(pref!="") { pref <- paste(pref,"_",sep="") }
    ofn <- cat.path(dir$res,fn="freqPlotPerCNVRegion",suf=type,ext="pdf",pref=pref)
    pdf(ofn)
    linev <- F
    allx <- vector("list",length(chr.list)); names(allx) <- paste(chr.list)
    for (cc in 1:length(chr.list)) {
      xx <- tt[tt[,1]==chr.list[cc],-c(1:2)]
      xx <- apply(xx,2,as.numeric); allx[[paste(chr.list[cc])]] <- xx
      cnt <- xx[,3]; 
      if(linev) {
        plot(xx[,1][cnt!=0]/10^6,cnt[cnt!=0],col="blue", ylim=range(xx[,2:3],na.rm=T)*c(1,1.1),type="l",
             ylab=paste(type,"count"),xlab="position (mb)",main=paste("chromosome",chr.list[cc]))
      } else {
        plot(xx[,1][cnt!=0]/10^6,cnt[cnt!=0],col="blue", ylim=range(xx[,2:3],na.rm=T)*c(1,1.1),pch="-",cex=0.75,
             ylab=paste(type,"count"),xlab="position (mb)",main=paste("chromosome",chr.list[cc]))
      }
      cnt <- xx[,2] ; 
      if(linev) {
        lines(xx[,1][cnt!=0]/10^6,cnt[cnt!=0]+.1,col="red") 
      } else {
        points(xx[,1][cnt!=0]/10^6,cnt[cnt!=0]+.1,col="red",pch="-",cex=0.75)
      }
      legend("top",legend=c("affected","unaffected"),col=c("red","blue"),pch="-",bty="n",ncol=2,pt.cex=2)
    }
    dev.off()
    cat("~wrote plot:",ofn,"\n")
  }
  return(allx)
}


init.dirs.fn <- function(dir,overwrite=F,ignore=c("raw","sup"),
                         silent=F,info.dir=NULL,update.bash=F,file.spec="file.spec.txt",make.snp.samp.files=T,
                         plate.info="plate.lookup.txt",pheno="pheno.lookup.txt",sexfile="sex.lookup.txt",
                         bash.file="getDataGS.sh",snp.info.fn="snpdata.map",id.files=NULL)
{
  # create the directory structure implied by the object 'dir'
  # provides options to delete existing directories if overwrite=T
  if(!is.list(dir)) { failnow <- T } else { failnow <- F }
  if(!all(sapply(lapply(dir,is),"[",1) %in% "character")) { failnow <- T }
  if(failnow) { stop("object 'dir' should be a list of directory locations") }
  ## don't try to create 'ignore' directories, ie., might be read-only, etc.
  if(!is.null(ignore)) {
    ignore <- ignore[ignore %in% names(dir)]
    for (dd in 1:length(ignore)) {
      dir[[paste(ignore[dd])]] <- NULL
    }
  }
  ## sort list by number of 'slash's, as this will mean that the
  # higher level directories will already exist when creating sub dirs
  countslash <- function(text) { length(text[text=="/"]) } 
  out1 <- sapply(dir,strsplit,split="")
  slashcounts <- sapply(out1,countslash)
  dirs.ordered <- unlist(dir)[order(slashcounts[match(names(dir),names(slashcounts))])]
  ldd <- length(dirs.ordered)
  passfail <- logical(ldd)
  ## iterate through creating each directory, taking note of whether it already exists
  # option to keep or erase contents if it does (eg, to create a fresh setup)
  for (cc in 1:ldd) {
    next.dir <- dirs.ordered[cc]
    if(!file.exists(next.dir)) {
      passfail[cc] <- dir.create(next.dir)
    } else {
      if(!silent) { cat("directory ",dir.force.slash(getwd()),next.dir,"already existed, ",sep="") }
      if(overwrite) {
        cat("contains",length(list.files(next.dir)),"files and/or sub-directories\n")
        choicez <- c("I am sure - DELETE","DO NOT delete")
        choice <- select.list(choicez,preselect=choicez[2],title=paste("delete contents of",next.dir,"?"))
        if(choice==choicez[1]) {
          cat("attempt to overwrite...")
          unlink(next.dir,recursive=T)
          passfail[cc] <- dir.create(next.dir)
          cat(c("succeeded","failed")[2-as.numeric(passfail[cc])],"\n")
        } else {
          cat("Directory deletion cancelled by user\n")
        }
      } else {
        if(!silent) { cat("was left unmodified\n") }
      }
    }
  }
  if(!silent) {
    cat(length(which(passfail)),"of",length(passfail),"directories created successfully\n") }
  if(update.bash) 
  {  
    scr.file <- cat.path(dir$scr,bash.file)
    if(!file.exists(scr.file)) {
      src <- cat.path(info.dir,bash.file)
      if(file.exists(src)) {
        file.copy(from=src,to=scr.file,recursive=T,overwrite=T)
      } else {
        gs <- gs.bash.http()
        if(!is.null(gs)) {
          # GITHUB WAY:
          writeLines(gs,file=scr.file)
          cat(" copied file",basename(scr.file),"\ninto:",dir$scr,"\nfrom github/plumbCNV/ \n")
        } else {
          ## RFORGE WAY:
          warning("Could not reach github, reverting to RForge version of the script which might be outdated")
          rforge.url <- "http://r-forge.r-project.org/scm/viewvc.php/*checkout*/scripts/getDataGS.sh?revision=2&root=plumbcnv"
          download.file(url=rforge.url,destfile=scr.file)
          cat(" copied file",basename(scr.file),"\ninto:",dir$scr,"\nfrom Rforge/plumbcnv/ \n")
        }
      }
      system(paste("chmod +x",scr.file))
    }
  }
  copy.from.aux <- function(src,targ) {
    if(file.exists(src)) {
      if(!file.exists(targ)) {
        file.copy(from=src, to=targ, overwrite = F, recursive = F, copy.mode = T)
        cat(" copied file",basename(src),"\ninto:",dirname(targ),"\nfrom:",dirname(src),"\n")
    } }
  }
  if(!is.null(info.dir)) {
    if(file.spec!="") {
      o_spcf <- cat.path(info.dir,file.spec)
      n_spcf <- cat.path(dir$lrr.dat,file.spec)
      copy.from.aux(o_spcf,n_spcf)
    }
    if(plate.info!="") {
      o_spcf <- cat.path(info.dir,plate.info)
      n_spcf <- cat.path(dir$ano,plate.info)
      copy.from.aux(o_spcf,n_spcf)
    }
    if(pheno!=""){
      o_spcf <- cat.path(info.dir,pheno)
      n_spcf <- cat.path(dir$ano,pheno)
      copy.from.aux(o_spcf,n_spcf)
    }
    if(sexfile!=""){
      o_spcf <- cat.path(info.dir,sexfile)
      n_spcf <- cat.path(dir$ano,sexfile)
      copy.from.aux(o_spcf,n_spcf)
    }
    if(snp.info.fn!=""){
      o_spcf <- cat.path(info.dir,snp.info.fn)
      n_spcf <- cat.path(dir$ano,snp.info.fn)
      copy.from.aux(o_spcf,n_spcf)
    }
    o_spcf <- cat.path(info.dir,"snpNames.txt")
    n_spcf <- cat.path(dir$ano,"snpNames.txt")
    if(file.exists(o_spcf) & !file.exists(n_spcf)) { 
      copy.from.aux(o_spcf,n_spcf) 
    } else {
      if(!file.exists(n_spcf)) {
        warning("couldn't find file with snp IDs, please create 'snpNames.txt' in",dir$ano)
      }
    }
    o_spcf <- cat.path(info.dir,"subIdsALL.txt")
    n_spcf <- cat.path(dir$ano,"subIdsALL.txt")
    if(file.exists(o_spcf) & !file.exists(n_spcf)) { 
      copy.from.aux(o_spcf,n_spcf) 
    } 
    if(!file.exists(n_spcf)) {
      if(make.snp.samp.files) {
        plinf <- cat.path(info.dir,plate.info)
        if(!file.exists(plinf)) { 
          plinf <- cat.path(info.dir,pheno)
        }
        if(file.exists(plinf)) {
          pl <- reader(plinf)
          rn <- rownames(pl)
          if(!is.null(rn)) {
            if(length(rn)==length(unique(rn))) {
              idz <- rn; writeLines(paste(idz),con=n_spcf)
            }
          } else {
            rn <- pl[,1]
            if(length(rn)==length(unique(rn))) {        
              idz <- rn; writeLines(paste(idz),con=n_spcf)
            }
          }
        } else {
          if(!file.exists(n_spcf)) {
            warning("couldn't find file with sample IDs, please create 'subIdsALL.txt' in",dir$ano) }
        }
      }
    }
    if(is.null(id.files)){
      o_spcf <- cat.path(dir$ano,"subIdsALL.txt")
      n_spcf <- cat.path(dir$ids,"subIdsALL.txt")
      if(length(list.files(dir$ids))<1) {
        copy.from.aux(o_spcf,n_spcf)
        warning("no separate ID files found/specified, so assuming 1 group, copied to:\n",dir$ids)
      }
    } else {
      o_spcf <- cat.path(info.dir,id.files)
      n_spcf <- cat.path(dir$ids,id.files)[1:length(o_spcf)]
      for (cc in 1:length(o_spcf)) {  copy.from.aux(o_spcf[cc],n_spcf[cc])  }
    }
  }
  return(dirs.ordered)
}


check.readiness <- function(dir=NULL,mode=2,snp.mode=1,penn.check=T,plink.check=T,
                            penn.path="/usr/local/bin/penncnv/",verbose=T,hmm="hh550.hmm") {
  ## check whether all appropriate files are in place to run plumbCNV
  ## Check for appropriate Penn CNV installation
  # mode 1: run the pipeline starting from a raw genome studio file, using file.spec.txt
  # mode 2: run the pipeline starting from pre-formatted LRR/BAF longfiles or matrices
  # mode 3: run the pipeline starting with big matrix object(s)
  # snp mode 1: run SNP QC using SNP STATS, must have genotype data
  # snp mode 2: skip SNP QC step (e.g, may have already excluded low quality snps)
  # snp mode 3: run SNP QC using plink (if from scratch) or have already run SNP-QC in plink manually
  load.all.libs(more.bio=c("snpStats","irlba"))
  
  if(.Platform$OS.type=="windows") { cat("\n\nWarning!: plumbCNV is not tested on windows\n");
                                     cat(" software such as PennCNV and Plink are required\n")
                                     cat(" as are bash commands. Pipeline may not succeed.\n\n")
  }
  if(plink.check) {
    try(system("plink -h > temp.txt ",intern=F))
    plink.test <- system("head temp.txt",intern=T)
    system("rm temp.txt")
    if(!length(plink.test)>0 | !is.character(plink.test)) {
      cat("Plink does not seem to be installed. If it is, please modify path so 'plink -h'",
          "shows help text\nMost steps in plumbCNV can be run without plink, but watch",
          "for errors in the final CNV\n",
          "calling if not installed. The PennCNV file 'raw.cnv' should still be usable\n")
    } else {
      cat("\nPlink installation detected\n\n"); #print(plink.test); cat("\n")
    }
  }
  if(penn.check) {
    penurl <- "http://www.openbioinformatics.org/penncnv/penncnv_download.html"
    if(!file.exists(penn.path)) {
      cat("Error: Make sure this open-source perl script is correctly installed to:",penn.path,"\n")
      cat("To install:")
      cat("(1) unzip the source tar.gz file \n(2) cd into the kext/ directory \n")
      cat("(3) optionally edit the “Makefile” to change options  (4) type 'make'\n")
      cat("(5) If there is no error message, the installation is done!\n")
      cat("Otherwise lperl might need to be installed by an admin;\n sudo apt-get install libperl-dev\n")
      cat("[NB: Standard Perl must also be installed, which should be fairly standard on linux systems]\n")
      stop(paste("Find PennCNV at:",penurl))
    }
    scn.scrpt <- "scan_region.pl"; cln.scrpt <- "clean_cnv.pl"
    conv.scrpt <- "penncnv_to_plink.pl"; force.scripts <- c(scn.scrpt,cln.scrpt)
    if(!file.exists(penn.path)) { stop("Error: PennCNV path incorrect or software not installed") }
    penn.available <- list.files(penn.path)
    if(!conv.scrpt %in% penn.available) { 
      scr.url <- "http://www.openbioinformatics.org/penncnv/download/penncnv_to_plink.pl"
      res <- download.file(scr.url,cat.path(penn.path,conv.scrpt),quiet=F)
    }
    hmm.file <- paste(penn.path,"lib/",hmm,sep="")
    if(!file.exists(hmm.file)) 
    { 
      stop(paste("expecting to find file",hmm.file,". Suggest reinstall PennCNV from",penurl)) 
    } else {
      cat("Found hidden markov parameter file lib/",hmm,".\n",sep="")
      cat("Modify parameters therein; B1_uf, B2_uf, B3_uf, from 0.01 to 0.03,")
      cat("if using Affymetrix arrays\n")
    }
    if(!check.linux.install(c("perl"))) {
      stop("perl installation not detected - please install perl\n")
    }
    cat("\nPennCNV installation detected\n")
  }
  # check dir
  if(is.null(dir)) { stop("Error: invalid directory 'dir' object") } else {
    dir <- validate.dir.for(dir,names(make.dir()),warn=T) # checks dir object is complete
  }
  #check for existence of compulsory files depending on where we plan to start from
  modez <- c("scratch","normal","big")
  snp.modez <- c("normal","skip","plink")
  if(is.character(mode)) { mode <- which(modez %in% tolower(mode)) }
  if(is.character(snp.mode)) { snp.mode <- which(snp.modez %in% tolower(snp.mode)) }
  if(!length(mode)==1 | !is.numeric(mode)) { mode <- 2 ; warning("invalid run mode, set to 2") }
  if(!length(snp.mode)==1 | !is.numeric(snp.mode)) { snp.mode <- 1 ; warning("invalid snp.mode, set to 1")}
  mode.lines <- c("Mode 1: run the pipeline starting from a raw genome studio file, using file.spec.txt",
                  "Mode 2: run the pipeline starting from pre-formatted LRR/BAF longfiles or matrices",
                  "Mode 3: run the pipeline starting with big matrix object(s)")
  snp.mode.lines <- c("Snp mode 1: run SNP-QC using snpStats package [NB: must have genotype data]",
                      "Snp mode 2: skip SNP-QC step (e.g, may have already excluded low quality snps)", 
                      "Snp mode 3: run SNP-QC using plink (if from scratch) or have already run SNP-QC in plink manually")
  # check standard linux commands used. only warn if not using bash script in scratch mode
  if(any(!check.linux.install(c("sed","getopts","grep","zcat","seq","awk","cat")))) {
    if(mode==1) {
      stop("Error: standard bash commands missing necessary when mode='scratch' - please install any with warnings above\n")
    } else {
      warning("standard bash commands missing: should be ok because 'mode' is not 'scratch'")      
    }
  }
  ## loop through validating required files depending on mode
  if(verbose) {
    cat("\nValidating required files for plumbCNV assuming:\n *    ",
        mode.lines[mode],"\n *",snp.mode.lines[snp.mode],"\n\n")
  }
  modez <- modez[mode]
  snp.modez <- snp.modez[snp.mode]
  # select appropriate options depending on mode chosen
  mode.select <- list(c(1,2,3),c(1,2,4,5,6,7,8),c(1,2,4,5,6,9))[[mode]]
  if(snp.modez=="plink" & modez!="scratch") { mode.select <- c(mode.select,10,11) }
  must.have.dirs <- c("ano","ano","lrr.dat","ano","ano","ids","col",
                      "baf.col","big","cr.plk","cr.plk")[mode.select]
  id.file <- list.files(dir$ids)[1]; if(is.na(id.file)) { id.file <- "<id_files!>" }
  lrr.file <- list.files(dir$col)[1]; if(is.na(lrr.file)) { lrr.file <- "<lrr_files!>" }
  baf.file <- list.files(dir$baf.col)[1]; if(is.na(baf.file)) { baf.file <- "<baf_files!>" }
  big.file <- list.files(dir$big)[1]; if(is.na(big.file)) { big.file <- "<big_files!>" }
  must.have.files <- c("plate.lookup.txt","pheno.lookup.txt","file.spec.txt",
                       "snpNames.txt","subIdsALL.txt",
                       id.file,lrr.file,baf.file,big.file,
                       "snpdata.irem", "snpdataout.imiss")[mode.select]
  must.have.descr <- c("sample IDs and plate IDs","sample IDs and phenotype code(s)",
                       "raw file type, group and column specifications",
                       "list of all snp IDs","list of all sample IDs",
                       "sample ids for each raw data file",
                       "LRR data in long or matrix format",
                       "BAF data in long or matrix format",
                       "bigMatrix description and backing files",
                       rep("plink snp-qc callrate results",2))[mode.select]
  ll <- min(c(length(must.have.dirs),length(must.have.files)))
  for (cc in 1:ll) {
    if(!is.file(must.have.files[cc],dir[[must.have.dirs[cc]]])) {
      cat("To run plumbCNV you should have file:",must.have.files[cc],"at path",dir[[paste(must.have.dirs[cc])]],"\n")
      if(is.file(must.have.files[cc],dir)) {
        warning(paste("it seems the file:",must.have.files[cc],"exists but is in the wrong directory: errors may ensue"))
      } else {
        stop(paste("Please ensure file(s) containing",must.have.descr[cc],"exists and restart the pipeline"))
      }
    }
  }
  # use auto mode if a valid 'file.spec.txt' file is found in the expected place dir$lrr
  if(is.null(get.file.specs(dir))) { auto.mode(F,set=T) } else { auto.mode(T,set=T) }
  if(auto.mode() & verbose) { cat("Found 'file.spec.txt', so will attempt to use information from",
                                  "this file to adjust import settings\n")}
  return(T)
}


#### SIMULATION FUNCTIONS

spike.in.cnv <- function(chr,pos,sample="",dir=NULL,big.lrr=NULL,big.baf=NULL,snp.info=NULL,add=T,type=1,noise=1) {
  ## add a simulated CNV to an existing dataset (bigMatrix)
  if(!is.null(dir)) { 
    dir <- validate.dir.for(dir,c("big","ano"))
    DT <- read.data.tracker(dir)
    snp.info <- read.snp.info(dir)
    if(!is.null(DT)) {
      if(is.null(big.lrr)) { big.lrr <- "big.lrr" }
      if(is.null(big.baf)) { big.baf <- "big.baf" }
      lst <- get.tracker.for.plot(DT=DT,dir=dir,PREPOSTPC=F,baf.file=big.baf,lrr.file=big.lrr,
                                  LRR=T,BAF=T,samples=sample)
      big.lrr <- lst$lrr.file; big.baf <- lst$baf.file
    } else {
      warning("invalid data tracker or 'dir'")
    }
  } else {
    if(is.null(big.lrr) | is.null(big.baf) | is.null(snp.info)) {
      stop("if dir is empty, then big.lrr, big.baf and snp.info arguments must be provided")
    }
  }
  if(is.null(big.lrr) | is.null(big.baf)) { stop("must have valid big.lrr and big.baf big.matrix objects or 'dir' which contains these")}
  if(is(snp.info)[1]!="RangedData") { stop("snp.info must be type RangedData") }
  lrr <- getBigMat(big.lrr,dir$big)
  baf <- getBigMat(big.baf,dir$big)
  incompletes <- which( (!sample %in% colnames(lrr)) | (!sample %in% colnames(baf)) )
  if(length(incompletes)>0) {
    sample <- sample[-incompletes]    
    if(length(sample)<1) { stop("no valid samples") } else {
      warning(length(incompletes),"sample(s) removed from list as not in both baf.file and lrr.file")
    }
  }
  pos.warn <- "pos must be a vector of length 2, or a matrix with 2 columns (start, end)"
  if(!is.null(dim(pos))) {
    if(dim(pos)[2]==2) {
      if(nrow(pos)!=length(chr)) {
        chr <- rep(chr[1],times=nrow(pos)); warning("chr, pos had different lengths, only chr[1] will be used")
      }
      st <- pos[,1]; en <- pos[,2]
    } else {
      stop(pos.warn)
    }
  } else {
    if(length(pos)!=2) { stop(pos.warn) } else {
      st <- pos[1]; en <- pos[2]; chr <- chr[1]
    }
  }
  sample.num.lrr <- match(sample,colnames(lrr))
  sample.num.baf <- match(sample,colnames(baf))
  rnlrr <- rownames(lrr)
  rnbaf <- rownames(baf)
  prv.chr <- 0; rnsnp <- character()
  for (cc in 1:length(chr)) {
    rst <- start(snp.info[chr[cc]]); ren <- end(snp.info[chr[cc]])
    force.chr.pos(Pos=c(st[cc],en[cc]),Chr=chr[cc],snp.info=snp.info)
    select <- (rst >= st[cc] & ren <= en[cc])
    if(length(which(select))<1) { warning("in set #",cc,"no snps in range"); next }
    if(chr[cc]!=prv.chr) { rnsnp <- rownames(snp.info[chr[cc]]) }
    snp.list <- rnsnp[select,]
    snp.list <- snp.list[(snp.list %in% rnlrr) & (snp.list %in% rnbaf)]
    n.snps <- length(snp.list)
    baf.y <- replicate(length(sample),sim.snp.baf(n.snps,type=type,noise=noise))
    lrr.y <- replicate(length(sample),sim.snp.lrr(n.snps,type=type,noise=noise))
    snp.poz.lrr <- match(snp.list,rnlrr)
    snp.poz.baf <- match(snp.list,rnbaf)
    if(add) {
      # add on top of current LRR (e.g, keeps organic noise, assuming already normal copy)
      lrr[snp.poz.lrr,sample.num.lrr] <- lrr[snp.poz.lrr,sample.num.lrr]+lrr.y
    } else {
      # completely replace current data with simulated data
      lrr[snp.poz.lrr,sample.num.lrr] <- lrr.y
    }
    baf[snp.poz.baf,sample.num.baf] <- baf.y
  }
  flush(lrr); flush(baf)
}



#spike.in.cnv(chr=1,pos=c(1.0e+07,1.2e+07),sample="5426924068_R04C02",dir=dir,add=F,type="normal",noise=0.5)

#cnv.plot(dir=dir,samples="5426924068_R04C02",LRR=T,BAF=T,Chr=1,Pos=c(1.0e+07,1.2e+07)+c(-10^6,10^6),tag.cnv=T,lrr.file="big.lrr",Cnv=c(1.0e+07,1.2e+07))




Xcn <- function(n,m,s,f,baf=rep(NA,n),cn=3,...) {
  #simulating a DUP with cn=3/4 copies
  # simulate n points with mean m, sd s
  # f - function type, baf= beta allele freq
  # retrieve p's for each genotype based on Baf
  cn <- length(f)-1
  baf.p <- baf.ps(baf,cn=cn,...) 
  X <- numeric(n)
  for (cc in 1:n)
  {
    r <- sample(cn+1,1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


Xh <- function(n,m,s,f,baf=rep(NA,n),p.fail=0,...) {
  #simulating a copy-1-DEL or LOHz/ROHm (no hets)
  # simulate n points with mean m, sd s, p=probability of erroneous het , 
  # f - function type, baf= beta allele freq 
  baf.p <- baf.ps(baf,...) # retrieve p's for each genotype based on Baf
  baf.p[,2] <- p.fail   # replace the middle (het) value with an error probability
  X <- numeric(n)
  for (cc in 1:n)
  {
    r <- sample(3,1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


X2 <- function(n,m,s,f,baf=rep(NA,n),...) {
  # simulate a BA.Freq for each snp, then place BAF-signal accordingly
  X <- numeric(n)
  baf.p <- baf.ps(baf,...) # retrieve p's for each genotype based on Baf
  for (cc in 1:n)
  {
    r <- sample(length(f),1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


Xn <- function(n,m,s) {
  # simulate simple normal distribution
  X <- rnorm(n,m,s)
  return(X)
}


baf.ps <- function(baf=NA,cn=2,mn=.13,stdev=.15) {
  # random BAF if NA, else use the value
  n <- length(baf); sel <- which(is.na(baf)) ; m <- length(sel)
  if(m>0) {
    baf[sel] <- pmin(.5,pmax(0.0001,rnorm(m,mn,stdev)))   # m,s derived from ichip
    TF <- as.logical(round(runif(m)))
    baf[sel][TF] <- 1-baf[sel][TF]
  } 
  baf <- as.numeric(baf)
  baf <- pmin(1,pmax(0,baf))
  if(cn==3) {
    baf.p <- cbind((baf^3), 3*((baf^2)*(1-baf)), 3*(baf*((1-baf)^2)), ((1-baf)^3))
  } else {
    if(cn==4) {
      baf.p <- cbind((baf^4), 4*((baf^3)*(1-baf)), 6*((baf^2)*((1-baf)^2)), 4*(baf*((1-baf)^3)),((1-baf)^4))
    } else { 
      #i.e, cn==2#
      baf.p <- cbind(baf*baf, 2*baf*(1-baf), (1-baf)*(1-baf))
    }
  }
  return(baf.p)
}


ceils <- function(m,s) { 1-abs(rnorm(1,m,s)) }
floors <- function(m,s) { abs(rnorm(1,m,s)) }
middles <- function(m,s) { rnorm(1,m,s) }


insert.sim <- function(x=c(1:length(lrr.vec)),lrr.vec,baf.vec,type="normal",noise=1){
  if(length(lrr.vec)!=length(baf.vec)) { stop("lrr and baf need to have the same length") }
  lrr.vec[x] <- sim.snp.lrr(n=length(x),type=type,noise=noise)
  baf.vec[x] <- sim.snp.baf(n=length(x),type=type,noise=noise)
  return(list(LRR=lrr.vec,BAF=baf.vec))
}


format.cnv <- function(type) {
  if(is.character(type)) { 
    type <- tolower(type) 
    if(type=="del1") { type <- "del" }
    if(type=="dup3") { type <- "dup" }
  } else {
    if(is.numeric(type)) {
      type <- max(0,min(type,4))
    } else {
      stop("invalid type")
    }
    type <- c("del0","del","normal","dup","dup4")[type+1]
  }
  return(type)
}


sim.snp.lrr <- function(n=1,type="normal",noise=1,mode="LRR",...) {
  ## allow lower/uppercase, del/dup without number defaults to copy 1/3,
  ## numeric values can also be entered as copy number 0,1,2,3,4
  type <- format.cnv(type=type)  
  gen.fn <- switch(type,normal=normal,loh=loh,del0=del0,del=del1,dup=dup3,dup4=dup4,...)
  out <- gen.fn(n=n,noise=noise,mode=mode)
  return(out)
}


sim.snp.baf <- function(n=1,type="normal",noise=1,...) {
  return(sim.snp.lrr(n=n,type=type,noise=noise,mode="BAF",...))
}


## states for simulation

del0 <- function(n=1,noise=1,mode="LRR"){
  if(mode=="LRR") {
    m1 <- c(-3.52)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(1.32)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    out <- runif(n,0,1)
  }
  return(out)
}


del1 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- c(-0.66)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.28)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.01,.5,0.01)
    s1 <- c(.025,.5,.025)*noise
    p.fail <- .0075*noise
    f1 <- c("floors","middles","ceils")
    out <- Xh(n,m1,s1,f1,baf=baf,p.fail=p.fail,...)
  }
  return(out)
}



normal <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...) {
  if(mode=="LRR") {
    m1 <- 0 + rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- .16*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.015,.5,0.015)
    s1 <- c(.0225,.025,.0225)*noise
    f1 <- c("floors","middles","ceils")
    out <- X2(n,m1,s1,f1,baf=baf,...)
  }
  return(out)
}


loh <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- 0 + rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- .21*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.01,.5,0.01)
    s1 <- c(.02,.5,.02)*noise
    p.fail <- .0075*noise
    f1 <- c("floors","middles","ceils")
    out <- Xh(n,m1,s1,f1,baf=baf,p.fail=p.fail,...)
  }
  return(out)
}


dup3 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- c(0.395)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.21)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0,.35,.65,0)
    s1 <- c(.0175,.025,.025,.0175)*noise
    #p1 <- c(1,3,3,1)/8
    f1 <- c("floors","middles","middles","ceils")
    out <- Xcn(n,m1,s1,f1,baf=baf,cn=3,...)
  }
  return(out)
}




dup4 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...) {
  if(mode=="LRR") {
    m1 <- c(.678)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.19)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0,.25,.5,.75,0)
    s1 <- c(.0175,.02,.025,.02,.0175)*noise
    #p1 <- c(1,4,6,4,1)/16
    f1 <- c("floors","middles","middles","middles","ceils")
    out <- Xcn(n,m1,s1,f1,baf=baf,cn=4,...)
  }
  return(out)
}


get.baf.lik <- function(dir="",n=100000,noise=2,recalc=F) {
  dir <- validate.dir.for(dir,"ano")
  store.path <- cat.path(dir$ano,"BAFlikelihood.RData")
  if(file.exists(store.path) & !recalc) {
    sim.bafs <- get(paste(load(store.path)))
    if(all(dim(sim.bafs)==c(101,5))) {
      return(sim.bafs)
    }
  }
  sim.bafs <- matrix(nrow=101,ncol=5); 
  rwz <- c(0:100/100)
  colnames(sim.bafs) <- paste(0:4); rownames(sim.bafs) <- paste(rwz)
  for (cc in 1:ncol(sim.bafs)) {
    sim.bafs[,cc] <- table(c(round(sim.snp.baf(n-101,type=cc-1,noise=noise),2),rwz))
  }
  likelihood.matrix <- log(sim.bafs / rowSums(sim.bafs))
  save(likelihood.matrix,file=store.path)
  return(likelihood.matrix)
}


#### ROH functions

extractROH <- function(dd,min.size=100,merge.gap.pc=.01,verbose=F) {
  missing <- which(dd==0)
  ld <- length(dd)
  #print(length(missing))
  for (cc in 1:length(missing)) {
    nxt <- missing[cc]; 
    rng <- dd[(max(1,nxt-3)):(min(nxt+3,ld))]
    md <- Mode(rng)
    dd[nxt] <- max(md,1)
  }
  dd[dd==2] <- 0
  dd[dd>0] <- 1
  yy1 <- rle(dd)
  poz <- c(1,1+cumsum(yy1[[1]])); poz <- poz[-length(poz)]
  starts <- poz[yy1[[1]]>min.size & yy1[[2]]]
  ends <- poz[which(yy1[[1]]>min.size & yy1[[2]])+1]-1
  pozzies <- cbind(starts,ends); Le <- length(ends)
  if(nrow(pozzies)>1) {
    ii <- which((starts[2:Le]-ends[1:(Le-1)])<=ceiling(merge.gap.pc*(ends[2:Le]-starts[1:(Le-1)])))
    if(verbose) { cat("ignoring",length(ii),"isolated heterozygous snps\n"); for (cc in 1:(length(ii))) { 
      cat(paste("merged",paste(pozzies[ii[cc],],collapse=","),"with",paste(pozzies[ii[cc]+1,],collapse=",")),"\n") 
    } }
    if(length(ii)>0) {
      if(any((ii+1)>nrow(pozzies))) { 
        #print(pozzies)
        stop(paste("ii",ii,"dim",dim(pozzies))) 
      }
      pozzies[ii,2] <- pozzies[ii+1,2]
      pozzies <- pozzies[-(ii+1),]
    }
  } else {
    if(length(pozzies)<1) { warning("empty range") ; return(pozzies) } else {
      if(starts[1] > ends[1]) { pozzies <- NA; warning("end > start") } else { dim(pozzies) <- c(1,2) }
    }
  }
  colnames(pozzies) <- c("start","end")
  return(pozzies)
}

convert.snp.indx.to.pos <- function(ROHlist,snpMat,snp.info) {
  must.use.package("genoset")
  snp.info <- toGenomeOrder(select.autosomes(snp.info),strict=T)
  valid.snps <- colnames(snpMat)[colnames(snpMat) %in% rownames(snp.info)]
  select <- which(rownames(snp.info) %in% valid.snps)
  if(length(valid.snps)<ncol(snpMat)) { warning("some snp names from snpMat were not found in snp.info - set to NA") }
  snp.info <- snp.info[select,]  
  st <- start(snp.info); en <- end(snp.info); ch <- chr(snp.info)
  #print(paste("st",st,"en",en,"chr",ch)[en>st])
  main.fun <- function(X,st,en,chr) {
    if(is.null(X) | length(which(!is.na(X)))<2) { warning("sample had no ROH regions exceeding min.size"); return(NULL) }
    naz <- which(is.na(X[,1]) | is.na(X[,2]))
    if(length(naz)>0) { X <- X[-naz,] }
    crossoverchr <- which(ch[X[,1]]!=ch[X[,2]]) # remove any instances where end of 1 chr joined to start of next
    if(length(crossoverchr)>0) {
      warning(paste("ROH region crossed from chr",ch[X[,1]][crossoverchr],":",st[X[,1]][crossoverchr],
                    " into chr",ch[X[,2]][crossoverchr],":",st[X[,2]][crossoverchr],
                    ", so was removed",{ if(length(crossoverchr)>1) {"\n"} else {""} },sep=""))
      X <- X[-crossoverchr,]
    }
    if(length(X)<1 | is.null(dim(X))) {
      if(length(X)==2) { 
        #print(st[X[1]]); print(en[X[2]]); print(ch[X[1]])
        newIR <- RangedData(ranges=IRanges(start=st[X[1]],end=en[X[2]]),space=ch[X[1]],nsnp=(1+(X[2]-X[1])))
        #newIR <- toGenomeOrder(newIR,strict=T) # this was commented out before??
      } else {
        newIR <- new("RangedData")
      }
    } else {
      if(any(en[X[,2]] < st[X[,1]])) { 
        stop("end position before start"); #print(paste("st",st[X[,1]],"en",en[X[,2]],"chr",ch[X[,1]])[en[X[,2]] > st[X[,1]]])
      }
      newIR <- RangedData(ranges=IRanges(start=st[X[,1]],end=en[X[,2]]),
                          space=ch[X[,1]],nsnp=(1+(X[,2]-X[,1])))
      #newIR <- toGenomeOrder(newIR,strict=T) # this was commented out before??
    } 
    return(newIR)
  }
  newList <- lapply(ROHlist,main.fun,st=st,en=en,ch=ch)
  names(newList) <- rownames(snpMat)  
  return(newList)
}


# extract runs of homozygosity for a SNP matrix (slow)
get.ROH.for.SnpMatrix <- function(snpMat,snp.info=NULL,snp.fail=NULL,dir=NULL,sample.list=NULL,...) {
  # either use 'dir' to automatically get snp.info and failures
  if(!is(snpMat)[1]=="SnpMatrix" & !is.raw(snpMat)) { warning("snpMat must be a SnpMatrix object") ; return(NULL) }
  if(!is.null(dir)) { dir <- validate.dir.for(dir,"ano") }
  if(is.null(snp.info) & is.null(dir)) { stop("at least one of snp.info or dir must be non-null to run this function") }
  doROHrow <- function(x,...) { dd <- as.numeric(x); out <- extractROH(dd,...); return(out) }
  if(!is.null(dir) & !is(snp.info)[1]=="RangedData") { snp.info <- read.snp.info(dir) }
  snp.info <- toGenomeOrder(select.autosomes(snp.info),strict=T)
  if(is.character(snp.fail)) { ii <- snp.fail } else { 
    if("QCfail" %in% colnames(snp.info)) {
      ii <- rownames(snp.info)[snp.info$QCfail!=0]
    } else {
      if(!is.null(dir)) { 
        ii <- get.all.snp.fails(dir) 
      } else { 
        warning("did not find any information on failing snps, all will be used") 
      }
    }
  }
  sf <- rownames(snp.info)
  sf <- sf[!sf %in% ii]
  dimmo <- dim(snpMat) ; snpRaw <- as.raw(snpMat); dim(snpRaw) <- dimmo
  rownames(snpRaw) <- rownames(snpMat); colnames(snpRaw) <- colnames(snpMat)
  if(is.character(sample.list)) {
    select.samps <- which(rownames(snpRaw) %in% sample.list)
    if(length(select.samps)<1) {
      stop("no samples from sample.list were found")
    }
  } else {
    select.samps <- 1:nrow(snpRaw)
  }
  select.snps <- narm(match(sf,colnames(snpRaw)))
  snpRaw <- snpRaw[select.samps,select.snps]
  if(length(select.samps)<1 | length(select.snps)<1 ) { warning("snp/sample selection gave empty SnpMatrix") ; return(NULL) }
  if(length(select.samps)==1 | length(select.snps)==1 ) {
    dim(snpRaw) <- c(length(select.samps),length(select.snps)) } # in case there is only 1 row or col 
  cat("processing",nrow(snpRaw),"samples for ROH ...")
  ROHlist <- vector("list",nrow(snpRaw)); uu <- system.time(ROHlist <- apply(snpRaw,1,doROHrow,...))
  save(ROHlist,file=cat.path(dir$ano,"ROH.list.pre.RData"))
  cat("took",round(uu[3]/60,1),"minutes\nnow converting to chr, pos format ... ")
  ROHpos <- convert.snp.indx.to.pos(ROHlist,snpRaw,snp.info) 
  ROHpos <- lapply(ROHpos,toGenomeOrder,strict=T)
  # names(ROHpos) <- rownames(snpRaw)[select.samps]
  cat("done\n")
  return(ROHpos)  
}


##
# extract plate information associated with the big.matrix
get.plate.info.for.big <- function(bigMat, dir) {
  dir <- validate.dir.for(dir,"big")
  plt <- get.plate.info(dir)
  plate.vec <- plt[[1]][,2][match(colnames(bigMat),plt[[1]][,1])]
  return(plate.vec)
}

# get plate means for a big matrix by using plate information
get.plate.snp.stats.for.big <- function(bigMat, dir, func=mean, n.cores=1,...) {
  dir <- validate.dir.for(dir,"big")
  plt <- get.plate.info.for.big(bigMat=bigMat,dir=dir)
  plate.func <- function(X,plt) { tapply(X,factor(plt),func,na.rm=T) }
  #  plate.sd <- function(X,plt) { tapply(X,factor(plt),sd,na.rm=T) }
  mnz.list <- bmcapply(bigMat,1,plate.func,plt=plt,n.cores=n.cores,dir=dir$big,...)
  #  sdz.list <- bmcapply(bigMat,F,plate.sd,plt=plt,n.cores=n.cores,dir=dir$big)
  return(mnz.list)
}

# for given genome ranges will try to find the closest snps to the end of the ranges
end.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(snp.info=snp.info,chr=chr,pos=pos,ranged=ranged,start=F,end=T,nearest=nearest))
}

# for given genome ranges will try to find the closest snps to the start and end of the ranges
range.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(snp.info=snp.info,chr=chr,pos=pos,ranged=ranged,start=T,end=T,nearest=nearest))
}

# for given genome ranges will try to find the closest snps to the start of the ranges
start.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,start=T,end=F,nearest=T) {
  # will preferably find an exact match but if nearest=T, will fall-back on nearest match
  must.use.package("genoset",T)
  nmz <- NULL
  if(!is(ranged)[1]=="RangedData") {
    if(!is.null(chr) & !is.null(pos)) {
      if(is.null(dim(pos))) { st <- pos[1]; en <- pos[2] } else {
        st <- pos[,1]; en <- pos[,2]
      }
      if(length(st)>length(chr)) { chr <- rep(chr[1],length(st)) } else { chr <- chr[1:length(st)] }
    } else {
      stop("if not using 'ranged', then chr and pos must be valid")
    }
  } else {
    st <- start(ranged); en <- end(ranged); chr <- chr(ranged)
    nmz <- rownames(ranged)
  }
  if(!is(snp.info)[1]=="RangedData") {
    stop("snp.info must be of type RangedData")
  } else {
    if(is.null(rownames(snp.info))) {  rownames(snp.info) <- paste(1:nrow(snp.info)) }
  }
  st.snps <- en.snps <- character(length(chr)) ; prch <- 0
  for (cc in 1:length(chr)) {
    if(chr[cc]!=prch) { 
      ref <- snp.info[paste(chr[cc])] ; st.ref <- start(ref); rnref <- rownames(ref)
      if(is.null(ref)) { stop(paste("snp.info did not contain chr",chr[cc])) }
    }
    #exact
    if(start) { 
      ind <- match(st[cc],st.ref)
      if(any(is.na(ind)) & nearest) {
        difs <- abs(st[cc]-st.ref)
        ind <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind)==0) { st.snps[cc] <- NA } else {
        st.snps[cc] <- rnref[ind]
      }
    }
    if(end){
      ind2 <- match(en[cc],st.ref)
      if(any(is.na(ind2)) & nearest) {
        difs <- abs(en[cc]-st.ref)
        ind2 <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind2)==0) { en.snps[cc] <- NA } else {
        en.snps[cc] <- rnref[ind2]
      }
    }
    prch <- chr[cc]
  }
  if(start & !end) {
    return(st.snps)
  }
  if(!start & end) {
    return(en.snps)
  }
  #otherwise looks like want both
  out <- cbind(st.snps,en.snps)
  if(!is.null(nmz)) { if(length(nmz)==nrow(out)) { rownames(out) <- nmz } }
  return(out)
}



## may want to add option to import from plink, BED, BIM, FAM !###
## annotBed <- gtfToBed(annotTrack)  ## install.packages("GeneticTools")
## bigmemoryExtras to use assaydata slot in eSet object
plumbCNV <- function(dir.base,dir.raw,snp.support="snpdata.map",gsf=gsf,delete.as.we.go=F,
                     dt.name="datatracker",run.mode="scratch",snp.run.mode="normal",
                     grps=NA,snp.fields=NULL,geno.file=NULL,big.lrr=NULL,big.baf=NULL,
                     aux.files.dir=NULL,plink.imp=F,
                     n.cores=1,q.cores=NA,grid.id="all.q",cluster.fn="q.cmd",low.ram=T,
                     start.at=0,pause.after=6,erase.previous=F,verbose=F,hide.penn.plink=T,
                     ucsc="hg18",
                     HD.mode=F,restore.mode=F,
                     callrate.samp.thr=.95,callrate.snp.thr=.95,hwe.thr=10^-8,
                     snp.grp.miss=F,grp.hwe.z.thr=4,grp.cr.thr=.001,het.lo=.1,het.hi=0.4,
                     nSD=3,mean.thr=c("LB","UB"),dlrs.thr=c(NA,"UB"),gc.thr=c("LB","UB"),
                     badPlateThresh=0.33,skip.chr.ab=F,lob=2,hib=2.5,pctile.bound=0.01,
                     cohort.pc.correct=F,num.pcs=9,
                     batch="plate",other.batch=list(),
                     lrr.report=T,chr.ab.report=T,plate.report=T,
                     pc.to.keep=.11,assoc=F,n.store=50,correct.sex=F,add.int=F,
                     comparison=T,comp.gc=F,comps="plate",use.penn.gc=F,
                     penn.path="/usr/local/bin/penncnv64/",hmm="hh550.hmm",relative=F,run.manual=F,print.cmds=F,
                     result.pref="cnvResults",out.format="Ranges",results="DT",print.summary.overlaps=F,
                     cnv.qc=T,rare.qc=T,plate.qc=T,pval=0.05,del.rate=0.4,dup.rate=0.18,thr.sd=3,plate.thr=3,
                     rmv.low.plates=F,min.sites=10,rare.olp=0.5,rare.pc=0.01,rmv.bad.reg=T,settings=NULL) 
{
  # run the entire plumbCNV pipeline with mostly default filenames.
  orig.dir <- getwd()
  section.vals <- c("Extraction from Genome Studio File","Import of LRR/BAF data",
                    "SNP-QC","SAMPLE-QC","PCA-CORRECTION","PENN CNV-CALLING","CNV-QC")
  start.at <- min(max(start.at,0),6); pause.after <- min(max(pause.after,0),6) # force legal values
  if(start.at!=0) { cat("plumbCNV pipeline attempting to (re)start at",section.vals[start.at+1],"\n\n") }
  if(pause.after!=6) { cat("NB: pipeline scheduled to pause after",section.vals[pause.after+1],"\n") }
  ### define directory structure
  proper.start <- ( (start.at==0 & run.mode=="scratch") | (start.at==1 & run.mode=="normal") )
  proper.start <- ( proper.start | (start.at==2 & run.mode=="big" & snp.run.mode!="skip") )
  proper.start <- ( proper.start | (start.at==3 & run.mode=="big" & snp.run.mode=="skip") )
  got.DT <- F
  if(!proper.start) {
    dir <- make.dir(dir.raw=dir.raw,dir.base=dir.base)
    DT <- read.data.tracker(dir,fn=dt.name,warn.only=T)
    if(is.data.tracker(DT)) {
      if(max(start.at-1,0) %in% as.numeric(unlist(getSlot(DT,"proc.done")))) {
        dir.DT <- getSlot(DT,"dir")
        if(dir$base!=dir.DT$base) {
          warning("plumbCNV session appears to have been moved - will attempt to smooth transition")
          cat("This session of plumbCNV analysis seems to have been moved at some point. ",
                  "Moving directories isn't recommended and could cause unpredictable results. ",
                  "However, if all parts are preserved it can work, will attempt to sync ",
                  "required files with the currently specified director. Cancel if you don't ",
                  "want this change, will wait 30 seconds",sep="")
          DT2 <- change.dir.data.tracker(dir.DT,dir.DT$base,dir$base)
          write.data.tracker(DT2,dir)
          DT <- DT2
          dir <- getSlot(DT,"dir")
        } else { dir <- dir.DT }
#         if(!is.null(settings)) { 
#           old.settings <- getSlot(DT,"settings")
#           if(!is.null(old.settings)) {
#            # settings <- update.list.with.list (old.settings, settings)
#             both.settings <- c(names(old.settings),names(settings))
#             both.dups <- duplicated(both.settings)
#             all.settings <- c(old.settings[!names(old.settings) %in% both.settings[both.dups]],settings)
#           }
#           DT <- setSlot(DT,settings=all.settings) 
#         }
        got.DT <- T
      } else {
        cat("Datatracker object did not contain indication of completion of the prior step:",
            section.vals[start.at],"\n")
        start.at <- switch(run.mode,scratch=0,normal=1,big=c(2,3)[1+(snp.run.mode=="skip")])
        cat("Changed settings to start at:",section.vals[start.at])
      }
    } else {
      cat("No valid datatracker object was present\n")
    }
  }
  if(!got.DT |  proper.start) {
    dir <- make.dir(dir.raw=dir.raw,dir.base=dir.base)
  }
  if(proper.start) {
    #prv(dir)
    init.dirs.fn(dir,overwrite=erase.previous,silent=T,update.bash=T,info.dir=aux.files.dir)
    erase.previous <- F  # change erase so it can't be done again below
    if(check.readiness(dir=dir,mode=run.mode,snp.mode=snp.run.mode)) { cat("plumbCNV file check successful\n") }
    if(!is.na(q.cores)) { if(!check.linux.install("qsub")) { q.cores <- NA } }
    # if file.spec.txt is missing, check for some essential parameters.
    tst.fst <- get.file.specs(dir)
    if(is.null(tst.fst)) { 
      cat(" file.spec.txt, not in use, files and settings must be entered manually\n")
      if(is.null(snp.fields) & snp.run.mode=="normal") { cat(" snp.fields is null, not valid without file.spec.txt\n") }
      if(is.null(geno.file) & snp.run.mode=="normal") { cat(" geno.file is null, not valid without file.spec.txt\n") }
      if(is.null(big.lrr) & snp.run.mode=="big") { cat(" big.lrr is null, not valid in run.mode='big'\n") }
      if(is.null(big.baf) & snp.run.mode=="big") { cat(" big.baf is null, not valid in run.mode='big'\n") }
    }
  }
  load.all.libs() ## load all libraries now
  ## 0. EXTRACT The DATA from Genome Studio Long Format file(s)
  if(run.mode=="scratch" & start.at==0) {
    Header("0. DATA EXTRACTION","#")
    sfiles <- init.DATA.read(dir,doLRR=T,doBAF=T,plink.imp=plink.imp,n.cores=n.cores,hwe.thr=hwe.thr,
                             callrate.snp.thr=callrate.snp.thr,callrate.samp.thr=callrate.samp.thr,
                             snp.info.sup=snp.support,genome.stud.file=gsf,combine.files=F) 
    proc.done <- 0
  } else {
    f01 <- list.files(dir$col,pattern="LRR.dat");  f02 <- list.files(dir$baf.col,pattern="LRR.dat")
    f1 <- list.files(dir$col); f2 <- list.files(dir$baf.col)
    if((length(f01)>0) & (length(f01)==length(f02)) & (length(f01)<length(f1))) {
      # i.e, if user has conformed to *.LRR.dat and *.BAF.dat suffixes select only these
      f1 <- f01; f2 <- f02  
    }
    if(length(f1)>0 & length(f2)>0) {
      sfiles <- list(lrr.files=f1, baf.files=f2)
      proc.done <- NULL
    } else {
      if(start.at<2) {
        cat("Error: Text data files (matrix or long format) missing from:\n LRR:",dir$col,
            "\n BAF:",dir$baf.col,"\n")
        stop("Please place LRR/BAF datafiles in the appropriate directories")
      }
    }
  }
  # unless starting after CNV-QC, these files should always be deleted if present
  if(start.at<7) { 
    if(is.list(dir)) {
      initialise.excl.files(dir,reset.files=c("BadReg.txt","CNVQC.txt")) 
    } else { warning("dir (directory list) object didn't load") }
  }
  if(!got.DT |  proper.start) {
    #?    dir <- make.dir(dir.raw=dir.raw,dir.base=dir.base)
    #?    if(!proper.start) { erase.previous <- F } # cancel erase if this is not a fresh start
    init.dirs.fn(dir,overwrite=erase.previous,silent=T,update.bash=T,info.dir=aux.files.dir)  
    # start a new datatracker, despite one existing, given we've started at the start
    #,snps="snpNames.txt",samples="subIdsALL.txt")
    DT <- init.data.tracker(dir,grps=grps,other.batch=other.batch,snp.fields=snp.fields,
                            geno.file=geno.file,big.lrr=big.lrr,big.baf=big.baf,ucsc=ucsc,verbose=verbose,
                            settings=settings) 
  }
  ## 1. IMPORT The DATA
  if(run.mode!="big" & start.at<2) {
    Header("1. DATA IMPORT","#")
    DT <- setSlot(DT,shape.lrr=sfiles$lrr.files, shape.baf=sfiles$baf.files,
                  proc.done=proc.done,grps=grps,settings=settings)
    write.data.tracker(DT,fn=dt.name)
    if(pause.after==0) {  return(DT) }
    DT <- import.DATA.big(DT) # convert these shaped datafiles in to bigmatrices
    write.data.tracker(DT,fn=dt.name)
    if(delete.as.we.go) { remove.for.step(DT,1,n.pcs=num.pcs) }
    if(pause.after==1) {  return(DT) }
  } else {
    big.fls <- unlist(getSlot(DT,"big.baf"))
    if(!is.file(big.fls,dir$big,dir)) {
      cat("Error: run mode was 'big' or start.at>1 but BAF big matrix isn't in:\n",
          dir$big,"\nor else the data.tracker object contains the wrong locations:\n")
      print(big.fls)
      stop("Please make sure big matrices are in the appropriate location, else use a different run.mode")
    }
    big.fls <- unlist(getSlot(DT,"big.lrr"))
    if(!is.file(big.fls,dir$big,dir)) {
      if(start.at<4) {
        cat("Error: run mode was 'big' or start.at>1 but LRR matrices aren't in:\n",
            dir$big,"\nor else the data.tracker object contains the wrong locations:\n")
        print(big.fls)
        stop("Please make sure big matrices are in the appropriate location, else use a different run.mode")
      }
    }
  }
  ## option to do own sample QC, then just pass callrate.txt lists to right location, 
  # or use file with these pre-excluded
  ## 2. Run SNQ-QC (SNP/Sample Call rate, HWE)
  # NB : only adds to SNP exlusion files, don't apply any filters to bigmats until (3)
  if(snp.run.mode!="skip" & start.at<3) { 
    Header("2. SNP QUALITY CONTROL","#")
    DT <- run.SNP.qc(DT=DT,import.plink=(snp.run.mode=="plink"),HD.mode=HD.mode,
                     restore.mode=restore.mode,ucsc=ucsc,
                     callrate.samp.thr=callrate.samp.thr,callrate.snp.thr=callrate.snp.thr, min.snps=min.sites,
                     hwe.thr=hwe.thr,group.miss=snp.grp.miss,grp.hwe.z.thr=grp.hwe.z.thr, grp.cr.thr=grp.cr.thr,
                     het.lo=het.lo,het.hi=het.hi,autosomes.only=T,n.cores=n.cores, low.ram=low.ram)
    DT <- setSlot(DT,settings=settings,proc.done=2)
    write.data.tracker(DT,fn=dt.name)
    if(delete.as.we.go) { remove.for.step(DT,2,n.pcs=num.pcs) }
    if(pause.after==2) {  return(DT) }
  } else  {
    cat("NB: SNP-QC skipped (due to parameter selection in snp.run.mode)\n")
    snp.info <- read.snp.info(dir)
    snp.info <- snp.update.qc.fail(snp.info,dir=dir)
    write.snp.info(snp.info,dir)
  }
  
  # 3. Run Sample QC
  # first thing is that it applies the filters derived in (2)
  if(start.at<4) {
    Header("3. SAMPLE QUALITY CONTROL","#")
    DT <- run.SAMPLE.qc(DT=DT,init=T,verbose=verbose,ucsc=ucsc,restore.mode=restore.mode,
                        nSD=nSD,mean.thr=mean.thr,dlrs.thr=dlrs.thr,gc.thr=gc.thr,
                        badPlateThresh=badPlateThresh,
                        skip.chr.ab=skip.chr.ab,lob=lob,hib=hib,pctile.bound=pctile.bound,
                        cohort.pc.correct=cohort.pc.correct,pc.to.keep=pc.to.keep,
                        num.pcs=num.pcs,batch=batch,n.cores=n.cores,
                        lrr.report=lrr.report,chr.ab.report=chr.ab.report,plate.report=plate.report)
    DT <- setSlot(DT,settings=settings,proc.done=3)
    write.data.tracker(DT,fn=dt.name)
    if(delete.as.we.go) { remove.for.step(DT,3,n.pcs=num.pcs) }
    if(pause.after==3) {  return(DT) }
  }
  
  # 4. Run PCA and PC-Correction to clean batch effects
  if(start.at<5) {
    Header("4. PCA AND PC-BATCH CORRECTION","#")
    if(assoc & cohort.pc.correct) { assoc <- F; warning("set 'assoc' to false as cohort.pc.correction removes mean differences between cohorts") }
    DT <- run.PCA.correct(DT=DT,pc.to.keep=pc.to.keep,assoc=assoc,num.pcs=num.pcs,n.store=n.store,correct.sex=correct.sex,add.int=add.int,
                          comparison=comparison,comp.gc=comp.gc,comps=comps,ucsc=ucsc,n.cores=n.cores,restore.mode=restore.mode)
    DT <- setSlot(DT,settings=settings,proc.done=4)
    write.data.tracker(DT,fn=dt.name)
    if(delete.as.we.go) { remove.for.step(DT,4,n.pcs=num.pcs) }
    if(pause.after==4) {  return(DT) }
  }
  
  # 5. Run PennCNV to call CNVs 
  if(start.at<6) {
    Header("5. RUN CNV CALLING WITH PENN-CNV","#")
    DT <- run.PENN.cnv(DT=DT,num.pcs=num.pcs,n.cores=n.cores,low.ram=T,q.cores=q.cores,grid.id=grid.id,
                       restore.mode=restore.mode,relative=relative,run.manual=run.manual,use.penn.gc=use.penn.gc,
                       hide.penn.out=hide.penn.plink,penn.path=penn.path,print.cmds=print.cmds,ucsc=ucsc,hmm=hmm)
    DT <- setSlot(DT,settings=settings,proc.done=5)
    write.data.tracker(DT,fn=dt.name)
    if(delete.as.we.go) { remove.for.step(DT,5,n.pcs=num.pcs) }
    if(pause.after==5) {  return(DT) }
  }
  
  # 6. Do CNV-QC
  if(start.at<7) {
    Header("6. CNV QC AND SUMMARY","#")
    DT <- run.CNV.qc(DT=DT,num.pcs=num.pcs,penn.path=penn.path,
                     ucsc=ucsc,hide.plink.out=hide.penn.plink,verbose=verbose,
                     out.format=out.format,result.pref=result.pref,
                     cnv.qc=cnv.qc,rare.qc=rare.qc,plate.qc=plate.qc,pval=pval,restore.mode=restore.mode,
                     del.rate=del.rate,dup.rate=dup.rate,thr.sd=thr.sd,plate.thr=plate.thr,
                     min.sites=min.sites,rare.olp=rare.olp,rare.pc=rare.pc,rmv.low.plates=rmv.low.plates,
                     rmv.bad.reg=rmv.bad.reg)
    DT <- setSlot(DT,settings=settings,proc.done=6)
  }
  write.data.tracker(DT,fn=dt.name)
  if(delete.as.we.go) { remove.for.step(DT,6,n.pcs=num.pcs) }
  prt <- sct <- CNVR <- NULL # initialise for 'return'
  cnvResults <- getSlot(DT,"cnv.calls",ret.obj=T,n.pcs=num.pcs)
  if(is.list(cnvResults) & length(cnvResults)==1) { cnvResults <- cnvResults[[1]] } #unlist
  if(is(cnvResults[[1]])[1]=="RangedData") {
    cnvResults <- lapply(cnvResults,toGenomeOrder) # ensure CNV results are in genome order
  }
  # temporary code to fiddle the phenotypes
#  sample.info <- read.sample.info(make.dir("/chiswick/data/ncooper/immunochipRunTest/"))
#  rn <- (cnvResults[[4]]$id)
#  whc <- match(rn,rownames(sample.info))
#  grpz <- sample.info$grp[whc]
#  grpz[grpz==-1] <- 0
#  print(unique(grpz))
#  cnvResults[[4]]$phenotype <- grpz
  ##
  if(print.summary.overlaps | !(tolower(results) %in% c("dt","ranges"))) {
    big.summary <- do.CNV.all.overlaps.summary(cnvResults,dir,comps=c(1:5),dbs=1:3,len.lo=1, len.hi=5000000,
                                               min.sites=min.sites,n.cores=1) #n.cores) 1 is faster?
    n.phenos <- length(unique(cnvResults[[1]][["phenotype"]]))
    if(print.summary.overlaps) { 
      #print("here")
      sct <- summary.counts.table(big.summary) ;
      #print("there") 
      if(n.phenos>1) {
        ofn <- cat.path(dir$res,"phenosummary",suf=result.pref,ext="txt")
        sink(ofn)
        prt <- print(pheno.ratios.table(dir,sum.table=sct))
        sink()
        #print(pheno.ratios.table(dir,sum.table=sct))
        cat("wrote to file:",ofn,"\n")
      }
    }
  }
  ofn <- cat.path(dir$res,"finalsummary",suf=result.pref,ext="txt")
  sink(ofn)
  print.snp.sample.summary(dir) # to file copy
  sink()
  print.snp.sample.summary(dir) # on screen copy
  if("id" %in% colnames(cnvResults[[1]])) { cat(length(unique(cnvResults[[1]]$id)),"samples had at least one CNV\n") }
  
  # PLOT COMPARITIVE CNVRs
  fn <- cat.path(dir$cnv.qc,pref="rare",fn=c("DEL","DUP"),ext="cnv.summary")
  if(any(file.exists(fn)) & n.phenos>1) {
    cnvr.del <- plot.pheno.cnvs(fn[1],type="DEL",pref=result.pref,dir=dir)
    cnvr.dup <- plot.pheno.cnvs(fn[2],type="DUP",pref=result.pref,dir=dir)
    CNVR <- list(cnvr.del,cnvr.dup) 
  }
  
  ## need to decide what to do about multiple sets of pca, eigen, cnv.result...? still allow it????
  setwd(orig.dir) # no matter where we've ended up, reset to the original directory  
  cat("\nplumbCNV pipeline complete\n")
  cat(" resulting *.cnv files: allCNV, allDel, allDup, rareDEL, rareDUP contain, respectively:\n")
  cat(" all CNVs, all deletions, all duplications, rare deletions, rare duplications. \n")
  cat(" These files can be processed and analysed using plink, or using R packages such as:\n")
  cat(" cnvGSA (association tests) or genoset (using the saved RangedData object, with \n")
  cat(" plumbCNV plot LRR/BAF functions and bioconductor overlap functions, etc)\n")
  cat(" thankyou for using plumbCNV.\n\n")
  
  if(tolower(results)=="dt") {
    return(DT)
  }
  if(tolower(results)=="overlaps") {
    return(big.summary)
  }
  if(tolower(results)=="ranges") {
    return(cnvResults)
  }
  return(list(ranges=cnvResults,overlaps=big.summary,table=sct,ratios=prt,cnvr=CNVR,DT=DT))
}


## version where any number of settings can be entered as a list ##
plumbcnv <- function(settings=list(),...) {
  other.args <- list(...)
  if(!is.list(settings)) { stop("settings should be a list() of plumbCNV arguments") }
  #settings <- update.list.with.list (settings, other.args)
  #both.settings <- duplicated(c(names(settings),names(other.args)))
  #all.settings <- c(settings[!names(settings) %in% both.settings],other.args)
  both.settings <- c(names(settings),names(other.args))
  both.dups <- duplicated(both.settings)
  all.settings <- c(settings[!names(settings) %in% both.settings[both.dups]],other.args)
  all.settings <- c(all.settings,settings=list(all.settings)) # add the settings themselves as an arg
  return(do.call(plumbCNV,args=all.settings))
}



# load other files with functions
#library(NCmisc)
#library(reader)
#source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/bigHelpers.R")
#source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/TempFunctions.R")
