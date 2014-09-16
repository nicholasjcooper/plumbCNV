


### FUNCTION INDEX ###
# apply.snp.thresholds - now remove samples failing HWE + callrate 95 regardless of source also optionally remove those failing group criteria if columns present in snp.info
# call.rate.summary - generate summary of callrate performance for snps/samples
# colSummary - do a snpStats::col.summary() only on selected set of snps
# het.density.plots make - density plot of Heterozygosity; can enter a file location/data.frame or vecs of heterozygosity
# hwe.density.plots make - density plot of HWE; can enter a file location/data.frame or vecs of hwe/
# hwe.vs.callrate.plots - make plots of HWE against callrate; can enter a file location/data.frame or vecs of hwe/
# hz.vs.callrate.plots - make plots of HWE against callrate; can enter a file location/data.frame or vecs of hwe/ 
###########################

## comment out,as now in SnpMatrixList.R
if(F) {
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
    MAF.cond <- snp.info$maf < maf
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
  if(!is.null(dir)) { dev.off(); cat("wrote file",ofn,"\n") }
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
  if(!is.null(dir)) { dev.off(); cat("wrote file:",ofn,"\n") }
}



hwe.vs.callrate.plots <- function(data="snpqc.txt",callrate.snp.thr=.95,hwe.thr=10^-5,
                                  Z.hwe=NULL,call.rate=NULL,zoom=T,full=T,dir=NULL,
                                  fn="HWEvsCallrate.pdf",excl=F, incl=F,jpg=FALSE) {
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
  if(!is.null(dir)) { 
   ofn <- cat.path(dir$cr,fn,ext=if(jpg) {"jpg"} else { "pdf" }); scl <- 1; 
    if(jpg) { jpeg(ofn) } else { pdf(ofn) }
  } else { 
    if(zoom & full) { par(mfrow=c(1,2)); scl <- 0.6 } else { scl <- 1 }
  }
  if(full) {
    plot(Z.hwe,call.rate,pch=".",xlab="HWE Z-score",ylab="call rate",
         main="full range",ylim=c(1,0),bty="n")
    add.boundary.and.legend(callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,scl=scl)
  }
  if(zoom) {
    if(zoom & full & jpg) { dev.off(); ofn <- cat.path("",ofn,suf="_zoom",ext="jpg"); jpeg(ofn) }
    plot(Z.hwe,call.rate,pch=".",xlim=c(-10,10),ylim=c(1,callrate.snp.thr-.01),
         xlab="HWE Z-score",ylab="call rate",main="cutoff range",bty="n")
    add.boundary.and.legend(callrate.snp.thr=callrate.snp.thr,hwe.thr=hwe.thr,scl=scl)
  }
  if(!is.null(dir)) { dev.off(); cat("wrote file:",ofn,"\n") }
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
  if(!is.null(dir)) { dev.off(); cat("wrote file:",ofn,"\n") }
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


