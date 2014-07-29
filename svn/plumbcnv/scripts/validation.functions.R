




# plumbCNV critical functions
  


# reads in the plink cnvr summary file and can either extract 'by.cnv', ie., a row for each
# CNV but with a column specifying which CNVR it belongs to, or else, just a list of CNVRs, with
# counts and FETs (if FET==T) for each
extract.cnv.regions <- function(dir, type="DEL", by.cnv=F, enriched=T, genes=T, lwr=.25,upr=4, FET=F, prt=F) {
  #probably only work when 2 phenotypes are present
  dir <- validate.dir.for(dir,"cnv.qc")
  must.use.package("genoset",T)
  type <- toupper(type); if(type!="DUP") { type <- "DEL" }
  cnv.ov <- cat.path(dir$cnv.qc,pref="rare",fn=type,ext="cnv.overlap")
  ## read in cnv overlap file
  #print(cnv.ov)
  tt <- read.table(cnv.ov,header=T,stringsAsFactors=FALSE)
  # tt <- reader(cnv.ov)
  inter <- tt[tt[,2]=="CON",]
  un <- tt[tt[,2]=="UNION",] # CNVRs
  pp <- strsplit(un[,4],":",fixed=T)
  case <- as.numeric(sapply(pp,"[",1)); ctrl <- as.numeric(sapply(pp,"[",2))
  ii <- ((case/ctrl)[case!=0 & ctrl!=0]) #[1:100]
  if(prt) { cat("Ratio summary:\n"); print(summary(ii)) }
  case0 <- case; case0[case0==0] <- 0.5; ctrl0 <- ctrl; ctrl0[ctrl0==0] <- 0.5;
  ii0 <- (case0/ctrl0)
  main <- tt[tt[,2]!="UNION" & tt[,2]!="CON",]  # results for individuals
  #print(table(tt[[1]]))
  #print(dim(main)); print(dim(tt)); 
  if(FET) {
    fet <- FET(case,ctrl,dir)
    CNVRs <- RangedData(IRanges(start=as.numeric(un[,6]),end=as.numeric(un[,7]),
                                names=un[,1]),space=as.numeric(un[,5]),cases=case,ctrls=ctrl,sig=fet)
  } else {
    CNVRs <- RangedData(IRanges(start=as.numeric(un[,6]),end=as.numeric(un[,7]),
                                names=un[,1]),space=as.numeric(un[,5]),cases=case,ctrls=ctrl)
  }
  
  if(enriched) {
    many.ctrl <- un[which(ii0<=lwr),]
    many.case <- un[which(ii0>=upr),]
    if(genes) {
      UNAFFECTED.REG <- many.ctrl[,"POOL"]; indx2 <- match(UNAFFECTED.REG,rownames(CNVRs))
      AFFECTED.REG <- many.case[,"POOL"]; indx <- match(AFFECTED.REG,rownames(CNVRs))
      CNVRs[["genes"]] <- (find.overlaps(CNVRs,db="gene",vec.out=T,delim=";",quiet=TRUE))
      many.ctrl <- cbind(many.ctrl,substr(CNVRs$genes[indx2],1,40))[,-match(c("FID","TYPE","SCORE"),colnames(many.ctrl))]
      many.case <- cbind(many.case,substr(CNVRs$genes[indx],1,40))[,-match(c("FID","TYPE","SCORE"),colnames(many.case))]
    } 
    if(prt) {
      cat("\nEnriched in unaffected summary:\n"); print(many.ctrl); 
      cat("\nEnriched in affected:\n"); print(many.case) # examples where case/ctrl outnumber each other >=4:1
    }
  }
  if(!by.cnv) {
    return(CNVRs)
  } else {
    cnvs_w_CNVRs <- RangedData(IRanges(start=as.numeric(main[,6]),
                                       end=as.numeric(main[,7])),space=as.numeric(main[,5]),
                               id=main[,3],phenotype=main[,4],cnvr=main[,1])
    return(cnvs_w_CNVRs)
  }
}


# taking the CNVR results for dels and dups, prints a summary table of top results
toptables <- function(oo1,oo2,pv=0.05,prt=T) {
  sig.dels <- which(oo1$sig<pv)
  sig.dups <- which(oo2$sig<pv)
  tt1 <- tt2 <- NULL
  if(prt) { cat("\nCNV-regions frequency by phenotype analysis: using Fishers Exact Test, showing results with p<",pv,"\n",sep="") }
  if(length(sig.dels)<1) { 
    warning("no significant DELs at p<",pv) 
  } else {
    tt1 <- as.data.frame(oo1[sig.dels,,drop=FALSE])
    colnames(tt1) <- c("chr","start (Mb)","end (Mb)","length (Kb)",
                       "n","cases","ctrls","sig","genes")
    tt1[["genes"]] <- as.character(tt1[["genes"]])
    rownames(tt1) <- tt1[,5]
    tt1 <- tt1[,-5] #-which(duplicated(tt1$genes)),-c(5)]
    tt1[,4] <- round(tt1[,4]/1000,1)
    tt1[,c(2,3)] <- round(tt1[,c(2,3)]/10^6,2)
    tt1[,7] <- format(tt1[,7],digits=3)
    if(any(is.na(tt1$genes))) { tt1$genes[is.na(tt1$genes)] <- "Unknown" }
    if(any(tt1$genes=="NA")) { tt1$genes[tt1$genes=="NA"] <- "Intergenic" }
    if(any(tt1$genes=="")) { tt1$genes[tt1$genes==""] <- "Intergenic" }
    DELtt <- cbind(tt1[,1:7],substr(tt1[,8],1,16))
    colnames(DELtt)[8] <- "genes"
    if(prt) {
      cat("DELETIONS:\n")
      print(DELtt[order(DELtt[,6]-DELtt[,5]),])
    }
  }
  if(length(sig.dups)<1) { 
    warning("no significant DUPs at p<",pv) 
  } else {
    tt2 <- as.data.frame(oo2[sig.dups,,drop=FALSE])
    colnames(tt2) <- c("chr","start (Mb)","end (Mb)","length (Kb)",
                       "n","cases","ctrls","sig","genes")
    tt2[["genes"]] <- as.character(tt2[["genes"]])
    rownames(tt2) <- tt2[,5]
    tt2 <- tt2[,-5] #-which(duplicated(tt2$genes)),-c(5)]
    tt2[,4] <- round(tt2[,4]/1000,1)
    tt2[,c(2,3)] <- round(tt2[,c(2,3)]/10^6,2)
    tt2[,7] <- format(tt2[,7],digits=3)
    if(any(is.na(tt2$genes))) { tt2$genes[is.na(tt2$genes)] <- "Unknown" }
    if(any(tt2$genes=="NA")) { tt1$genes[tt2$genes=="NA"] <- "Intergenic" }
    if(any(tt2$genes=="")) { tt2$genes[tt2$genes==""] <- "Intergenic" }
    DUPtt <- cbind(tt2[,1:7],substr(tt2[,8],1,16))
    colnames(DUPtt)[8] <- "genes"
    if(prt) {
      cat("DUPLICATIONS:\n")
      print(DUPtt[order(DUPtt[,6]-DUPtt[,5]),])
    }
  }
  return(list(DEL=tt1,DUP=tt2))
}


# creates summary of CNVs with different percentages passing quality thresholds
# could use a clean up
full.cnvr.summary <- function(cnvs1,oo1,dir,thr=c(.5,.75,.9)) {
  sample.info <- read.sample.info(dir)
  p.c <- table(sample.info$phenotype,sample.info$QCfail)[,1]
  case.d <- p.c[2]; cont.d <- p.c[1]
  all.reg <- unique(cnvs1$cnvr)
  #tapply(cnvs1$score[cnvs1$cnvr=="S20"],cnvs1$phenotype[cnvs1$cnvr=="S20"])
  #tapply(cnvs1$score[cnvs1$cnvr=="S20"],cnvs1$phenotype[cnvs1$cnvr=="S20"])
  nreg <- length(all.reg)
  ch <- character(nreg); nm <- numeric(nreg); int <- integer(nreg)
  results <- data.frame(t_ct=int,t_cs=int,t_sig=nm,
                        qs_ct=nm,qs_cs=nm,
                        qc50_ct=int,qc50_cs=int,qc50_sig=nm,
                        qc75_ct=int,qc75_cs=int,qc75_sig=nm,
                        qc90_ct=int,qc90_cs=int,qc90_sig=nm,
                        r1=ch,r1_ct=int,r1_cs=int,r1_sig=nm,
                        r2=ch,r2_ct=int,r2_cs=int,r2_sig=nm,
                        r3=ch,r3_ct=int,r3_cs=int,r3_sig=nm,stringsAsFactors=FALSE)
  typ <- sapply(sapply(results,is),"[",1)
  
  for (cc in 1:nreg) {
    # total counts by pheno
    ph.tl <- tapply(cnvs1$score[cnvs1$cnvr==all.reg[cc]],cnvs1$phenotype[cnvs1$cnvr==all.reg[cc]],length)
    ph.tl <- full.cnts12(ph.tl)
    ph.tl <- c(ph.tl,FET(ph.tl[2],ph.tl[1],case.d=case.d,cont.d=cont.d))      
    # scores by pheno
    ph.mn <- tapply(cnvs1$score[cnvs1$cnvr==all.reg[cc]],cnvs1$phenotype[cnvs1$cnvr==all.reg[cc]],mean)
    ph.mn <- full.cnts12(ph.mn)
    # counts above thresh
    ph.thr <- vector("list",length(thr)); names(ph.thr) <- paste(round(thr*100,0))
    for (jj in 1:length(thr)) {
      ph.thr[[jj]] <- tapply(cnvs1$score[cnvs1$cnvr==all.reg[cc]],cnvs1$phenotype[cnvs1$cnvr==all.reg[cc]],function(X) { length(which(X>thr[jj])) })
      ph.thr[[jj]] <- full.cnts12(ph.thr[[jj]])
      ph.thr[[jj]] <- c(ph.thr[[jj]],FET(ph.thr[[jj]][2],ph.thr[[jj]][1],case.d=case.d,cont.d=cont.d))
      names(ph.thr[[jj]]) <- c("ctrl","case","sig")
    }
    pvec <- cnvs1$phenotype[cnvs1$cnvr==all.reg[cc]]
    wdz <- tapply(width(cnvs1)[cnvs1$cnvr==all.reg[cc]],pvec,c)
    wdz <- full.cnts12(wdz)
    wds1 <-  wdz[[1]]
    wds2 <-  wdz[[2]]
    sets1.2 <- cluster.wd.set(wds1,wds2)
    sws <- summary.wd.set(sets1.2,wds1,wds2)
    
    cp <- function(x) { paste(x,collapse="-") }
    rngs <- paste(sapply(sws$ranges,cp))
    ctrl <- (sapply(sws$counts,"[",1))
    case <- (sapply(sws$counts,"[",2))
    case[is.infinite(case)] <- 0; ctrl[is.infinite(ctrl)] <- 0
    sum.info <- cbind(rngs,paste(ctrl),paste(case),paste(format(FET(case,ctrl,case.d=case.d,cont.d=cont.d),digits=3)))
    results[cc,] <- c(ph.tl,ph.mn,unlist(ph.thr),as.vector(t(sum.info)))
    loop.tracker(cc,nreg)
  }
  for(cc in 1:length(typ)) {
    results[[cc]] <- as(results[[cc]],typ[cc])
  }
  rownames(results) <- paste(unique(cnvs1$cnvr))
  print(head(results))
  results <- results[match(rownames(oo1),rownames(results)),]
  return(results)
}

## internal ; support function for full.cnvr.summary
full.cnts12 <- function(X) {
  # ensure there are counts for phenos 1 & 2
  inz <- which(paste(1:2) %in% names(X))
  if(length(inz)<2) { 
    if(length(inz)>0) {
      if(all(inz==1)) { 
        X <- c(X,0) 
      } else {
        if(all(inz==2)) { X <- c(0,X) }
      }
    } else {
      X <- c(0,0)
    }
  }
  return(X)
}

## internal ; support function for full.cnvr.summary
summary.wd.set <- function(X,wds1,wds2) {
  ngrp <- length(X)
  nlvl <- length(X[[1]])
  wds1.2 <- list(wds1,wds2)
  len.sum <- ran.sum <- vector("list",nlvl)
  for (cc in 1:nlvl) {
    rng <- gr <- vector("list",ngrp)
    for(dd in 1:ngrp) {
      gr[[dd]] <- X[[dd]][[cc]]
      the.dat <- wds1.2[[dd]][X[[dd]][[cc]]]
      rng[[dd]] <- c("","")
      if(length(the.dat)>0) {
        if(any(!is.na(the.dat))) {
          rng[[dd]] <- range(the.dat,na.rm=T)
        }
      }
    }
    ran.sum[[cc]] <- range(unlist(rng))
    if(any(is.infinite(ran.sum[[cc]]))) { ran.sum[[cc]] <- c("","") }
    len.sum[[cc]] <- sapply(gr,length)
  }
  return(list(counts=len.sum,ranges=ran.sum))
}
  
## internal ; support function for full.cnvr.summary
cluster.wd.set <- function(wds1,wds2) {
  if(length(wds1)<length(wds2)) {
    sw <- T; wds3 <- wds1; wds1 <- wds2; wds2 <- wds3
  } else { sw <- F }
  # do we need to worry about log10 existing?
  wds1 <- round(log10(wds1),1); wds2 <- round(log10(wds2),1)
  if((length(which(wds1==(Mode(wds1)[1])))/length(wds1))>.4) {
    cat1 <- which(wds1==Mode(wds1)[1])  } else { cat1 <- 1:length(wds1) }
  u <- (mean(wds1[cat1]));# print(u)
  zz <- ((wds1)-u)/sd((wds1))
  x1 <- which(zz==0)
  y1 <- which(zz>0)
  z1 <- which(zz<0)
  zz <- ((wds2)-u)/sd((wds2))
  x2 <- which(zz==0)
  y2 <- which(zz>0)
  z2 <- which(zz<0)
  wdsS1 <- list(x1,y1,z1); wdsS2 <- list(x2,y2,z2)
  if(sw) { wdsS3 <- wdsS1; wdsS1 <- wdsS2; wdsS2 <- wdsS3 }
  return(list(controls=wdsS1,cases=wdsS2))
}

# x.y.for.range.dat <- function(ranges,range.dat,genome=T) {
#   snp.info <- read.snp.info(dir)
#   snp.info <- snp.info[snp.info$QCfail==0,]
#   require(genoset)
#   if(nrow(ranges)!=length(range.dat)) { warning("mismatch between data lengths"); return(NULL) }
#   if(is(ranges)[1]!="RangedData") { stop("ranges must be RangedData")}
#   if(any(!names(range.dat) %in% ranges$id)) { warning("some range.dat names weren't in ranges$id"); return(NULL) }
#   y.dat <- range.dat
#   
# }





## functions very useful to others 




## take a cnvrs object (oo) and make a manhatten plots
# yl is the ylim (log10) for the graph
# lab01 is whether to label regions significant <.01
# lab05 is whether to label regions significant <.05
# regions below bonferroni will always be labelled
my.manhat <- function(oo,yl=8,lab05=F,lab01=F) {
  pch.case <- c(as.numeric(oo$cases < oo$ctrls)+15); pch.case[pch.case==15] <- 17
  print(table(pch.case))
  bnf <- .05/nrow(oo)
  thrs1 <- -log10(bnf); thrs2 <- 2 ;thrs3 <- -log10(.05) #2
  l10ps <- -log10(oo$sig)
  l10ps[l10ps==-Inf] <- -324
  l10ps[l10ps==Inf] <- 324
  gene.splt <- strsplit(oo$genes,";",fixed=T)
  lns <- sapply(gene.splt,length)
  shrt <- sapply(gene.splt,"[",1)
  shrt[lns>1] <- paste(shrt[lns>1],(lns-1)[lns>1],sep="+")
  shrt[is.na(shrt)] <- "interg."  #"Chr6: p25.3" #"interg."
  shrt[shrt==""] <- "interg." #"Chr6:p25.3"
  plot(10^-6*genoPos(oo),l10ps,pch=".",ylim=c(0,yl),xlim=c(0,2800),
       cex=2,xlab="genome position (Mb)",
       ylab="-log10 pvalues for case/control difference (F.E.T)")
  abline(h=-log10(bnf),lty="dotted") # bonferonni line
  sel <- which(l10ps>thrs3 & l10ps<=thrs2)
  dbls <- which(duplicated(shrt[sel]))
  if(length(dbls)>0) { sel1 <- sel[dbls]; sel <- sel[-dbls] } else { sel1 <- c() }
  if(length(sel1)>0) {
    points(x=(10^-6*genoPos(oo))[sel1], y=l10ps[sel1],cex=.3,col="white",pch=19)
  }
  if(length(sel)>0) {
    points(x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=.5,col="green",pch=pch.case[sel],bg="red")
    if(lab05) { text(labels=shrt[sel],x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=0.5,pos=3) }
  }
  sel <- which(l10ps>thrs2 & l10ps<=thrs1)
  dbls <- which(duplicated(shrt[sel]))
  if(length(dbls)>0) { sel1 <- sel[dbls]; sel <- sel[-dbls] } else { sel1 <- c() }
  if(length(sel1)>0) {
    points(x=(10^-6*genoPos(oo))[sel1], y=l10ps[sel1],cex=.3,col="white",pch=19)
  }
  if(length(sel)>0) {
    points(x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=.6,col="orange",pch=pch.case[sel],bg="red")
    if(lab01) { text(labels=shrt[sel],x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=0.6,pos=3) }
  }
  sel <- which(l10ps>thrs1)
  dbls <- which(duplicated(shrt[sel]))
  if(length(dbls)>0) { sel1 <- sel[dbls]; sel <- sel[-dbls] } else { sel1 <- c() }
  if(length(sel1)>0) {
    points(x=(10^-6*genoPos(oo))[sel1], y=l10ps[sel1],cex=.3,col="white",pch=19)
  }
  if(length(sel)>0) {
    points(x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=.7,col="red",pch=pch.case[sel],bg="red")
    text(labels=shrt[sel],x=(10^-6*genoPos(oo))[sel], y=l10ps[sel],cex=0.7,pos=3)
  }
  legend("top",legend=c("p < bonferroni","p < .01",  "p < .05", "p > .05",
                        "cases > controls", "controls > cases",
                        paste("bonferroni threshold (N=",nrow(oo),")",sep="")),
         col=c("red","orange","green",rep("black",4)),
         pch=c(19,19,19,15,17,16,NA),pt.cex=c(0.7,0.6,0.5,0.3,0.6,0.6,NA), 
         lty=c(NA,NA,NA,NA,NA,NA,"dotted"),ncol=2,cex=0.9)
}




#oo2 <- extract.cnv.regions(dir,type="dup",by.cnv=F,lwr=0.25,upr=4,FET=T,prt=F)
#cnvs2 <- extract.cnv.regions(dir,type="dup",by.cnv=T,lwr=0.25,upr=4,FET=T)
# might use with trios.analysis()
# filter CNV counts using a threshold, e.g, quality score
filt.counts.cnvr <- function(cnvrs,cnvs,col="QS",thresh=0.9,excl.less.than=T) {
  #cnvs <- toGenomeOrder2(cnvs,strict=T)
  if(excl.less.than) {
    cnvs <- cnvs[cnvs[[col]]>=thresh,]
  } else {
    cnvs <- cnvs[cnvs[[col]]<=thresh,]
  } 
  # return(cnvs)
  ii <- tapply(cnvs[["phenotype"]],cnvs[["cnvr"]],table)
  ctrls <- sapply(ii,"[","1")
  cases <- sapply(ii,"[","2")
  #print(head(ii)); print(head(CNVRs))
  # leave values as they were if new value is NA
  ##CNVRs[["cases"]][!is.na(match(rownames(CNVRs),names(ii)))] <- cases[narm(match(rownames(CNVRs),names(ii)))]
  ##CNVRs[["ctrls"]][!is.na(match(rownames(CNVRs),names(ii)))] <- ctrls[narm(match(rownames(CNVRs),names(ii)))]
  # any CNVS with all excluded go to NA ==> 0
  cnvrs[["cases"]] <- cases[(match(rownames(cnvrs),names(ii)))]
  cnvrs[["ctrls"]] <- ctrls[(match(rownames(cnvrs),names(ii)))]
  cnvrs[["cases"]][is.na(cnvrs[["cases"]])] <- 0
  cnvrs[["ctrls"]][is.na(cnvrs[["ctrls"]])] <- 0
  return(list(cnvr=cnvrs,cnv=cnvs))
}




#thisCNVR <- cnvs1[which(cnvs1$cnvr %in% rownames(tts[[1]])[16]),]
#plot.a.CNVR(thisCNVR,dir)
# for a CNV, make a plot (horizontal red and blue lines) comparing case control freq and extent
plot.a.CNVR <- function(thisCNVR,dir,genes=T) {
  ## assuming phenotypes 1,2 for controls,T1D respectively
  ## ... further arguments for getting the data
  require(genoset)
  thisCNVR <- thisCNVR[rev(order(width(thisCNVR))),] # sort by CNV length
  dat <- get.dat.for.ranges(thisCNVR,dir,pcCor=T,snp.pos=T,genome=F)
  ccc <- chr(thisCNVR)
  if(length(unique(ccc))>1) { stop("CNVs must belong to a single chromosome") }
  ccc <- ccc[1]
  mult <- if(genes) { 1.4 } else { 1.2 }
  pheno.list <- unique(thisCNVR$phenotype);
  if(length(pheno.list)<4) { col.list <- c("blue","red","green") } else { col.list <- get.distinct.cols(22) }
  cl <- col.list[match(thisCNVR$phenotype,pheno.list)]
  xl <- (range(unlist(dat$x))/10^6)  # *c(1,1.02)
  for(cc in 1:nrow(thisCNVR)) {
    xx <- unlist(dat$x[[cc]])
    cnt <- rep(cc,length(xx))
    if(thisCNVR$phenotype[cc]==1) { cl <- "blue" } else { cl <- "red" }
    if(cc==1) {
      plot(xx/10^6,cnt,col=cl, xlim=xl,type="l", ylim=c(0,nrow(thisCNVR)*mult),
           ylab=paste("count"),xlab="position (mb)",main=paste("Chromosome",ccc))
    }
    lines(xx/10^6,cnt+.1,col=cl) 
  }
  legend("top",legend=c("T1D","Controls"),col=c("red","blue"),pch="-",bty="n",ncol=2,pt.cex=2)
  if(genes) {
    gn <- plot.gene.annot(chr=ccc[1], pos=xl*10^6, x.scl=10^6, y.ofs=cc*1.15, width=cc/5, txt=T, chr.pos.offset=0,
                          build="hg18", dir=dir, box.col="green", txt.col="black", join.col="red")
    return(gn)
  }
  return(NULL)
}

# for a CNV, make a plot (dots or lines) comparing case control signal intensity, BAF coherence
plot.lrr.CNVR <- function(thisCNVR,dir,genes=T,BAF=F,type="l",DEL=T,cex=1) {
  ## assuming phenotypes 1,2 for controls,T1D respectively
  ## ... further arguments for getting the data
  require(genoset)
  thisCNVR <- thisCNVR[rev(order(width(thisCNVR))),] # sort by CNV length
  dat <- get.dat.for.ranges(thisCNVR,dir,pcCor=T,snp.pos=T,genome=F,BAF=BAF)
  ccc <- chr(thisCNVR)
  if(length(unique(ccc))>1) { stop("CNVs must belong to a single chromosome") }
  ccc <- ccc[1]
  mult <- if(genes) { 1.4 } else { 1.2 }
  pheno.list <- unique(thisCNVR$phenotype);
  if(length(pheno.list)<4) { col.list <- c("blue","red","green") } else { col.list <- get.distinct.cols(22) }
  cl <- col.list[match(thisCNVR$phenotype,pheno.list)]
  xl <- (range(unlist(dat$x))/10^6)  # *c(1,1.02)
  yl <- (range(unlist(dat$y)))  # *c(1,1.02)
  if(BAF) { yl <- c(0,.4) } else { yl[which(yl==max(yl))] <- max(yl)+(0.7*diff(yl)) }
  if(type!="l") { type <- "p" }
  tiny <- diff(xl)/10000
  for(cc in 1:nrow(thisCNVR)) {
    xx <- unlist(dat$x[[cc]])
    yy <- unlist(dat$y[[cc]])
    if(BAF) { yy <- abs(yy-.5)^2 } # to jitter on border; yy[yy>.24] <- yy[yy>.24]+runif(length(yy[yy>.24]),-0.005,0.005)}
    cnt <- rep(cc,length(xx))
    if(thisCNVR$phenotype[cc]==1) { cl <- "blue" } else { cl <- "red" }
    if(cc==1) {
      plot(xx/10^6,yy,col=cl, xlim=xl, ylim=yl, type=type,
           ylab=(if(BAF) { "BAF (squared difference to 0.5)" } else {"LRR" }),
           xlab="position (mb)",main=paste("Chromosome",ccc),pch=".",cex=cex)
    }
    if(type=="l") {
      lines((xx/10^6)+(cc*tiny),yy,col=cl) 
    } else {
      points((xx/10^6)+(cc*tiny),yy,col=cl, pch=".",cex=cex)
    }
  }
  legend("top",legend=c("T1D","Controls"),col=c("red","blue"),pch="-",bty="n",ncol=2,pt.cex=2)
  if(genes) {
    if(BAF) { y.ofsx=max(yl)*0.8; wd<-max(yl)/5 } else { y.ofsx=.5*max(yl); wd<-diff(yl)/5 }
    gn <- plot.gene.annot(chr=ccc[1], pos=xl*10^6, x.scl=10^6, y.ofs=y.ofsx, width=wd, txt=T, chr.pos.offset=0,
                          build="hg18", dir=dir, box.col="green", txt.col="black", join.col="red")
    return(gn)
  }
  return(NULL)
}

## functions a bit useful to others 


#re-does analysis using data passing quality score thresholds
#a bit messy but may be important to add this?
process.quality.scores <- function(DT,suffix,dir,restore=T) {
  # dodgy hack of dodgy code to extract quality scores.. quite some unecessary stuff here
  # works its way up to doing the case-control fischer analysis after extracting Quality scores,
  # and does it for several thresholds.
  thr <- 0 # set this to 0, leave filters for a later stage
  res.fn <- getSlot(DT,"cnvresult")
  cnvResults <- x.result <- get(paste(load(cat.path(dir$res,res.fn))))  #cnvResultsPCAAssoc9.RData"))))
  exp.names <-  c("allCNV","allDel","allDup","rareDEL","rareDUP")
  if(length(cnvResults)==5) { if(!all(names(cnvResults) %in% exp.names)) { names(cnvResults) <- exp.names } }
  #x.result$rareDEL <- x.result$rareDEL[!x.result$rareDEL$roh,]
  ofn <- cat.path(dir$res,"qs.del.backup",suf=suffix,ext="RData")
  if(F & file.exists(ofn) & restore) { print(load(ofn)) } else {
    qs.sc <- get.quality.scores(x.result$rareDEL,dir)
    save(qs.sc,file=ofn)
  }
  qs.mat <- make.qs.table(qs.sc)
  ofn <- cat.path(dir$res,"qs.del.results",suf=suffix,ext="txt")
  write.table(qs.mat,file=ofn,quote=F);   cat("wrote:",ofn,"\n")
  scrs <- qs.mat[[2]]; scrs[is.na(scrs)] <- qs.mat[[3]][is.na(scrs)]
  cat(length(scrs[scrs>=.95]),"/",length(scrs)," DELs pass .95 threshold\n",sep="")
  print(summary(scrs))
  # save QS to file
  #print(Dim(cnvResults))
  #cnvResults$rareDEL <- remove.duplicated.id.ranges(x.result$rareDEL); 
  cnvResults$rareDEL <- x.result$rareDEL;
  #print(Dim(cnvResults))
  cnvResults$rareDEL[["score"]] <- scrs;
  save(cnvResults,file=res.fn)
  x.result$rareDEL <- x.result$rareDEL[scrs>=thr,]
  ofn <- cat.path(dir$res,"qs.dup.backup",suf=suffix,ext="RData")
  if(F & file.exists(ofn) & restore) { print(load(ofn)) } else {
    qs.sc2 <- get.quality.scores(x.result$rareDUP,dir)
    save(qs.sc2,file=ofn)
  }
  qs.mat2 <- make.qs.table(qs.sc2)
  ofn <- cat.path(dir$res,"qs.dup.results",suf=suffix,ext="txt")
  write.table(qs.mat2,file=ofn,quote=F);   cat("wrote:",ofn,"\n")
  scrs2 <- qs.mat2[[2]]; scrs2[is.na(scrs2)] <- qs.mat2[[3]][is.na(scrs2)]
  cat(length(scrs2[scrs2>=.75]),"/",length(scrs2)," DUPs pass .75 threshold\n",sep="")
  print(summary(scrs2))
  # save QS to file
  #cnvResults$rareDUP <- remove.duplicated.id.ranges(x.result$rareDUP); 
  cnvResults$rareDUP <- x.result$rareDUP;
  cnvResults$rareDUP[["score"]] <- scrs2;
  save(cnvResults,file=res.fn)
  cat("wrote file results with quality scores inserted to:",res.fn,"\n")
  x.result$rareDUP <- x.result$rareDUP[scrs2>=thr,] # what is the different to cnvResults???
  
  ## add QS to cnvrs
  oo1 <- extract.cnv.regions(dir,type="del",by.cnv=F,lwr=0.25,upr=4,FET=T,prt=F)
  cnvs1 <- extract.cnv.regions(dir,type="del",by.cnv=T,lwr=0.25,upr=4,FET=T)
  oo2 <- extract.cnv.regions(dir,type="dup",by.cnv=F,lwr=0.25,upr=4,FET=T,prt=F)
  cnvs2 <- extract.cnv.regions(dir,type="dup",by.cnv=T,lwr=0.25,upr=4,FET=T)
  qs1 <- qs.mat; qs2 <- qs.mat2
  #print(load(getSlot(read.data.tracker(dir),"cnvresults")))
  cnvResults[[4]][["score"]] <- qs1[,1]
  cnvResults[[5]][["score"]] <- qs2[,1]
  cnvs1 <- add.scores.to.cnvrs(cnvs1,cnvResults[[4]],"score")
  cnvs2 <- add.scores.to.cnvrs(cnvs2,cnvResults[[5]],"score")
  results1 <- full.cnvr.summary(cnvs1,oo1,dir,thr=c(.5,.75,.9))
  results2 <- full.cnvr.summary(cnvs2,oo2,dir,thr=c(.5,.75,.9))
  return(list(DEL=results1,DUP=results2,cnvs1=cnvs1,cnvs2=cnvs2))
}


add.scores.to.cnvrs <- function(cnvs.r,cnvs.d,colname="score") {
  # for a dataset of individual samples with CNVRs,
  # add scores from the cnvResults to this, e.g, cnvs1, cnvResults[[4]]
  cnvs.r[[colname]] <- rep(NA,nrow(cnvs.r))
  all.sids <- unique(cnvs.r$id)
  for (dd in 1:length(all.sids)) {
    sset1 <- which(cnvs.r$id==all.sids[dd])
    sset2 <- which(cnvs.d$id==all.sids[dd])
    if(length(sset1)==length(sset2)) {
      cnvs.r[[colname]][sset1] <- cnvs.d[[colname]][sset2]
    } else {
      new <- (cnvs.r[sset1,])
      old <- (cnvs.d[sset2,])
      new.id <- paste(chr2(new),start(new),end(new),sep=".")
      old.id <- paste(chr2(old),start(old),end(old),sep=".")
      sset0 <- match(new.id,old.id)
      cnvs.r[[colname]][sset1] <- old[[colname]][sset0]
    }
  }
  return(cnvs.r)
}



# internal, quite specific to plot.all.samples.for.cnv()
is.cnvr.in.cnv.result.object <- function(cnvResult) {
  expect.names <- c("ranges","overlaps","table","ratios","cnvr","DT")
  if(!all(names(cnvResult) %in% expect.names)) { 
    if(length(cnvResult)==6) {
      names(cnvResult) <- expect.names
      cnvr <- which(expect.names=="cnvr")
    } else {
      all.names <- sapply(cnvResult,names); troo <- logical(length(all.names))
      for(cc in 1:length(all.names)) { 
        troo[cc] <- all(all.names[[cc]] %in% c("deletions","duplications")) & !is.null(all.names[[cc]])
      }
      if(length(troo)==1) { cnvr <- which(troo) } else { 
        stop("couldn't find good (single) candidate for cnvr sublist in the cnvResult object")
      }
    }
  } else { cnvr <- which(expect.names=="cnvr") }
  return(cnvr)
}

# take a CNV region from the top-table for DELs or DUPs and plot LRR,BAF all the individuals comprising it
# e,g
# 
# plot.all.samples.for.cnv(dir,reg="S1",dup=T,suffix="96")
plot.all.samples.for.cnv <- function(dir,reg="S1",dup=F,cnvResult=NULL,suffix="",LRR=T,BAF=T,PREPOSTPC=F,...) {
  DT <- read.data.tracker(dir)
  if(is.null(cnvResult)) { cnvResult <- reader(cat.path(dir$res,"fullresult",suf=suffix,ext="RData")) }
  if(length(cnvResult)<2) { if(is.character(cnvResult)) { if(cnvResult=="") { stop("cnv.result file not found in DT") }}}
  if(is.null(cnvResult)) { stop("object cnvResult must be valid, not NULL, to use this function")}
  cnvr <- is.cnvr.in.cnv.result.object(cnvResult)
  if(dup){ cnv.txt <- "DUP" } else { cnv.txt <- "DEL" }
  cnv.row <- cnvResult[[cnvr]][[(1+as.numeric(dup))]][reg,]
  cnv.lim <- c(chr2(cnv.row),start(cnv.row),end(cnv.row))
  cnv.ranges <- cnvResult[[1]][[(4+as.numeric(dup))]][cnv.lim[1]]
  cnv.ranges <- cnv.ranges[(start(cnv.ranges)>=cnv.lim[2] & end(cnv.ranges)<=cnv.lim[3]),]
  # do the plots - may be slow
  plot.all.ranges(cnv.ranges,DT=DT,file=cat.path("",cnv.txt,suf=reg,ext="pdf"),dir=dir,pc.flank=5,
                  LRR=LRR,BAF=BAF,PREPOSTPC=PREPOSTPC,...)   
  return(NULL)
}


## functions totally specific to me

## take a set of tables produced by the below commands and convert to independent counts/stats
# L1 <- 200000; fun <- logiti
# Rrr40 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap40")
# Rrr55 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap55")
# Rrr70 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap70")
# 
# L1 <- 200000; fun <- FET
# rr40 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap40")
# rr55 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap55")
# rr70 <- length.over.hap(del[del$score>=.95,],c(0,L1),T,FUN=fun,col="Hap70")
# 
# tab1 <- rbind(rr40,rr55,rr70)[c(1,3,5,2,4,6),]
# tab2 <- rbind(Rrr40,Rrr55,Rrr70)[c(1,3,5,2,4,6),]
#
# n1b <- 82; n2b <- 90; n1a <- 529; n2a <- 799
# independize(tab1,n1a,n2a,n1b,n2b,confs=c(0,55,70), labs=c("0-200 kb","> 200 kb"),FUN=FET)  
# n1b <- 6524; n2b <- 9238; n1a <- 6524; n2a <- 9238
# independize(tab1,n1a,n2a,n1b,n2b,confs=c(0,55,70), labs=c("0-200 kb","> 200 kb"),FUN=FET) 

independize <- function(tab,n1a,n2a,n1b,n2b,confs=c(0,55,70), labs=c("0-200 kb","> 200 kb"),FUN=FET,digits=6) { 
  new.tab <- as.data.frame(tab,stringsAsFactors=FALSE)
  for (cc in 1:ncol(new.tab)) {
    if(is(new.tab[[cc]])[1]=="factor") { new.tab[[cc]] <- paste(new.tab[[cc]]) }
  }
  lab.col <- 1
  hf <- nrow(tab)/2
  cas.con.cols <- c(2,3) # which columns have the case and control counts respectively
  change.rows <- (1:nrow(tab))[-c(nrow(tab),(hf))] # rows to change to difs
  if(tail(confs,1)!=100) { confs <- c(confs,100) }
  HapLabs <- character(hf)
  for (cc in 1:(length(confs)-1)) {
    HapLabs[cc] <- paste(confs[cc],"%-",confs[cc+1],"%",sep="")
  }
  for (ee in change.rows) {
    for (dd in cas.con.cols) {
      new.tab[ee,dd] <- tab[ee,dd]-tab[ee+1,dd] 
    }
  }
  new.tab[[lab.col]] <- rep(labs,each=hf)
  Cs <- cas.con.cols[1]; Ct <- cas.con.cols[2]
  n1 <- c(rep(n1a,hf),rep(n1b,hf))
  n2 <- c(rep(n2a,hf),rep(n2b,hf))
  for (cc in 1:nrow(tab)) {
    cS <- new.tab[cc,Cs]; cT <- new.tab[cc,Ct] # case and control counts
    #prv(cS,cT)
    new.tab$CI[cc] <- FUN(cS,cT,cont.d=n2[cc],case.d=n1[cc],stat="conf.int")
    new.tab$OR[cc] <- round(FUN(cS,cT,cont.d=n2[cc],case.d=n1[cc],stat="estimate"),3)
    new.tab$p.value[cc] <- round(FUN(cS,cT,cont.d=n2[cc],case.d=n1[cc]),digits)
    new.tab$n1[cc] <- n1[cc];     new.tab$n2[cc] <- n2[cc]; 
  }
  new.tab[["HapConfidence"]] <- rep(HapLabs,2)
  return(new.tab)
}
                        
more.than.0 <- function(X) { length(narm(X)[narm(X)>0]) }


# eg
# del[["Hap40"]] <- hap.mean(del,JS,FUN=num.more.than.55,n=.4)
# r40 <- length.over.hap(del[del$score>=.95,],LL2,T,rel.hap=T,col="Hap40")
# r50 <- length.over.hap(del[del$score>=.95,],LL2,T,rel.hap=T,col="Hap50")
# r60 <- length.over.hap(del[del$score>=.95,],LL2,T,rel.hap=T,col="Hap60")
# r70 <- length.over.hap(del[del$score>=.95,],LL2,T,rel.hap=T,col="Hap70")
# length.over.hap(del[del$score>=.95,],LL2,F,rel.hap=T,col="Hap55")
# length.over.hap(del[del$score>=.95,],LL2,F,rel.hap=F,col="Hap55")
# length.over.hap(del[del$score>=.95,],LL2,T,rel.hap=F,col="Hap55")
length.over.hap <- function(X,LL=0,gt=T,rel.hap=T,col="Hap55",FUN=FET) {
  nn <- length(LL)
  ll <- cnt <- matrix(nrow=nn,ncol=2)
  LL <- c(LL,250*10^6)
  CI <- character(nn)
  rat <- p <- numeric(nn)
  X <- remove.duplicated.id.ranges(X)
  for (cc in 1:nn) {
    if(gt) {
#      del.sh <- X[width(X) >= LL[cc],]      
      del.sh <- X[width(X) >= LL[cc] & width(X) <= LL[cc+1],]
      #print(cc) ;print(nrow(del.sh))
    } else {
      del.sh <- X[width(X) <= LL[cc],]
    }
    if(rel.hap) {
      ll[cc,] <- tapply(del.sh[[col]],factor(del.sh[["phenotype"]]),length)
    } else {
      ll[cc,] <- c(9238,6524) # cases and controls passing QC
    }
    if(nrow(del.sh)>0){
      cnt[cc,] <- tapply(del.sh[[col]],factor(del.sh[["phenotype"]]),more.than.0)
    } else { cnt[cc,] <- c(0,0) }
    CI[cc] <- FUN(cnt[cc,2],cnt[cc,1],cont.d=ll[cc,1],case.d=ll[cc,2],stat="conf.int")
    rat[cc] <- FUN(cnt[cc,2],cnt[cc,1],cont.d=ll[cc,1],case.d=ll[cc,2],stat="estimate")
    p[cc] <- FUN(cnt[cc,2],cnt[cc,1],cont.d=ll[cc,1],case.d=ll[cc,2])
  }
  LL <- head(LL,nn)
  results <- data.frame(length=LL,cases=cnt[,2], controls=cnt[,1], n1=ll[,2], n2=ll[,1], 
                        OR=rat, CI=CI, p.value=p)
  return(results)
}

num.more.than.55 <- function(X,n=.55) {
  if(!is.numeric(X)) { X <- as.numeric(X) }
  out <- length(X[X>n])
  return(out)
}

## retrieve a mean/median/max haploinsufficiency probability vector
hap.mean <- function(X,ref,FUN=mean,...,gene.col="gene",build=36) {
  if(!is(X)[1] %in% c("GRanges","RangedData")) { stop("X needs to be RangedData") }
  if(!is.data.frame(ref))  { stop("ref should be a dataframe containing haploinsufficiency scores by GENE id")}
  X <- Gene.pos(ranges=X,bioC=T,build=build)
  gn.l <- X[[gene.col]]
  gen.list <- paste(unique(gn.l))
  gen.list <- gen.list[gen.list!=""]
  LL <- length(gen.list)
  ll <- nrow(X)
  hapMean <- rep(NA,ll)
  for (cc in 1:LL) {
    gnz <- gen.list[cc]  # X[cc,gene.col]
    if(is.character(gnz)) { 
      if(length(gnz)>0) {  ZZZ <- ref[strsplit(gnz,";",fixed=T)[[1]],1] } 
      #prv(ZZZ)
      ZZ <- narm(as.numeric(ZZZ))
      if(length(ZZ)>0) { hapMean[gn.l %in% gnz] <- FUN(ZZ,...) }
    }
    #loop.tracker(cc,LL)
  }
  return(hapMean)
}

## retrieve a mean/median/max haploinsufficiency probability vector
hap.gene.list <- function(X,ref,gene.col="gene",thr=.55) {
  if(!is(X)[1] %in% c("GRanges","RangedData")) { stop("X needs to be RangedData") }
  if(!is.data.frame(ref))  { stop("ref should be a dataframe containing haploinsufficiency scores by GENE id")}
  X <- remove.duplicated.id.ranges(X)
  X <- Gene.pos(ranges=X,bioC=T)
  ph <- X$phenotype
  gn.l <- X[[gene.col]]
  gen.list <- paste(unique(gn.l))
  gen.list <- gen.list[gen.list!=""]
  LL <- length(gen.list)
  ll <- nrow(X)
  outlist <- vector("list",ll)
  hapMean <- rep(NA,ll)
  for (cc in 1:LL) {
    gnz <- gen.list[cc]  # X[cc,gene.col]
    if(is.character(gnz)) { 
      if(length(gnz)>0) {  
        ZZZ <- ref[strsplit(gnz,";",fixed=T)[[1]],1]
        YYY <- rownames(ref[strsplit(gnz,";",fixed=T)[[1]],]) 
        XXX <- narm(YYY[ZZZ>=thr & !is.na(ZZZ)])
      } 
      #prv(ZZZ)
      if(length(XXX)>0) { 
        #cat(XXX)
        band <- Band(XXX[1])
        times <- length(which(gn.l %in% gnz))
        cnts <- table(c(ph[gn.l %in% gnz],1,2))-1        
        OR <- round((cnts[["2"]]/6500) / (cnts[["1"]]/9500),3)
        cat(paste("T1D: ",cnts[2],", Ctrls: ",cnts[1],", OR: ",OR, ", Locus: ",band," | Haplogenes: ",comma(XXX),"\n",sep=""))
        cat("\n")
        #print(times)
        for (dd in 1:times) {
          outlist[gn.l %in% gnz][[dd]] <- XXX 
        }
      }
    }
    #loop.tracker(cc,LL)
  }
  return(outlist)
}

## add up number of case and control fails in DUPs and DELs
extract.qs.qc <- function(XX,DEL.thr=.9,DUP.thr=.9) {
  tfails1=which(XX[[4]]$phenotype==1 & XX[[4]]$score<DEL.thr)
  cfails1=which(XX[[4]]$phenotype==0 & XX[[4]]$score<DEL.thr)
  tfails2=which(XX[[5]]$phenotype==1 & XX[[5]]$score<DUP.thr)
  cfails2=which(XX[[5]]$phenotype==0 & XX[[5]]$score<DUP.thr)
  all.sel <- 1:nrow(XX[[4]])
  n.case.fails1 <- length(unique(XX[[4]]$id)) - length(unique(XX[[4]][all.sel[-tfails1],]$id)) 
  n.ctrl.fails1 <- length(unique(XX[[4]]$id)) - length(unique(XX[[4]][all.sel[-cfails1],]$id)) 
  all.sel <- 1:nrow(XX[[5]])
  n.case.fails2 <- length(unique(XX[[5]]$id)) - length(unique(XX[[5]][all.sel[-tfails2],]$id)) 
  n.ctrl.fails2 <- length(unique(XX[[5]]$id)) - length(unique(XX[[5]][all.sel[-cfails2],]$id)) 
  flist <- list(case.fails1=length(tfails1),n.case.fails1=n.case.fails1,
                ctrl.fails1=length(cfails1),n.ctrl.fails1=n.ctrl.fails1,
                case.fails2=length(tfails2),n.case.fails2=n.case.fails2,
                ctrl.fails2=length(cfails2),n.ctrl.fails2=n.ctrl.fails2)
}





## look at overall quality scores for RUNs and see whether it corresponds to quality of individual calls?
# not sure where this is/was used
get.correspondence <- function(all.runs.d,plots=F) {
  CNVID <- paste(all.runs.d$id,"_",chr(all.runs.d),":",
                 round(start(all.runs.d)/10^6), "-",round(end(all.runs.d)/10^6),sep="")
  
  #tapply(all.runs.d$score,factor(CNVID),function(x) { round(sum(x,na.rm=T),2) })
  mnz <- tapply(all.runs.d$score,factor(all.runs.d$RUN),function(x) { round(mean(x,na.rm=T),3) })
  
  #tapply(all.runs.d$score,factor(CNVID),function(x) { round(sum(x,na.rm=T),2) })
  #tapply(all.runs.d$RUN,factor(CNVID),function(x) { mnz[match(x,mnz)] })
  
  run.mnz <- mnz[match(paste(all.runs.d$RUN),names(mnz))]
  
  rnz <- tapply(all.runs.d$RUN,factor(CNVID),c)
  
  # is quality of the run related to quality of ind. calls, consistency, etx
  #  x1 <- tapply(run.mnz,factor(CNVID),function(x) { round(sum(x,na.rm=T),2) })
  y1 <- tapply(run.mnz,factor(CNVID),function(x) { round(mean(x,na.rm=T),2) })
  z1 <- tapply(run.mnz,factor(CNVID),function(x) { length(x) })
  #  z2 <- tapply(all.runs.d$score,factor(CNVID),function(x) { round(mean(x,na.rm=T),2) })
  if(plots) {
    smoothScatter(z2,z1)
    smoothScatter(z2,y1)
    smoothScatter(z2,x1)
  }
  sc.per.cnv <- y1[match(CNVID,names(y1))]
  cor(sc.per.cnv,all.runs.d$score)
  ## function to convert these CNVids to chr, pos format, useful for graphing in order
  cnv.id.to.pos <- function(X) {
    ii <- strsplit(X,"_",fixed=T)
    poz <- sapply(ii,tail,1)
    #idz <- character(length(ii))
    #for (cc in 1:length(ii)) { idz[cc] <- ii[[cc]][-length(ii[[cc]])] }
    pp <- convert.textpos.to.data(poz)
    #p[["id"]] <- idz
    return(pp)
  }
  #CNVID2 <- sapply(CNVID,cnv.id.to.pos)
  CNVID2 <- cnv.id.to.pos(CNVID)
  length.tab <- cbind(cnv.id.to.pos(names(z1)),z1)
  length.tab <- length.tab[order(length.tab[,2]),]
  length.tab <- length.tab[order(length.tab[,1]),]
  rownames(length.tab) <- NULL
  range.tab <- data.frame.to.ranged(as.data.frame(length.tab),stringsAsFactors=FALSE)
  nr <- length(z1); tc <- table(all.runs.d$RUN); nc <- length(tc)
  mm <- matrix(nrow=nr,ncol=nc); colnames(mm) <- names(tc); rownames(mm) <- names(rnz)
  colvec <- 1:ncol(mm)
  for(cc in 1:nrow(mm)) {
    mm[cc,] <- as.numeric(colvec %in% rnz[[cc]])
  }
  ### now cluster by columns!!
  pc.sim <- function(mm,rr,cc,zz) { length(which(as.logical(mm[,rr]) & as.logical(mm[,cc])))/(zz[cc]) }
  zz <- colSums(mm); res <- matrix(nrow=length(zz),ncol=length(zz))
  colnames(res) <- rownames(res) <-  names(tc)
  ## get similarity matrix between CNVs in all runs
  for (ii in 1:length(zz)) {
    for (jj in 1:length(zz)) {
      res[ii,jj] <- pc.sim(mm,ii,jj,zz)
    }
  }
  return(list(res=res,mnz=mnz))
}


## adjust sample fail counts for known inconsistency at time of running
update.samplefail.counts <- function(RUNS,cc) { 
  samp <- RUNS[[cc+1]][[2]]$SAMP
  denom.adj1 <- abs((RUNS[[cc]][[1]][["cases"]][[1]])[5,]-(RUNS[[cc+1]][[1]][["cases"]][[1]])[5,])
  denom.adj2 <- abs((RUNS[[cc]][[1]][["controls"]][[1]])[5,]-(RUNS[[cc+1]][[1]][["controls"]][[1]])[5,])
  denom.adj <- denom.adj1[1] + denom.adj2[1]
  if("samples passed pre-CNV QC" %in% names(samp)) {
    if(!"samples failed CNV-QC" %in% names(samp)) {
      # was here
      RUNS[[cc+1]][[2]]$SAMP[["samples failed CNV-QC"]] <- denom.adj
    } 
    all.pass <- samp[["samples passed pre-CNV QC"]]-RUNS[[cc+1]][[2]]$SAMP[["samples failed CNV-QC"]]
  } else {
    if("samples passed QC" %in% names(samp)) {
      all.pass <- RUNS[[cc+1]][[2]]$SAMP[["samples passed QC"]]
    } else {
      warning("insufficient SAMP data")
      all.pass <- NA
    }
  }
  if(!"samples passed QC" %in% names(samp)) {
    RUNS[[cc+1]][[2]]$SAMP[["samples passed QC"]] <- all.pass
  }
  if(!"samples passed pre-CNV QC" %in% names(samp)) {
    RUNS[[cc+1]][[2]]$SAMP[["samples passed pre-CNV QC"]] <- all.pass+denom.adj
  }
  ph <- RUNS[[cc+1]][[2]]$PH
  ph[1,] <- ph[1,]-c(denom.adj2[1],denom.adj1[1])
  ph[2,] <- ph[2,]+c(denom.adj2[1],denom.adj1[1])
  RUNS[[cc+1]][[2]]$PH <- ph
  ## fix for the CNV-QC off ones
  if(!"samples passed QC" %in% names(RUNS[[cc]][[2]]$SAMP)) {
    RUNS[[cc]][[2]]$SAMP[["samples passed QC"]] <- RUNS[[cc]][[2]]$SAMP[["samples passed pre-CNV QC"]]
  } 
  if(!"samples passed pre-CNV QC" %in% names(RUNS[[cc]][[2]]$SAMP)) {
    RUNS[[cc]][[2]]$SAMP[["samples passed pre-CNV QC"]] <- RUNS[[cc]][[2]]$SAMP[["samples passed QC"]]
  }
  RUNS[[cc]][[2]]$SAMP[["samples failed CNV-QC"]] <- 0
  
  return(RUNS)
}


## create the ratios taking into account pass/fail on QS criteria
update.divisors <- function(RUNS,cc,flist) {
  #'divisors'      #to run by rDELs,rDUPs
  if(cc %% 2==0) { 
    denom.adj1 <- abs((RUNS[[cc]][[1]][["cases"]][[1]])[5,]-(RUNS[[cc-1]][[1]][["cases"]][[1]])[5,])
    denom.adj2 <- abs((RUNS[[cc]][[1]][["controls"]][[1]])[5,]-(RUNS[[cc-1]][[1]][["controls"]][[1]])[5,])
    ## could add corrections to # passing all QC here
    ## + at least 1 rare CNV + at least 1 rare DUP
    ## + rDUPs
  } else { 
    denom.adj1 <-  denom.adj2 <- 0
  }
  PH <- RUNS[[cc]][[2]]$PH
  rnms <- c("total samples", "passing pre-CNV-QC", "passing ALL QC", 
            "at least 1 CNV", "at least 1 rare DEL", "at least 1 rare DUP",
            "pass QS rDELs Ss", "pass QS rDUP Ss","rDELs","rDUPs",
            "pass QS rDELs","pass QS rDUPs")
  divisors <- data.frame(controls=numeric(length(rnms)),cases=numeric(length(rnms)))
  ### if patch.on=T, then some of these counts should be modified (load if not in mem)
  if(patch.on & !exists("run.count.correction"))  { load(patch.file) }
  #run.samp.correction, run.count.correction
  if(cc %% 2==0 & patch.on) {
    ## get the numbers of controls, cases which need extra dups(add1), and the number of DUPs(add2)
    add1 <- c(run.samp.correction[[1]][paste(cc-1)],run.samp.correction[[2]][paste(cc-1)])
    add2 <- c(run.count.correction[[1]][paste(cc-1)],run.count.correction[[2]][paste(cc-1)])
  } else {
    add1 <- add2 <- 0
  }
  divisors[1,] <- c(sum(PH[1:2,1]),sum(PH[1:2,2]))
  divisors[2,] <- c(PH[1,1]+denom.adj2[1],PH[1,2]+denom.adj1[1])
  divisors[3,] <- c(PH[1,1],PH[1,2])+add1
  cnt <- RUNS[[cc]][[1]][["controls"]][[1]]
  cse <- RUNS[[cc]][[1]][["cases"]][[1]]
  if(!is.null(dim(cnt)) & !is.null(dim(cse))) {
    divisors[4,] <- c(cnt[5,1], cse[5,1])+add1
    divisors[5,] <- c(cnt[5,4], cse[5,4])
    divisors[6,] <- c(cnt[5,5], cse[5,5])+add1
    divisors[7,] <- c(flist$n.ctrl.fails1, flist$n.case.fails1)
    divisors[8,] <- c(flist$n.ctrl.fails2, flist$n.case.fails2)
    divisors[9,] <- c(cnt[4,4], cse[4,4])
    divisors[10,] <- c(cnt[4,5], cse[4,5])+add2
    divisors[11,] <- c(flist$ctrl.fails1, flist$case.fails1)
    divisors[12,] <- c(flist$ctrl.fails2, flist$case.fails2)
  }
  rownames(divisors) <- rnms
  divisors[["ratio"]] <- divisors[,1]/divisors[,2]
  divisors[["delORs"]] <- divisors[,3]/divisors[5,3]
  divisors[["dupORs"]] <- divisors[,3]/divisors[6,3]
  divisors[["QSdelORs"]] <- divisors[,3]/divisors[7,3]
  divisors[["QSdupORs"]] <- divisors[,3]/divisors[8,3]
  divisors[["delORc"]] <- divisors[,3]/divisors[9,3]
  divisors[["dupORc"]] <- divisors[,3]/divisors[10,3]
  divisors[["QSdelORc"]] <- divisors[,3]/divisors[11,3]
  divisors[["QSdupORc"]] <- divisors[,3]/divisors[12,3]
  RUNS[[cc]][["divisors"]] <- divisors
  return(RUNS)
}


# function to make the nice graphs of 54 runs, given a set of starting params
# runconfig contains the information for all the runs
## NB: it's a nasty function that uses globals
plot.dgv.result <- function(qs=T,cnv="DEL",db="ALL",an="SENS",stat="PC",nsnp=6) {
  if(an=="SENS") { xl <- "dgv_sensitivity" } else { xl <- "dgv_specificity" }
  if(qs) { yy <- mnz; yl <- "whole run mean quality score" }
  if(!qs) { yy <- narm(vec); yl <- "case:control CNV-ratio" }
  dat <- get.dgv.result(an=an,cnv=cnv,db=db,stat=stat,nsnp=nsnp)
  plot(dat,yy,pch=c(1,4,3)[runconfig[,5]+2], main=paste(cnv," [",db,", ",nsnp,"]",sep=""),
       col=c("red","blue","green")[runconfig[,2]+2],
       cex=c(1,2,3)[runconfig[,3]+2],bty="l",lwd=lwdz[runconfig[,4]+2],
       xlab=xl,ylab=yl) # sample-qc
}


### FUNCTION to extract data for sens/spec plotting , uses RUNS globally
get.dgv.result2 <- function(cnv="DEL",db="ALL",an="SENS",stat="PC",nsnp=6) {
  ## get data from the RUNS spec/sens data (which was filtered 
  # with q.score cutoffs but not for rare freq)
  dbs <- c(1,9,17)[which(c("CGHSEQ","BEAD","ALL") %in% toupper(db))]
  ns.skip <- (10-nsnp)/2 # rows to skip from starting point
  if(toupper(an)=="SENS") { ss <- 4 } else { ss <- 0 }
  if(toupper(cnv)=="DEL") { dd <- 0 } else { dd <- 1 }
  if(toupper(stat)=="PC") { pc <- 4 } else { pc <- 2 }
  dat <- unlist(sapply(RUNS,function(X) { if(length(X)>3) { X[["dgv"]][ss+ns.skip+dbs,dd+pc] } else { NULL } }))
  return(dat)
}

### FUNCTION to extract data for sens/spec plotting , uses res.list globally ?(better %s)
get.dgv.result <- function(cnv="DEL",db="ALL",an="SENS",stat="PC",nsnp=6) {
  ## get data from the res.list spec/sens data (which was filtered for rare freq,
  # and can be readily run with various filters, etc).
  dbs <- which(c("CGHSEQ","BEAD","ALL") %in% toupper(db))
  sel <- (toupper(paste(res.list$src))==toupper(db)) & (paste(res.list$db)==paste(nsnp))
  if(toupper(an)=="SENS") { ss <- 1 } else { ss <- 5 }
  if(toupper(cnv)=="DEL") { dd <- 0 } else { dd <- 1 }
  if(toupper(stat)=="PC") { pc <- 2 } else { pc <- 0 }
  dat <- unlist(tapply(res.list[sel,ss+dd+pc],factor(res.list[sel,]$run),c))
  ii <- pad.left(names(dat),"0")
  names(dat) <- ii; dat <- dat[order(names(dat))]
  #res.list[]
  return(dat)
}


### FUNCTIONS TO PASS EACH KIND OF OUTPUT FILE FROM THE RUNS ###
suck.metaboverlap.file <- function(fnm) {
  if(!file.exists(fnm)) { return(NULL) }
  fil <- readLines(fnm)
  all <- strsplit(fil," "); all <- lapply(all,function(X) { X[X!=""] })
  all <- unlist(all)
  all <- as.numeric(all)
  outlist <- narm(all)
  names(outlist) <- c("ichip DELs overlap with mchip",
                      "mchip DEL hits","mchip DEL unique hits", "ichip DUPs overlap with mchip",
                      "mchip DUP hits","mchip DUP unique hits",
                      "mchip DELs detectable on ichip","DELs found",
                      "mchip DELs likely detectable on ichip","DELs likely found",
                      "mchip DUPs detectable on ichip","DUPs found",
                      "mchip DUPS likely detectable on ichip","DUPS likely found")
  return(outlist)
}


suck.pheno.file <- function(fnp) {
  get.tab.with.anc <- function(fil,anc) {
    hdz <- tail(strsplit(fil[anc[1]+1]," ")[[1]],5)
    tab1 <- strsplit(fil[anc[1]+2:7]," ")
    tab1 <- lapply(tab1,function(X) { X[X!=""] })
    tabmain <- t(sapply(tab1,tail,5)); tabmain <- apply(tabmain,2,as.numeric)
    tabrn <- lapply(tab1,function(X) { paste(X[-c(length(X):(length(X)-4))],collapse=" ") })
    colnames(tabmain) <- hdz
    rownames(tabmain) <- tabrn
    return(tabmain)
  }
  if(!file.exists(fnp)) { return(NULL) }
  tab0 <- tab1 <- vector("list",3)
  fil <- readLines(fnp)
  anc <- grep("$raw.counts$`0`",fil,fixed=T)
  if(length(anc)>0) {
    for (cc in 1:length(anc)) {
      tab0[[cc]] <- get.tab.with.anc(fil,anc[cc])
    }
  } else { warning("no 0 counts found") }
  anc <- grep("$raw.counts$`1`",fil,fixed=T)
  if(length(anc)>0) {
    for (cc in 1:length(anc)) {
      tab1[[cc]] <- get.tab.with.anc(fil,anc[cc])
    }
  } else { warning("no 1 counts found") }
  names(tab0) <- names(tab1) <- c("Gene","Exon","DGV")
  return(list(controls=tab0,cases=tab1))
}



suck.final.file <- function(fnf) {
  get.list.with.anc <- function(fil,anc,sep,st,en) {
    tab1 <- strsplit(fil[anc[1]+st:en],sep)
    tabmain <- sapply(tab1,head,1); tabmain <- as.numeric(tabmain)
    tabrn <- paste(sep,sapply(tab1,tail,1))
    tabrn <- rmv.spc(tabrn); tabrn <- gsub("  "," ",tabrn)
    names(tabmain) <- tabrn
    return(tabmain)
  }
  get.list.with.anc2 <- function(fil,anc,sep,st,en,nval) {
    fil <- gsub("pheno ","pheno",fil)
    fil <- gsub("grp ","grp",fil)
    hdz <- tail(strsplit(fil[anc[1]+2]," ")[[1]],nval)
    tab1 <- strsplit(fil[anc[1]+st:en],sep)
    tab1 <- lapply(tab1,function(X) { X[X!=""] })
    tabmain <- t(sapply(tab1,tail,nval)); tabmain <- apply(tabmain,2,as.numeric)
    tabrn <- paste(sep,sapply(tab1,head,1),sep="")
    tabrn <- rmv.spc(tabrn); tabrn <- gsub("  "," ",tabrn)
    colnames(tabmain) <- hdz
    rownames(tabmain) <- tabrn
    return(tabmain)
  }
  if(!file.exists(fnf)) { return(NULL) }
  tabSNP <- vector("list",3)
  tabSAMP <- vector("list",4)
  tabCOH <- tabPH <- vector("list",2)
  fil <- readLines(fnf)
  anc <- grep("=== SNP SUMMARY ===",fil,fixed=T)
  tabSNP <- get.list.with.anc(fil,anc," SNPs",1,3)
  anc <- grep("=== SAMPLE SUMMARY ===",fil,fixed=T)
  tabSAMP <- get.list.with.anc(fil,anc," samples",1,4)
  anc <- grep("=== COHORT SUMMARY ===",fil,fixed=T)
  tabCOH <- get.list.with.anc2(fil,anc," ",3,4,3)
  anc <- grep("=== PHENOTYPE SUMMARY ===",fil,fixed=T)
  tabPH <- get.list.with.anc2(fil,anc," ",3,4,2)
  return(list(SNP=tabSNP,SAMP=narm(tabSAMP),COH=tabCOH,PH=tabPH)) 
}


suck.dgv.file <- function(fnd) {
  get.list.with.anc <- function(fil) {
    hdz <- strsplit(fil[1]," ")[[1]]; hdz <- hdz[hdz!=""]
    tab1 <- strsplit(fil[2:25]," ")
    tab1 <- lapply(tab1,function(X) { X[X!=""] })
    nval <- length(hdz)
    tabmain <- t(sapply(tab1,tail,nval)); tabmain <- as.data.frame(tabmain,stringsAsFactors=FALSE)
    tabmain[,2:5] <- apply(tabmain[,2:5],2,as.numeric)
    tabrn <- lapply(tab1,function(X) { paste(X[-c(1,length(X):(length(X)-nval+1))],collapse=" ") })
    colnames(tabmain) <- hdz
    tabmain[[1]] <- tabrn
    return(tabmain)
  }
  if(!file.exists(fnd)) { return(NULL) }
  fil <- readLines(fnd)
  tab <- get.list.with.anc(fil)
  return(tab) 
}


get.nxt.run <- function(nxt.r,rdir) {
  if(length(list.files(paste(rdir,"RUN",nxt.r,sep="")))<7) {
    return(paste("RUN",nxt.r," does not exist yet",sep=""))
  }
  fnp <- paste(rdir,"RUN",nxt.r,"/phenosummaryRUN",nxt.r,".txt",sep="")
  fnf <- paste(rdir,"RUN",nxt.r,"/finalsummaryRUN",nxt.r,".txt",sep="")
  fnm <- paste(rdir,"RUN",nxt.r,"/metaboverlapRUN",nxt.r,".txt",sep="")
  fnd <- paste(rdir,"RUN",nxt.r,"/dgv.validation.resultsRUN",nxt.r,".txt",sep="")
  
  stdz <- c(paste("qs.immuno.resultsRUN",nxt.r,".txt",sep=""), 
            paste("qs.dup.resultsRUN",nxt.r,".txt",sep=""), 
            paste("qs.del.resultsRUN",nxt.r,".txt",sep=""))
  stdz <- paste(rdir,"RUN",nxt.r,"/",stdz,sep="")
  
  ### extract each file type from next run ###
  next.pheno <- suck.pheno.file(fnp)
  next.final <- suck.final.file(fnf) 
  next.dgv <- suppressWarnings(suck.dgv.file(fnd))
  next.over <- suppressWarnings(suck.metaboverlap.file(fnm))
  regz  <- lapply(stdz,reader) ; names(regz) <- c("immuno","dup.qs","del.qs")
  outlist <- list(pheno=next.pheno,final=next.final,dgv=next.dgv,overlap=next.over)
  outlist <- c(outlist,regz)
  return(outlist)
}



## functions mostly specific to me



#internal
# gives the percent of a series failing a bonferroni test
percent.fail.bonf <- function(X,bonf=.05/length(X)) {
  if(max(X,na.rm=T)>1 | min(X,na.rm=T)<0) { stop("X must be p-values (and was <0 or >1)") }
  dd <- length(X)  
  nn <- length(X[X<bonf])
  return(nn/dd)
}


# take a DIL style file with IDs that are 'subject' or 'alternative' ids, and convert to sample ids
convert.to.sampid <- function(orig,alt,dir) {
  tt <- reader(orig)
  ss <- reader(alt)
  tt[["sampleid"]] <- ss[,1][match(tt[,1],ss[,2])]
  tt <- tt[!is.na(tt$sampleid),]
  if(length(which(duplicated(tt$sampleid)))>0) {
    tt <- tt[-which(duplicated(tt$sampleid)),]
  }
  rownames(tt) <- tt$sampleid
  tt <- tt[,-which(colnames(tt) %in% "sampleid")]
  sample.information <- tt
  write.csv(sample.information,file=cat.path(dir,basename(orig),suf="idconv"))
}


age.analyse <- function(dir,ph,corrected=F,select.col="phenotype",
                        select.val="no",cor.var="ageatbleed",cor.not.p=FALSE,n.cores=20) {
  # select.col is the column containing phenotype values
  # select.val is the value to analyse, e.g, 'ctrls', 't1d', 'cases', etc, depending on coding
  # cor.var is the variable to perform association with
  # corrected is whether to use the PC-corrected lrr matrix or the raw
  if(corrected) { mat.nm <- "big.pcc" } else { mat.nm <- "big.qc" }
  DT <- read.data.tracker(dir)
  big.pcc <- get.big.matrix(getSlot(DT,mat.nm)[[1]],dir$big) # raw
  #prv.big.matrix(big.pcc) #snp.info <- read.snp.info(dir)
  if(!cor.var %in% colnames(ph)) { stop("ph was missing column name",cor.var) }
  age <- ph[[cor.var]]
  names(age) <- ph[[1]]
  age.vec <- (age[(match(colnames(big.pcc),names(age)))])
  #print(head(ph))
  pheno <- (ph[[select.col]][(match(colnames(big.pcc),ph[[1]]))])
  #print(length(pheno))
  #ph[["phenotype"]] <- as.numeric(ph[[select.col]]!=select.val)
  select <- pheno==select.val
  select[is.na(select)] <- F
  bs1 <- big.select(big.pcc,select.cols=colnames(big.pcc)[select],dir=dir$big, pref="sel1")
  big.bs1 <- get.big.matrix(bs1,dir$big) # load descriptor as big.matrix
  if(is.null(rownames(ph)) & !anyDuplicated(ph[[1]])) { rownames(ph) <- ph[[1]] }
  PH1 <- ph[match(colnames(big.bs1),ph[[1]]),]
  #print(head(PH1));print(table(PH1$t1d));print(table(PH1$t1d,PH1$sample_type))
  #prv(pheno,select,age.vec,age,PH1,big.bs1,cor.var)
  if(cor.not.p) { a1 <- a2 <- FALSE } else { a1 <- a2 <- TRUE }
  qpa1 <- quick.pheno.assocs(big.bs1,PH1,cor.var,F.values=a1,p.values=a2,dir=dir$big,verbose=T,n.cores=n.cores)
  # run qpa with F.values=FALSE and p.values=FALSE to get correlations #
  if(cor.not.p) { 
    X <- qpa1; Y <- X[abs(X)>.1]
    qpa1 <- data.frame(r=as.numeric(X)); rownames(qpa1) <- names(X)
  } else { 
    X <- qpa1$p; Y <- X[abs(X)<(.05/nrow(qpa1))]
  }
  print(summary(X))
  cat(out.of(length(Y),length(X)),"\n")
  return(qpa1)
}

analyse.age.effect <- function(QPA,stat="percent.fail.bonf",scl=1000000,use.cor=FALSE,pref="") {
  if(scl==1000) { scl.lab <- "Kb" } else { scl.lab <- "base-pair"} 
  if(scl==1000000) { scl.lab <- "Mb" }
  if(!use.cor) {
    log10.1 <- (c(-1*log10(QPA$p)))
    log10.1[log10.1==-Inf] <- -324
    log10.1[log10.1==Inf] <- 324
    QPA[["log10"]] <- log10.1
  }
  QPA[["Pos"]] <- Pos(rownames(QPA))
  QPA[["Chr"]] <- Chr(rownames(QPA))
  QPA[["Mb"]] <- round(QPA[,"Pos"]/scl,0)
  QPA[["Chr.Mb"]] <- paste(QPA[["Chr"]],round(QPA[,"Pos"]/scl,0),sep=".")
  chrz <- gtools:::mixedsort(unique(QPA[,"Chr"]))
  mito <- which(chrz %in% c("MT","Mt","M","m","mt","XY","xy"))
  if(length(mito)>0) { chrz <- paste(chrz[-mito]) }
  tlm <- get.telomere.locs(bioC=T,kb=100)
  ctm <- get.centromere.locs(bioC=T)
  chrl <- get.chr.lens()/scl
  if(stat=="percent.fail.bonf") { 
    colm <- "p"  } else { 
      if(use.cor) { 
        colm <- "r" 
      } else { 
        colm <- "log10"
      } }
  ofn <- cat.path(getwd(),fn=stat,pref=pref,suf=colm,ext="pdf")
  pdf(ofn)
  statz <- vector("list",length(chrz))
  for (cc in 1:length(chrz)) {
    sel <- which(QPA[,"Chr"]==chrz[cc])
    if(stat=="percent.fail.bonf") {
      statz[[cc]] <- tapply(QPA[[colm]][sel],QPA$Mb[sel],function(X) { do.call(stat,args=list(X,bonf=.05/150000)) })
    } else {
      statz[[cc]] <- tapply(QPA[[colm]][sel],QPA$Mb[sel],function(X) { do.call(stat,args=list(X)) })
    }
    poz <- tapply(QPA$Mb[sel],QPA$Mb[sel],median)
    xl <- extend.pc(c(0,chrl[cc]),pc=.05)
    # main plot of analysis
    loess.scatter(poz,statz[[cc]],main=paste("Chromosome",chrz[cc]),type="l",
                  xlab=scl.lab,xlim=xl,ylab=paste(stat,colm,"value"))
    # add telomere and centromere to plot
    plot.ranges(rbind(tlm[chrz[cc]],ctm[chrz[cc]]),skip.plot.new=T,
                scl=scl.lab,ylim=c(0,.5),lty="dotted",full.vertical=TRUE)
  }
  dev.off()
  cat("wrote file:",ofn,"\n")
  return()            
}



### LINEAR CODE #####

if(F) {
  
  source("~/github/iChip/iFunctions.R")
  source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
  
  pheno.file <- "/home/ncooper/Documents/necessaryfilesICHIP/ichip.age.lookup.csv"
  
  setwd("/chiswick/data/ncooper/ImmunochipFamilies/")
  setwd("/chiswick/data/ncooper/immunochipRunTwo/")
  dir <- make.dir(getwd())
  DT <- read.data.tracker(dir)
  #ph <- reader("~/Documents/necessaryfilesICHIPFam/pheno.lookup.txt") # for famz
  ph <- reader(pheno.file)
  # ph <- shift.rownames(ph,T)
  #big.pcc <- get.big.matrix(getSlot(DT,"big.pcc"),dir$big) # corrected
  
  # define column names in pheno file
  pheno.var <- "t1d"
  age.var <- "age_ab_1" # "ageatbleed"
  
  #bs2 <- big.select(big.pcc,select.cols=colnames(big.pcc)[!ctrl],dir=dir$big, pref="sel2")
  #PH2 <- ph[match(colnames(get.big.matrix(bs2,dir$big)),rownames(ph)),]
  #qpa2 <- quick.pheno.assocs(get.big.matrix(bs2,dir$big),PH2,"ageatbleed",dir$big,n.cores=20)
  
  # correlations with sample type factor - check working correctly
  qpa1y <- NULL # age.analyse(dir,ph,corrected=F,select.col="t1d",select.val=2,cor.var="sample_type") # cases all same sample type
  qpa2z <- age.analyse(dir,ph,corrected=F,select.col="t1d",select.val=1,cor.var="sample_type")
  
  # this should work correctly now with new varnames , although phenotype should be t1d!
  # note that t1d are ALL cell-line and only have age for 1958bc, who are all same age
  
  qpa1a <- age.analyse(dir,ph,corrected=F,select.col=pheno.var,select.val=2,cor.var=age.var)
  qpa2a <- age.analyse(dir,ph,corrected=F,select.col=pheno.var,select.val=1,cor.var=age.var)
  qpa1b <- age.analyse(dir,ph,corrected=F,select.col="sample_type",select.val="CL",cor.var=age.var)
  qpa2b <- age.analyse(dir,ph,corrected=F,select.col="sample_type",select.val="wholeblood",cor.var=age.var)
  qpa3a <- age.analyse(dir,ph,corrected=F,select.col=pheno.var,select.val=2,cor.var=age.var,cor.not.p=T)
  qpa4a <- age.analyse(dir,ph,corrected=F,select.col=pheno.var,select.val=1,cor.var=age.var,cor.not.p=T)
  qpa3b <- age.analyse(dir,ph,corrected=F,select.col="sample_type",select.val="CL",cor.var=age.var,cor.not.p=T)
  qpa4b <- age.analyse(dir,ph,corrected=F,select.col="sample_type",select.val="wholeblood",cor.var=age.var,cor.not.p=T)
  
  qpa5a <- age.analyse(dir,ph,corrected=T,select.col=pheno.var,select.val=2,cor.var=age.var)
  qpa6a <- age.analyse(dir,ph,corrected=T,select.col=pheno.var,select.val=1,cor.var=age.var)
  qpa5b <- age.analyse(dir,ph,corrected=T,select.col="sample_type",select.val="CL",cor.var=age.var)
  qpa6b <- age.analyse(dir,ph,corrected=T,select.col="sample_type",select.val="wholeblood",cor.var=age.var)
  qpa7a <- age.analyse(dir,ph,corrected=T,select.col=pheno.var,select.val=2,cor.var=age.var,cor.not.p=T)
  qpa8a <- age.analyse(dir,ph,corrected=T,select.col=pheno.var,select.val=1,cor.var=age.var,cor.not.p=T)
  qpa7b <- age.analyse(dir,ph,corrected=T,select.col="sample_type",select.val="CL",cor.var=age.var,cor.not.p=T)
  qpa8b <- age.analyse(dir,ph,corrected=T,select.col="sample_type",select.val="wholeblood",cor.var=age.var,cor.not.p=T)
  
  save(qpa1a,qpa2a,qpa3a,qpa4a,qpa5a,qpa6a,qpa7a,qpa8a,qpa2z,
       qpa1b,qpa2b,qpa3b,qpa4b,qpa5b,qpa6b,qpa7b,qpa8b,file="QPA_New.RData")
  #save(qpa1,qpa2,qpa3,qpa4,qpa5,qpa6,qpa7,qpa8,file="QPA.RData")
  
  analyse.age.effect(qpa3a,stat="mean",scl=1000000,pref="t1dUncorrected",use.cor=T) #  cases
  analyse.age.effect(qpa5a,stat="mean",scl=1000000,pref="t1dCorrected",use.cor=T) #  cases
  
  
  qpa1 <- cbind(qpa1,qpa3) #  cases
  qpa2 <- cbind(qpa2,qpa4) #  ctrls
  
  analyse.age.effect(qpa1,stat="percent.fail.bonf",scl=1000000,pref="cases") #  cases
  analyse.age.effect(qpa2,stat="percent.fail.bonf",scl=1000000,pref="ctrls") #  ctrls
  
  analyse.age.effect(qpa1,stat="mean",scl=1000000,pref="cases",use.cor=T) #  cases
  analyse.age.effect(qpa2,stat="mean",scl=1000000,pref="ctrls",use.cor=T) #  ctrls
  analyse.age.effect(qpa1,stat="max",scl=1000000,pref="cases",use.cor=T) #  cases
  analyse.age.effect(qpa2,stat="max",scl=1000000,pref="ctrls",use.cor=T) #  ctrls
  analyse.age.effect(qpa1,stat="median",scl=1000000,pref="cases",use.cor=T) #  cases
  analyse.age.effect(qpa2,stat="median",scl=1000000,pref="ctrls",use.cor=T) #  ctrls
  
  # tables showing cor directions
  for (cc in 15:3) { bnf <- 10^-(cc); with(qpa2,print(paste(length(which(p<bnf)), (length(which(p<bnf & r>0))),(length(which(p<bnf & r<0)))))) }
  for (cc in 15:3) { bnf <- 10^-(cc); with(qpa1,print(paste(length(which(p<bnf)), (length(which(p<bnf & r>0))),(length(which(p<bnf & r<0)))))) }
  
  
  # working to merge pheno files
  
  orig <- "~/Documents/necessaryfilesMCHIP/SubjectPlateRegionTypeSexEthnic.txt"
  alt <- "~/Documents/necessaryfilesICHIP/alt.id.lookup.txt"
  convert.to.sampid(orig,alt,dir="~/Documents/necessaryfilesICHIP/")
  ## now go into CSV file and add 'ID/sampleid', etc as first header instead of ""
  dd1 <- reader(file.choose())
  dd2 <- reader(file.choose())
  dd <- gtools:::smartbind(dd1,dd2)
  rownames(dd) <- dd[[1]]
  dd <- dd[,-1]
  write.csv(dd,file=cat.path("~/Documents/necessaryfilesICHIP/","ichip.age.lookup",ext="csv"))
  
  tt <- reader(file.choose())
  
  
}



### LINEAR CODE ######

if(F) {
  #cbind(oo[[1]],oo[[2]],substr(oo[[3]],1,10))[order(oo[[1]]/oo[[2]]),]
  
  oo2 <- extract.cnv.regions(dir,type="dup",by.cnv=F,lwr=0.25,upr=4,FET=T)
  oo1 <- extract.cnv.regions(dir,type="del",by.cnv=F,lwr=0.25,upr=4,FET=T)
  pdf(cat.path(dir$res,"qqDEL.pdf")); qqplot(x=runif(nrow(oo1)),y=oo1$sig,type="l",ylab="p.value",xlab="uniform distribution 0,1"); abline(a=0,b=1,lty="dotted") ;dev.off()
  pdf(cat.path(dir$res,"qqDUP.pdf")); qqplot(x=runif(nrow(oo2)),y=oo2$sig,type="l",ylab="p.value",xlab="uniform distribution 0,1"); abline(a=0,b=1,lty="dotted") ;dev.off()
  
  cnvs2 <- extract.cnv.regions(dir,type="dup",by.cnv=T,lwr=0.25,upr=4,FET=T)
  cnvs1 <- extract.cnv.regions(dir,type="del",by.cnv=T,lwr=0.25,upr=4,FET=T)
  
  #pdf("manhatDEL8.pdf") ; my.manhat(oo1,9) ; dev.off()
  #pdf("manhatDUP8.pdf") ; my.manhat(oo2,6) ; dev.off()
  
  pdf("manhatDEL8b.pdf") ; my.manhat(oo1,19) ; dev.off()
  pdf("manhatDUP8b.pdf") ; my.manhat(oo2,7.5) ; dev.off()
  
  pdf("manhatDEL20.pdf") ; my.manhat(oo1,9.5) ; dev.off()
  pdf("manhatDUP20.pdf") ; my.manhat(oo2,7.5) ; dev.off()
  
  #save(oo1,oo2,cnvs1,cnvs2,file="RUN8_CNVR.RData")
  save(oo1,oo2,cnvs1,cnvs2,file="RUN20_CNVR.RData")
  sample.info <- read.sample.info(dir)
  p.c <- table(sample.info$phenotype,sample.info$QCfail)[,1]
  case.d <- p.c[1]; cont.d <- p.c[2];
  #file.choose(); setwd("/Users/nick/Dropbox/FYR/goodmaterials/")
  print(load("RUN19_CNVR.RData"))
  tts <- toptables(oo1,oo2,.01)
  ## SET tt = 1/2 to make table of comparitive cont-cont counts for dels/dups
  tt <- 1
  res.tab <- matrix(nrow=nrow(tts[[tt]]),ncol=4); colnames(res.tab) <-
    c("T1D","Control-UVA","Control-Sanger")
  for (cc in 1:nrow(tts[[tt]])) {
    if (tt==1) {
      thisCNVR <- cnvs1[which(cnvs1$cnvr %in% rownames(tts[[tt]])[cc]),]
    } else {
      thisCNVR <- cnvs2[which(cnvs2$cnvr %in% rownames(tts[[tt]])[cc]),]
    }
    controls <- thisCNVR[thisCNVR$phenotype==1,]
    cases <- thisCNVR[thisCNVR$phenotype==2,]
    ctrl.grp <- 1+as.numeric(substr(controls$id,1,1)=="7")
    uva.c <- length(which(ctrl.grp==1))
    san.c <- length(which(ctrl.grp==2))
    print(paste(rownames(tts[[1]])[cc],tts[[tt]][[1]][cc]))
    res.tab[cc,] <-c(nrow(cases),uva.c,san.c,
                     (fisher.test(x=matrix(c(uva.c,5461,san.c,4537),nrow=2))$p.value))
  }
  # rownames(res.tab) <- rownames(tts[[tt]]); res.tab <- cbind((tts[[tt]][[1]]),res.tab)
  # print(res.tab)
  
  cl1 <- table(width(controls))
  cl2 <- table(width(cases))
  print(cl1); print(cl2)
  #}
  if(length(cl1)>1) {
    while((as.numeric(names(cl1)[1])/as.numeric(names(cl1))[length(cl1)])<.33) {
      cl1 <- cl1[-1]; cat(".")
    }
  }
  if(length(cl2)>1) {
    while((as.numeric(names(cl2)[1])/as.numeric(names(cl2)[length(cl2)]))<.33) {
      cl2 <- cl2[-1]; cat(".")
    }
  }
  print(cl1); print(cl2)
  
  print(fisher.test(x=matrix(c(sum(cl2),case.d,sum(cl1),cont.d),nrow=2))$p.value)
  # print(chisq.test(cl1,cl2))
  #}
  #plot.all.ranges(controls)
  nxt <- controls
  rng <- c(min(start(nxt)),max(end(nxt)))
  cnv.plot(dir=dir,samples=nxt$id,BAF=T,PREPOST=T,cnvPlotFileName="ctrlschr22b.pdf",
           Chr=chr(nxt)[1],Cnv=rng,Pos=rng+(c(-3,3)*diff(rng)),tag.cnvs=T)
  nxt <- cases
  rng <- c(min(start(nxt)),max(end(nxt)))
  cnv.plot(dir=dir,samples=nxt$id,BAF=T,PREPOST=T,cnvPlotFileName="caseschr22b.pdf",
           Chr=chr(nxt)[1],Cnv=rng,Pos=rng+(c(-1,1)*diff(rng)),tag.cnvs=T)
  
  DT <- read.data.tracker(dir)
  bigPCC <- getBigMat(getSlot(DT,"big.pcc"),dir$big)
  bigLRR <- getBigMat(getSlot(DT,"big.qc"),dir$big)
  snp.info <- read.snp.info(dir)
  sample.info <- read.sample.info(dir)
  # cnv <- big.extract.snp.ranges(range.snp(snp.info,nxt),samples=nxt$id,bigPCC)
  bigD <- bigLRR
  
  #range.snp(snp.info,nxt)
  st.en <- range.snp(snp.info,nxt)[19,]
  #"rs7284771","rs2283798" biggest 450k
  #"rs8140331","rs7284345" common 7k
  #"rs7284771","rs7284345" common big 300k
  st.en <- c("imm_6_369915","imm_6_372152")   #13 0, chr 6
  st.en <- c("imm_20_61838983","imm_20_61850194")   #13 0, chr 20
  
  st.en <- c("rs7284771","rs7284345")
  st.en <- c("rs8140331","rs7284345")
  st.en <- c("imm_5_150137553","imm_5_150295407")
  rr <- match(st.en[1],rownames(bigD)):match(st.en[2],rownames(bigD))
  
  #match("imm_1_8110710",rownames(bigD))
  #match("imm_1_8113791",rownames(bigD))
  ####rr <- 1448:1455
  #match("imm_1_8110710",rownames(toGenomeOrder(snp.info)))
  tids <- which(colnames(bigD) %in% rownames(sample.info[sample.info$phenotype==1,]))
  cids <- which(colnames(bigD) %in% rownames(sample.info[sample.info$phenotype==0,]))
  #rowMeans(bigD[rr,cids],na.rm=T)
  #rowMeans(bigD[rr,cids],na.rm=T)
  #summary(bigD[rr,cids])
  #summary(bigD[rr,cids])
  #plot(bigD[1448,cids]) # :1455
  #boxplot(t(bigD[1448:1455,cids])) # :1455
  jitz1 <- runif(length(cids),-.5,.5); jitz2 <- runif(length(tids),-.5,.5)
  bigD <- bigPCC# bigLRR
  dat1 <- (colMeans(bigD[rr,cids],na.rm=T)); 
  dat2 <- (colMeans(bigD[rr,tids],na.rm=T)); 
  dat1[dat1< -2] <- -2; dat2[dat2< -2] <- -2
  cx1 <- (20*abs(dat1)); cx1[cx1<0.75] <- .75
  cx2 <- (20*abs(dat2)); cx2[cx2<0.75] <- .75
  yl <- range(c(dat1,dat2))+c(0,.5)
  plot(rep(1,length(cids))+jitz1,dat1,pch=".",
       col="blue",xlim=c(0.5,2.5),cex=cx1,ylim=yl,,xaxt="n",main="LIME1: Chr20:61.83Mb-61.85Mb",
       ylab="LRR-intensity (PC-corrected)",xlab="samples (by phenotype)",bty="l",cex.lab=1.5)
  points(rep(2,length(tids))+jitz2,dat2,pch=".",col="darkgreen",cex=cx2)
  legend("top",legend=c("controls","cases"),col=c("blue","green"),
         pch=".",pt.cex=6, bty="n",ncol=2,cex=2)
  
  ## add QS to cnvrs
  qs.fn1 <- "/chiswick/data/ncooper/immunochipRunTwo/RESULTS/Runz/RUN19/qs.del.resultsRUN19.txt"
  qs.fn2 <- "/chiswick/data/ncooper/immunochipRunTwo/RESULTS/Runz/RUN19/qs.dup.resultsRUN19.txt"
  qs1 <- reader(qs.fn1); qs2 <- reader(qs.fn2)
  print(load(getSlot(read.data.tracker(dir),"cnvresults")))
  cnvResults[[4]][["score"]] <- qs1[,1]
  cnvResults[[5]][["score"]] <- qs2[,1]
  cnvs1 <- add.scores.to.cnvrs(cnvs1,cnvResults[[4]])
  cnvs2 <- add.scores.to.cnvrs(cnvs2,cnvResults[[5]])
  results1 <- full.cnvr.summary(cnvs1,dir,thr=c(.5,.75,.9))
  results2 <- full.cnvr.summary(cnvs2,dir,thr=c(.5,.75,.9))
  
}
