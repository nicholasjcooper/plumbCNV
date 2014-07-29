

# convert Z scores to p values
# this should be in iFunctions
# Z.to.p <- function(Z) { p <- 2*pnorm(-abs(Z)); p[!is.finite(p)] <- NA; return(p) }

# conduct a test of proportions
t.o.p <- function(p1,p2,n1,n2) {
 p <- (p1 * n1 + p2 * n2) / (n1 + n2)
 SE <- sqrt((p*(1 - p)) * ((1/n1) + (1/n2)))
 z <- (p1 - p2) / SE
 return(list(Z=z,p=Z.to.p(z)))
}


## make a modified version of the top table that takes into account quality scores and 'decline'
compile.qs.results.to.cc.toptable <- function(qs.results,dir,suffix,cnvResult,decline.thresh=-.1) {
  results1 <- qs.results$DEL
  results1 <- within(results1,{decline <- rowMeans(cbind(c(qc90_cs-qc50_cs)/qc50_cs, c(qc90_ct-qc50_ct)/qc50_ct),na.rm=T)})
  print(head(results1[order(results1$qc75_sig),],15))
  results2 <- qs.results$DUP
  results2 <- within(results2,{decline <- rowMeans(cbind(c(qc90_cs-qc50_cs)/qc50_cs, c(qc90_ct-qc50_ct)/qc50_ct),na.rm=T)})
  print(head(results2[order(results2$qc75_sig),],15))
  save(qs.results,file=cat.path(dir$res,"qs.filt.results",suf=suffix,ext="RData"))

  tt1 <- (results1[order(results1$qc90_sig),])
  tt2 <- (results2[order(results2$qc75_sig),])
  oo1 <- cnvResult[[5]][[1]][which(rownames(cnvResult[[5]][[1]]) %in% (rownames(tt1[tt1$decline>= decline.thresh,]))),]
  oo2 <- cnvResult[[5]][[2]][which(rownames(cnvResult[[5]][[2]]) %in% (rownames(tt2[tt2$decline>= decline.thresh,]))),]
  oo1$sig  <- results1$qc90_sig[match(rownames(oo1),rownames(results1))]
  oo2$sig  <- results2$qc75_sig[match(rownames(oo2),rownames(results2))]
  oo1$cases  <- results1$qc90_cs[match(rownames(oo1),rownames(results1))]
  oo1$ctrls  <- results1$qc90_ct[match(rownames(oo1),rownames(results1))]
  oo2$ctrls  <- results2$qc75_ct[match(rownames(oo2),rownames(results2))]
  oo2$cases  <- results2$qc75_cs[match(rownames(oo2),rownames(results2))]
  toptables(oo1,oo2)
  return(list(DEL=oo1,DUP=oo2))
}

# internal function called by family.check.and.validate
# prints statistics checking transmission/denovo rates in mothers/ fathers, affected/unaffected, etc
# input parameters are quite specific objects
print.family.check <- function(tt1,tt2,mothrate,fathrate,denov_ctrl,denov_case,trans.del,trans.dup,ZZ1,ZZ2,ZZ3,ZZ4,ZZ5,ZZ6) {
  cat("Controls, denovo rate [vs trans]: DELs, DUPs",denov_ctrl,"\n") # ctrl denovo rate
  cat("Cases, denovo rate [vs trans]: DELs, DUPs",denov_case,"\n") # overall denovo rate for affected
  cat("affected vs control transmission rate for DELs",trans.del[1:2],"\n") # affected vs control transmission rate for dels
  cat("affected vs control transmission rate for DUPs",trans.dup[1:2],"\n") # '' '' for dups
  ZZ1 <- t.o.p(as.numeric(trans.del[1]),as.numeric(trans.del[2]),sum(tt1[17,]),sum(tt1[9,]-tt1[17,]))$Z
  ZZ2 <- t.o.p(as.numeric(trans.dup[1]),as.numeric(trans.dup[2]),sum(tt2[17,]),sum(tt2[9,]-tt2[17,]))$Z
  ZZ1a <- trans.del[3:4]; ZZ2a <- trans.dup[3:4]
  ZZ3 <- t.o.p(denov_case[1],denov_ctrl[1],sum(tt1[17,]),sum(tt1[9,]-tt1[17,]))$Z #DEL
  ZZ4 <- t.o.p(denov_case[2],denov_ctrl[2],sum(tt1[17,]),sum(tt1[9,]-tt1[17,]))$Z #DUP
  ZZ5 <- t.o.p(mothrate[1],fathrate[1],sum(tt1[3,]),sum(tt1[4,]))$Z #DEL
  ZZ6 <- t.o.p(mothrate[2],fathrate[2],sum(tt1[3,]),sum(tt1[4,]))$Z #DUP
  cat(paste("Transmissions (case v ctrl) rDELs: Z=",round(ZZ1,4),",p=",round(Z.to.p(ZZ1),5)),"\n")
  cat(paste("Transmissions (case v ctrl) rDUPs: Z=",round(ZZ2,4),",p=",round(Z.to.p(ZZ2),5)),"\n")
  cat(paste("Transmissions (case v ctrl) rDELs: chi=",round(ZZ1a[1],4),",p=",round(Z.to.p(ZZ1a[2]),5)),"\n")
  cat(paste("Transmissions (case v ctrl) rDUPs: chi=",round(ZZ2a[1],4),",p=",round(Z.to.p(ZZ2a[2]),5)),"\n")
  cat(paste("Denovos (case v ctrl) rDELs: Z=",round(ZZ3,4),",p=",round(Z.to.p(ZZ3),5)),"\n")
  cat(paste("Denovos (case v ctrl) rDUPs: Z=",round(ZZ4,4),",p=",round(Z.to.p(ZZ4),5)),"\n")
  cat(paste("Transmissions (mums v dads) rDELs: Z=",round(ZZ5,4),",p=",round(Z.to.p(ZZ5),5)),"\n")
  cat(paste("Transmissions (mums v dads) rDUPs: Z=",round(ZZ6,4),",p=",round(Z.to.p(ZZ6),5)),"\n")
}


summary.of.cnv.change.with.trios <- function(result.trio,result.no.trio,DEL=TRUE,ids=NULL,silent=FALSE,thr=NULL) {
  if(DEL) { list.el <- 4 } else { list.el <- 5 }
  if(length(result.trio)!=5) { result.trio <- result.trio[[1]] }
  if(length(result.no.trio)!=5) { result.no.trio <- result.no.trio[[1]] }
  tr <- (result.trio[[list.el]])
  ntr <- (result.no.trio[[list.el]])

  if(!is.null(thr)) { 
    tr <- tr[tr$score>=thr[list.el-3],] 
    ntr <- ntr[ntr$score>thr[list.el-3],]
  }
  nids <- ntr$id
  trids <- tr$id

  if(!is.null(ids)) { 
    nids <- nids[nids %in% ids]
    trids <- trids[trids %in% ids]
  } 
  n.id.cnt <- tapply(nids,nids,length)
  tr.id.cnt <- tapply(trids,trids,length)

  all.ids <- unique(c(nids,trids))
  n.idz <- length(all.ids)
  all.res <- matrix(0,nrow=length(all.ids),ncol=2)
  colnames(all.res) <- c("no.trios","trios")
  rownames(all.res) <- all.ids
  all.res[match(names(n.id.cnt),rownames(all.res)),1] <- n.id.cnt
  all.res[match(names(tr.id.cnt),rownames(all.res)),2] <- tr.id.cnt
  #prv(all.res) 
  dz <- table(all.res[,1],all.res[,2])

  total.pre <- sum(all.res[,1])
  total.post <- sum(all.res[,2])
  incr.tp <- sum(dz[upper.tri(dz)])  # more TP
  decr.fn <- sum(dz[lower.tri(dz)])  # less FN
  rel.tp <- sum(diag(dz))  # reliable TP
  if(!silent){
   cat("samples with increase in true positives:",incr.tp," ", incr.tp/n.idz,"%\n")
   cat("samples with decrease in false positives:", decr.fn," ", decr.fn/n.idz,"%\n")
   cat("samples with consistent true positives:", rel.tp," ",rel.tp/n.idz,"%\n")
   cat("total CNVs before:",total.pre,"after: ",total.post,"\n")
  }
  add <- cut <- keep1 <- keep2 <- 0
  for (rr in 1:nrow(dz)) {
    for (cc in 1:ncol(dz)) {
      if(cc>rr) { add <- add+((cc-rr)*dz[rr,cc]) ; keep1 <- keep1 + ((rr-1)*dz[rr,cc]) }
      if(rr>cc) { cut <- cut+((rr-cc)*dz[rr,cc]) ; keep2 <- keep2 + ((cc-1)*dz[rr,cc])  }
    }
  }
  rel.tp2 <- sum((-1+(1:length(diag(dz))))*diag(dz))
  keep <- (rel.tp2+keep1+keep2)
  if(!silent){
   cat("increase in true positives:", add," ", add/total.pre,"%\n")
   cat("decrease in false positives:", cut," ", cut/total.pre,"%\n")
   cat("persistent true positives:", keep," ", keep/total.pre,"% of pre",keep/total.post,"% of post\n")
  }
  return(c(incr.tp,incr.tp/n.idz,decr.fn,decr.fn/n.idz,rel.tp,rel.tp/n.idz,total.pre,
   total.post,add,add/total.pre,cut,cut/total.pre,keep,keep/total.pre,keep/total.post))
}



full.analysis.of.withandwithout <- function(trios=c(52,54,50,47,45),no.trios=c(51,53,49,48,46),
              labels=c("RDUP", "RDEL", "RDEL+", "RDUP+", "normal"),silent=TRUE,thr=NA) {
  pp <- read.ped.file("~/Documents/necessaryfilesICHIPFam/t1dgc-pedfile-2011-08-05.tab")
  qq <- pp[pp$father!=0 & pp$mother!=0,]
  rr <- pp[pp$father==0 & pp$mother==0,]
  kidz <- qq$sample; adultz <- rr$sample
  if(is.null(labels)) { labels <- rep("",length(trios)) }
  big.table <- matrix(nrow=length(trios)*6,ncol=22)
  colnames(big.table) <- c("HMM","Type","Who","tr.tr","ntr.tr","tr.dnv","ntr.dnv","incr.tp.samp","its.pc","decr.fp.samp","dfs.pc",
                           "consis.tp.samp","cts.pc","CNV.bef","CNV.aft",
                            "incr.tp","it.pc","decr.fp","df.pc","pers.tp","pc.pre","pc.post")
  big.table <- as.data.frame(big.table,stringsAsFactors=FALSE)
  big.table[[1]] <- rep(labels,each=6)
  big.table[[2]] <- c(rep("rDELs",3),rep("rDUPs",3))
  big.table[[3]] <- rep(c("ALL","Kids","Parents"),2)
  cl <- 8:(ncol(big.table))
  for (cc in 1:length(trios)) {
    if(!silent) { Header(paste(labels[cc],trios[cc],"vs",no.trios[cc])) }
    load(paste("RESULTS/fullresult",trios[cc],".RData",sep=""))
    cnvResult.trio <- cnvResult
    load(paste("RESULTS/fullresult",no.trios[cc],".RData",sep=""))
    cnvResult.raw <- cnvResult
    del.qs.tr <- reader(cat.path("RESULTS","qs.del.results",suf=trios[cc],ext="txt"))
    del.qs.ntr <- reader(cat.path("RESULTS","qs.del.results",suf=no.trios[cc],ext="txt"))
    dup.qs.tr <- reader(cat.path("RESULTS","qs.dup.results",suf=trios[cc],ext="txt"))
    dup.qs.ntr <- reader(cat.path("RESULTS","qs.dup.results",suf=no.trios[cc],ext="txt"))
    cnvResult.trio[[1]][[4]][["score"]] <- del.qs.tr[[1]]
    cnvResult.trio[[1]][[5]][["score"]] <- dup.qs.tr[[1]]
    cnvResult.raw[[1]][[4]][["score"]] <- del.qs.ntr[[1]]
    cnvResult.raw[[1]][[5]][["score"]] <- dup.qs.ntr[[1]]
    ccm <-  c(6291,1514)/sum(c(6291,1514))  #c(.37,.63) #cae-control ratio
    # get denovo and transmission rates
    tt1.tr <- trans.tests(reader(cat.path("RESULTS","famstats_trios",suf=trios[cc],ext="RData"))$tt1)[c(1,2,7,8)]
    tt1.ntr <- trans.tests(reader(cat.path("RESULTS","famstats_trios",suf=no.trios[cc],ext="RData"))$tt1)[c(1,2,7,8)]
    tt2.tr <- trans.tests(reader(cat.path("RESULTS","famstats_trios",suf=trios[cc],ext="RData"))$tt2)[c(1,2,7,8)]
    tt2.ntr <- trans.tests(reader(cat.path("RESULTS","famstats_trios",suf=no.trios[cc],ext="RData"))$tt2)[c(1,2,7,8)]     
    tr.change1 <- round(c(sum(tt1.tr[1:2]*ccm) , sum(tt1.ntr[1:2]*ccm)),3)
    tr.change2 <- round(c(sum(tt2.tr[1:2]*ccm) , sum(tt2.ntr[1:2]*ccm)),3)
    dn.change1 <- round(c(sum(tt1.tr[3:4]*ccm) , sum(tt1.ntr[3:4]*ccm)),3)
    dn.change2 <- round(c(sum(tt2.tr[3:4]*ccm) , sum(tt2.ntr[3:4]*ccm)),3)
    big.table[((cc-1)*6)+1,4:5] <- tr.change1; big.table[((cc-1)*6)+1,6:7] <- dn.change1 
    big.table[((cc-1)*6)+4,4:5] <- tr.change2; big.table[((cc-1)*6)+4,6:7] <- dn.change2
    #return(cnvResult.raw)
    if(!silent) { print("rDELs") }
    big.table[((cc-1)*6)+1,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=T,silent=silent,thr=thr)
    if(!silent) { print("children") }
    big.table[((cc-1)*6)+2,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=T,ids=kidz,silent=silent,thr=thr)
    if(!silent) { print("parents") }
    big.table[((cc-1)*6)+3,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=T,ids=adultz,silent=silent,thr=thr)
    if(!silent) { print("rDUPs") }
    big.table[((cc-1)*6)+4,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=F,silent=silent,thr=thr)
    if(!silent) { print("children") }
    big.table[((cc-1)*6)+5,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=F,ids=kidz,silent=silent,thr=thr)
    if(!silent) { print("parents") }
    big.table[((cc-1)*6)+6,cl] <- summary.of.cnv.change.with.trios(cnvResult.trio, cnvResult.raw,DEL=F,ids=adultz,silent=silent,thr=thr)
  }
  return(big.table)
}



# calculate proper sum level stats for transmissions, etc
trans.tests <- function(tt) {
    tt <- rowSums(tt)
    trios.with <- (tt[12]+tt[13]-tt[20]) 
    trios.with1 <- (tt[3]+tt[4]-tt[21]) 
    trios.with2 <- trios.with1 - trios.with
    did.trans <-  (tt[16]) 
    did.trans1 <- tt[7]-tt[16]
    did.not <- trios.with-did.trans
    did.not1 <- trios.with2 - did.trans1
    case.chi <- ((did.not-did.trans)^2)/(did.not+did.trans)
    control.chi <- ((did.not1-did.trans1)^2)/(did.not1+did.trans1)
    case.p <- pchisq(case.chi,1,lower.tail=F)
    control.p <- pchisq(control.chi,1,lower.tail=F)
    case <- did.trans/(did.not+did.trans)
    control <- did.trans1/(did.not1+did.trans1)
    denovo.cont <- (tt[8]-tt[14])/((tt[7]+tt[8])-(tt[14]+tt[16]))
    denovo.case <- tt[14]/(tt[14]+tt[16])
    out <- c(case, control, case.chi, case.p, control.chi, control.p,denovo.case,denovo.cont)
    names(out) <- c("case.tr.rate","control.tr.rate","case.chi", "case.p",
       "control.chi", "control.p","denovo.case","denovo.control")
    return(out)
}

# runs then print statistics checking transmission/denovo rates in mothers/ fathers, affected/unaffected, etc
# famstats is part of a name of the file to store family analysis statistics in
family.check.and.validate <- function(dir,fam.stats="famstats",top.dels=NULL,top.dups=NULL,suffix="") {
  load(cat.path(dir$res,"TDT_results",suf=suffix,ext="RData"))
  fam.stats <- cat.path(dir$res,fam.stats,ext="RData")
  tdt3 <- add.all.ids(tdt3,ped,dir)
  tdt4 <- add.all.ids(tdt4,ped,dir)
  tt1 <- get.CNV.wise.inheritance.counts(tdt3,ped=ped) # del table
  tt2 <- get.CNV.wise.inheritance.counts(tdt4,ped=ped) # dup table
  ## display the TDT information for significant hits ##
  if(length(top.dels)>0){
    print(tt1[-c(23:28),top.dels])
    print(round(tt1[c(23:28),top.dels],5))
  }
  if(length(top.dups)>0){
    print(tt2[-c(23:28),top.dups])
    print(round(tt2[c(23:28),top.dups],5))
  }
  ###
  print(tt0 <- cbind(rowSums(tt1),rowSums(tt2))[-c(23:28),]) # display the table
  #ratios of interest - assumption checking.. looks like dups over-called?
  cat("                                        DELS,   DUPs\n")
  cat("overall transmission ratio for mothers:",mothrate <- tt0[1,]/tt0[3,],"\n") # 
  cat("overall transmission ratio for fathers:",fathrate <- tt0[2,]/tt0[4,],"\n") # 
  cat("overall transmission ratio to affected for mothers:",tt0[10,]/tt0[12,],"\n") # 
  cat("overall transmission ratio to affected for fathers:",tt0[11,]/tt0[13,],"\n") # 
  trans.dup <- trans.tests(tt2) #rowMeans( tt2[c(27:28),] ,na.rm=T)
  trans.del <- trans.tests(tt1) #rowMeans( tt1[c(27:28),] ,na.rm=T)
  denov_ctrl <- (tt0[8,]-tt0[14,])/((tt0[7,]+tt0[8,])-(tt0[14,]+tt0[16,]))
  denov_case <- tt0[14,]/(tt0[14,]+tt0[16,])
  save(tt1,tt2,mothrate,fathrate,denov_ctrl,denov_case,trans.del,trans.dup,file=fam.stats)
  cat("wrote objects to file",fam.stats,"\n")
  print.family.check(tt1,tt2,mothrate,fathrate,denov_ctrl,denov_case,trans.del,trans.dup)
}


# in the current directory, generate plots of all members of each family that has a CNV
# e.g, 
# dir <- make.dir("/chiswick/data/ncooper/ImmunochipFamilies")
# family.check.and.validate(dir,"famstats_trios",suffix="48")
# plot.each.family.for.cnv(dir,reg="S1",chromo=1,cnv.bounds=c(197158752,197170596),suffix=48)
plot.each.family.for.cnv <- function(dir,ped,reg="S1",chromo=1,cnv.bounds=c(197158752,197170596),DEL=T,suffix=48) {
  print(load(cat.path(dir$res,"TDT_results",suf=suffix,ext="RData")))
  #S272"  # "S6" --> tdt3
  #cnv.bounds <- c(197158752,197170596) # c(31468234,31561619)  #c(20656322,21229324)
  cnv.w <- diff(cnv.bounds)
  plot.window <- c(cnv.bounds[1]-cnv.w,cnv.bounds[2]+cnv.w) 
  if(DEL) {
    tt1 <- get.CNV.wise.inheritance.counts(tdt3,ped=ped)
    ttt <- (as.numeric(tdt3[,reg]))
    names(ttt) <- rownames(tdt3)
  } else {
    tt2 <- get.CNV.wise.inheritance.counts(tdt4,ped=ped)
    ttt <- (as.numeric(tdt4[,reg]))
    names(ttt) <- rownames(tdt4)
  }
  tt <- tt1
  tt[,which(tt[8,]>5)]
  ped[["cnv.reg"]] <- ttt[match(rownames(ped),names(ttt))]
  ped <- ped.interp(ped,F)  
  with(ped , table(father,affected,cnv.reg,exclude=NULL))
  fam.grps <- tapply(rownames(ped),factor(ped$familyid),c)
  s6.grps.list <- tapply(ped[,"cnv.reg"],factor(ped$familyid),c)
  s6.who.list <- tapply(ped$who,factor(ped$familyid),c)
  s6.grps <- tapply(ped[,"cnv.reg"],factor(ped$familyid),max)
  ff <- which(as.numeric(s6.grps)>1)
  repz <- length(fam.grps)
  sstt <- proc.time()
  for(cc in 1:repz) {
    if(!cc %in% ff) { next }
    labl <- paste(substr(fam.grps[[cc]][1],1,6),paste(s6.grps.list[[cc]],collapse=""),
                  paste(s6.who.list[[cc]],collapse="_"),sep="_")
    cnv.plot(Cnv=cnv.bounds,PREPOSTPC=T,Chr=chromo,Pos=plot.window,samples=fam.grps[[cc]],
             LRR=T,BAF=T,show.fail=T,dir=dir,cnvPlotFileName=paste(labl,"pdf",sep="."))
    loop.tracker(cc,repz,st.time=sstt)
  }
  return(NULL)
}

# add column to a ped file, showing in plain english who is who within families, mum dad, boy, girl, etc
ped.interp <- function(ped,long=TRUE) {
  want <- c("familyid","member","father","mother","sex","affected")
  sample.info <- read.sample.info(dir)
  if(!all(want %in% colnames(ped))) { stop("invalid ped file frame [use 'get.pedData()']") }
  long.codes <- c("father","mother","boy","girl","control","t1d")
  short.codes <- c("F","M","B","G","Ct","T1")
  if(long) { codes <- long.codes } else { codes <- short.codes }
  ped[["who"]] <- ""
  ped[["who"]][(is.na(ped$father) & ped$sex==1)] <- codes[1]
  ped[["who"]][(is.na(ped$father) & ped$sex==2)] <- codes[2]
  ped[["who"]][(!is.na(ped$father) & ped$sex==1)] <- codes[3]
  ped[["who"]][(!is.na(ped$father) & ped$sex==2)] <- codes[4]
  ped[["who"]] <- paste(ped[["who"]],codes[5:6][ped$affected],sep=".")
  return(ped)
}

# the core function to just calculate the TDT for a set of CNVs/CNVRs
core.tdt <- function(dir,cnvrs,cnvs,ped,double.table,DEL=TRUE) {
  cnvs <- update.cnvrs.with.cn(cnvs,double.table,DEL=DEL)
  tdt3 <- make.cnv.reg.snp.matrix(cnvs)
  tdt3 <- add.all.ids(tdt3, ped, dir)
  ii <- tdt.snp(data=ped,snp.data=tdt3)
  tt <- p.value(ii,1)
  cnvrs[["tdt"]][match(names(tt),rownames(cnvrs))] <- tt
  sum.del <- ranged.to.data.frame(cnvrs,T)[order(cnvrs[[5]]),]
  colnames(sum.del)[c(7,9)] <- c("p.fet","p.tdt")                              
  ## remove NAs..
  sum.del2 <- sum.del[!is.na(sum.del$p.tdt),]
  sum.del2[["genes"]] <- substr(sum.del2[["genes"]],1,10);
  return(list(cnvrs=cnvrs,sum.del=sum.del,tdt3=tdt3,ped=ped))
}


## main function to conduct an analysis using trios, TDT, etc
trio.analysis <- function(dir=NULL, cnvResults, ped.file, result.pref="",quality.scores=FALSE,restore=TRUE) {
  dir <- validate.dir.for(dir,c("res","cnv.qc"))
  if(!is.list(cnvResults) | length(cnvResults)!=5 | !is(cnvResults[[1]])[1]=="RangedData") { 
    warning("cnvResults object should be a list of RangedData objects, length 5"); return(NULL) }
  cat("\nRunning TDT analysis for family trios\n")
  oo1 <- extract.cnv.regions(dir,type="del",by.cnv=FALSE,lwr=0.25,upr=4,FET=T,prt=F) # regionlist
  oo2 <- extract.cnv.regions(dir,type="dup",by.cnv=FALSE,lwr=0.25,upr=4,FET=T,prt=F) # regionlist
  oo1[["tdt"]] <- NA; oo2[["tdt"]] <- NA
  
  DT <- read.data.tracker(dir)
  double.del.table <- cnvResults[[4]][which(cnvResults[[4]]$cn==0),] # rareDELs
  double.dup.table <- cnvResults[[5]][which(cnvResults[[5]]$cn==4),] # rareDUPs

  ped <- get.pedData(ped.file)

  if(quality.scores) {
    results <- process.quality.scores(DT,suffix=result.pref,dir,restore=restore)
    oo3 <- results[[3]]
    oo4 <- results[[4]]
    oo3 <- remove.duplicated.id.ranges(oo3)
    oo4 <- remove.duplicated.id.ranges(oo4)
    outlist1 <- core.tdt(dir,cnvrs=oo1,cnvs=oo3,ped,double.del.table,DEL=TRUE)
    sum.del <- outlist1$sum.del
    outlist2 <- core.tdt(dir,cnvrs=oo2,cnvs=oo4,ped,double.dup.table,DEL=FALSE)
    sum.dup <- outlist2$sum.del
    keep.lev.dup <- 2; keep.lev.del <- 2
    q.levs <- c(0.5,0.75,0.95)
    colz <- c("cases", "ctrls", "p.tdt")
    cn <- matrix(ncol=length(colz),nrow=length(q.levs))
    for (jj in 1:length(q.levs)) {
      cc <- q.levs[jj]
      res1a <- filt.counts.cnvr(oo1,oo3,col="score",thresh=cc,excl.less.than=T)
      res2a <- filt.counts.cnvr(oo2,oo4,col="score",thresh=cc,excl.less.than=T)
      oo1a <- res1a$cnvr; oo3a <- res1a$cnv; oo2a <- res2a$cnvr; oo4a <- res2a$cnv
      outlist1a <- core.tdt(dir,cnvrs=oo1a,cnvs=oo3a,ped,double.del.table,DEL=TRUE)
      outlist2a <- core.tdt(dir,cnvrs=oo2a,cnvs=oo4a,ped,double.dup.table,DEL=FALSE) #sum.del even for Dups!
      if(jj==keep.lev.del) {  oo1 <- oo1a;  tdt3 <- outlist1a$tdt3 }
      if(jj==keep.lev.dup) {  oo2 <- oo2a ; tdt4 <- outlist2a$tdt3 }
      for (dd in 1:length(colz)) {
         cn[jj,dd] <- paste(colz[dd],cc,sep="_")
         # rint(cn)
         sum.del[[cn[jj,dd]]][match(rownames(outlist1a$sum.del),rownames(sum.del))] <- outlist1a$sum.del[[colz[dd]]] 
         sum.dup[[cn[jj,dd]]][match(rownames(outlist2a$sum.del),rownames(sum.dup))] <- outlist2a$sum.del[[colz[dd]]]     
      }
    }
    cs.pc <-(sum.del[,cn[3,1]] - sum.del[,cn[1,1]])/sum.del[,cn[1,1]]
    ct.pc <-(sum.del[,cn[3,2]] - sum.del[,cn[1,2]])/sum.del[,cn[1,2]]
    sum.del[["decline"]] <- rowMeans(cbind(cs.pc,ct.pc),na.rm=T) 
    cs.pc <-(sum.dup[,cn[3,1]] - sum.dup[,cn[1,1]])/sum.dup[,cn[1,1]]
    ct.pc <-(sum.dup[,cn[3,2]] - sum.dup[,cn[1,2]])/sum.dup[,cn[1,2]]
    sum.dup[["decline"]] <- rowMeans(cbind(cs.pc,ct.pc),na.rm=T) 
    del.p <- cn[keep.lev.del,3]; dup.p <- cn[keep.lev.dup,3]
  } else {
    oo3 <- extract.cnv.regions(dir,type="del",by.cnv=TRUE,lwr=0.25,upr=4,FET=T,prt=F) # cnv list
    oo4 <- extract.cnv.regions(dir,type="dup",by.cnv=TRUE,lwr=0.25,upr=4,FET=T,prt=F) # cnv list
    outlist <- core.tdt(dir,cnvrs=oo1,cnvs=oo3,ped,double.del.table,DEL=TRUE)
    oo1 <- outlist$cnvrs; sum.del <- outlist$sum.del; tdt3 <- outlist$tdt3
    #(list(cnvrs==cnvrs,sum.del=sum.del,tdt3=tdt3,ped=ped))
    outlist <- core.tdt(dir,cnvrs=oo2,cnvs=oo4,ped,double.dup.table,DEL=FALSE)
    oo2 <- outlist$cnvrs; sum.dup <- outlist$sum.del; tdt4 <- outlist$tdt3
    dup.p <- del.p <- "p.tdt"
  } 
  ofn=cat.path(dir$res,"TDT_results",suf=result.pref,ext="RData")
  save(oo1,oo2,sum.del,sum.dup,tdt3,tdt4,ped,file=ofn)
  CNVR <- list(deletions=oo1,duplications=oo2)
  cat(" wrote primary TDT family analysis results to",ofn,"\n")
  ofn <- cat.path(dir$res,"TDTsummary",suf=result.pref,ext="txt")
  sum.del2 <- sum.del[!is.na(sum.del[,del.p]),]
  sum.dup2 <- sum.dup[!is.na(sum.dup[,dup.p]),]
  sum.del2[["genes"]] <- substr(sum.del2[["genes"]],1,14)
  sum.dup2[["genes"]] <- substr(sum.dup2[["genes"]],1,14)
  sink(ofn)
   cat("\nDeletions TDT Results\n") ; print(sum.del2[sum.del2[,del.p]<.05,])
   cat("\nDuplications TDT Results\n") ; print(sum.dup2[sum.dup2[,dup.p]<.05,])  # to file copy
   top.dels <- rownames(cnvResults$cnvr[[1]][narm(which(cnvResults$cnvr[[1]]$tdt<.01)),])
   top.dups <- rownames(cnvResults$cnvr[[2]][narm(which(cnvResults$cnvr[[2]]$tdt<.01)),])
   #load(cat.path(dir$res,"TDT_results",suf=suffix,ext="RData"))
   fam.fn <- cat.path(dir$res,"famstats_trios",suf=result.pref,ext="RData")
   cat("\nTDT Counts\n") #; print(sum.del[sum.del[,del.p]<.05,])
   family.check.and.validate(dir,fam.fn,top.dels=top.dels,top.dups=top.dups,suffix=result.pref)
  sink()
  cat("\nDeletions TDT Results\n") ; print(sum.del2[sum.del2[,del.p]<.05,])
  cat("\nDuplications TDT Results\n") ; print(sum.dup2[sum.dup2[,dup.p]<.05,]) 
  load(fam.fn)
  print.family.check(tt1,tt2,mothrate,fathrate,denov_ctrl,denov_case,trans.del,trans.dup)
  cat("wrote TDT count objects to file",fam.fn,"\n")
  # on screen copy
  cat("wrote summary to file:",ofn,"\n")

  return(CNVR)
}


## add all ids in the study to CNV/snpMatrix object containing only some (e.g, passing qc, have CNV, etc)
add.all.ids <- function(tdt.cnv, ped, dir) {
  dir <- validate.dir.for(dir,"ano")
  want <- c("familyid","member","father","mother","sex","affected")
  sample.info <- read.sample.info(dir)
  if(!all(c("phenotype","QCfail") %in% colnames(sample.info))) { stop("invalid sample.info file in 'dir$ano'") }
  if(!all(want %in% colnames(ped))) { stop("invalid ped file frame [use 'get.pedData()']") }
  cur.ids <- rownames(tdt.cnv)
  all.ids <- rownames(sample.info)
  all.ids <- all.ids[all.ids %in% rownames(ped)]
  #new.ids.pass <- all.ids[!all.ids %in% cur.ids & sample.info$QCfail==0]
  new.ids.fail <- all.ids[(!all.ids %in% cur.ids) & sample.info$QCfail==1] # likely empty
  old.ids.fail <- all.ids[(all.ids %in% cur.ids) & sample.info$QCfail==1] # these should already be NAs
  CNVRs <- colnames(tdt.cnv)
  out <- as.data.frame(matrix(0,nrow=length(all.ids),ncol=length(CNVRs)))
  rownames(out) <- paste(all.ids); colnames(out) <- paste(CNVRs) 
  if(length(new.ids.fail)>0) { out[new.ids.fail,] <- NA }
  out <- data.frame.to.SnpMatrix(out)
  out[rownames(tdt.cnv),] <- tdt.cnv[rownames(tdt.cnv),]
  #assume all new.ids.pass == 1 [no need to do it explicitly though]
  return(out)
}


#make a snp matrix that represents a set of CNVs
make.cnv.reg.snp.matrix <- function(X) {
  # imagine a CNV region... call the CNV (copy 1/3) as heterozygous = 2
  #                               Normal (copy 2)   as homozygous  =  1
  #                                  CNV (copy 0/4+) as homozygous = 3
  if(is(X)[1]=="list") { stop("X was a list, should be matrix/data.frame/RangedData... etc") }
  if(is.null(dim(X))) { stop("X must be 2 dimensional") }
  if(length(Dim(X))>2) { stop("X must be 2 dimensional") }
  colnames(X) <- tolower(colnames(X))
  if(!all(c("id","phenotype","cnvr") %in% colnames(X))) { stop("Invalid X, must have columns named, 'id', 'phenotype' and 'cnvr'") }
  all.ids <- sort(unique(X$id))
  pheno <- X$phenotype[match(all.ids,X$id)]
  CNVRs <- unique(X$cnvr)
  CNVRs <- CNVRs[order(as.numeric(gsub("S","",CNVRs)))]
  out <- as.data.frame(matrix(0,nrow=length(all.ids),ncol=length(CNVRs)))
  rownames(out) <- paste(all.ids); colnames(out) <- paste(CNVRs) 
  id.lists <- tapply(X$id,factor(X$cnvr),c)
  if(!all(c("copy") %in% colnames(X))) {
    for (cc in 1:length(id.lists)) {
      out[match(id.lists[cc],rownames(out)),which(colnames(out)==names(id.lists)[cc])] <- 2
    }
  } else {
    if(any(names(table(X$copy))==0)) { doubler <- 0 } else { doubler <- 4 } # which copy del/dup?
    #print(doubler)
    Xd <- X[X$copy==doubler,] # which sample/CNV combos have copy 0/4
    #prv(Xd)
    for (cc in 1:length(id.lists)) {
      samps <- narm(match((id.lists[[cc]]),rownames(out)))
      cnvz <- which(colnames(out)==names(id.lists)[cc])
      out[samps,cnvz] <- 1 # heterozygous = 1 dup/dep
      #print(Dim(out[samps,cnvz]))
      #if(cc>500) { prv(samps,cnvz,out) }
      doubleidsinsamps <- Xd$id %in% rownames(out)[samps]
      doublecnvsinlist <- Xd$cnvr %in% names(id.lists)[cc]
      rel.rows <- which(doubleidsinsamps & doublecnvsinlist)
      #prv(doubleidsinsamps,doublecnvsinlist,rel.rows)
      if(length(rel.rows)>0) {
        samps <- narm(match(Xd$id[rel.rows],rownames(out)))
        out[samps,cnvz] <- 2 # homozygous = 2*dup/2*del
      }
      #loop.tracker(cc,length(id.lists))
    }
  }
  return(data.frame.to.SnpMatrix(out))
}



## extract family links from a ped file
get.ped.linked.sets <- function(ped,tdt3) {
  want <- c("familyid","member","father","mother","sex","affected")
  if(!all(want %in% colnames(ped))) { stop("invalid ped file frame [use 'get.pedData()']") }
  #which are the route cases
  ind <- which(with(ped,father==1 & mother==2 & affected==2))
  #which of these cases actually in the cnv dataset
  sel <- which(rownames(ped)[ind] %in% rownames(tdt3))
  #which family numbers to look at for potential valid tdt parents
  famz <- ped[ind[sel],"familyid"]
  #get dataset of key families
  ww <- ped[ped$familyid %in% famz,]
  #choose just the mothers and fathers
  parz <- which(is.na(ww$father) & is.na(ww$mother))
  #choose just the children
  kidz <- which(!is.na(ww$father) & !is.na(ww$mother))
  # fix errors in ped file coding
  badz <- with(ped,which(father==1 & mother==1))
  if(length(badz)>0) { ped[badz,"mother"] <- 2 }
  badz <- with(ped,which(father==2 & mother==2))
  if(length(badz)>0) { ped[badz,"father"] <- 1 }
  #create dataset of kids with link to datasets of parents
  kk <- ww[kidz,] # just the children
  mm <- ww[which(is.na(ww$father) & is.na(ww$mother) & ww$sex==2),]  #mothers
  ff <- ww[which(is.na(ww$father) & is.na(ww$mother) & ww$sex==1),]  #fathers
  kmlink <- match(kk$familyid,mm$familyid) # make an index to match each child to its mother in mm
  kflink <- match(kk$familyid,ff$familyid) # make an index to match each child to its father in ff
  matchlist <- tapply(1:nrow(kk),factor(kk$familyid),c) # make a list of rownumbers for children in each family
  mklink <- matchlist[match(mm$familyid,names(matchlist))] # make an index to match mother to her children in kk
  fklink <- matchlist[match(ff$familyid,names(matchlist))] # make an index to match father to his children in kk
  return(list(kk=kk,mm=mm,ff=ff,kmlink=kmlink,kflink=kflink,mklink=mklink,fklink=fklink))
}


# internal, just until NCmisc library is updated
list.to.env <- function(list) {
  if(!is.list(list)) { stop("this function's sole parameter must be a list object")}
  if(is.null(names(list))) { stop("list elements must be named") }
  if(length(list)>1000) { warning("list contains over 1000 elements, this operation will crowd the workspace") }
  for(cc in 1:length(list)) {
    assign(x=names(list)[cc],value=list[[cc]],pos=parent.frame())
  }
  return(NULL)
}


## generate table of counts, transmissions, etc for each CNV that is used for TDT and other transmission analysis
get.CNV.wise.inheritance.counts <- function(tdt.snp,ped=NULL,only.doubles=FALSE,
                                            replace.na=FALSE,replace.with=c(M=0.5,D=0.5,C=0.5)) {
  if(is(tdt.snp)[1]=="SnpMatrix") {  TDT <- SnpMatrix.to.data.frame(tdt.snp) } else { TDT <- tdt.snp }
  if(!is.data.frame(TDT)) { stop("tdt.snp must be a SnpMatrix or data.frame with snp-ids as rownames") }
  if(is.null(ped)) { stop("a valid 'pedData' object must be inputted (see snpStats:tdt.snp documentation)")}
 # prv(ped)
  ped.list <- get.ped.linked.sets(ped,tdt.snp)
  if(replace.na) { if(length(replace.with)!=3 | !all(is.numeric(unlist(replace.with)))) { 
    stop("replace.with must be a vector of three scalars for mother, father, child missing values")} }
  list.to.env(ped.list) # take the variables from the list and assign them into the local environment
  if(!exists("kk") | !exists("fklink")) { stop("failure to import variables from get.ped.linked.sets() function") }
  RN  <- c("mother-pass","father-pass","mother-kidcount","father-kidcount","total-pass",
                           "total-kidcount","in-a-parent","denovo","parent-count",
                           "aff-mother-pass","aff-father-pass","affect-mum-kidcount","affect-dad-kidcount",
                           "affected-denovos","affected-kidcount2","aff-in-a-parent","aff-parent-count",
                           "in.both","aff-in-both","double-counted-aff","double-counted-kids","total-aff-kidcount","chisq.tdt","p.tdt","chisq.ctrl","p.ctrl","aff-pc-trans","ctrl-pc-trans")
  countmat <- matrix(0,nrow=length(RN),ncol=ncol(TDT))
  rownames(countmat) <- RN
  colnames(countmat) <- colnames(TDT)
  affect.samps <- rownames(ped)[ped$affected==2]
  if(only.doubles) { thresh <- 1 } else { thresh <- 0 }
  counted.kids <- counted.affs <- vector("list",ncol(TDT))
  ## MUMS ##
  for(cc in 1:nrow(mm)) {
    loop.tracker(cc,2*nrow(mm)+nrow(kk))  # allows loop tracker to span all three loops
    dels <- which(TDT[rownames(mm)[cc],]>thresh) # which DELs does the next mum have
    if(length(dels)<1) { next }  # skip if this mother has none
    their.kids <- rownames(kk)[mklink[[cc]]]
    if(length(their.kids)<1) { next }
    mat <- TDT[their.kids,dels,drop=FALSE]  # extract her children for these DEL snps
    mat[mat>1] <- 1
    mat2 <- mat[rownames(mat) %in% affect.samps,,drop=FALSE] # only with affected kids
   # if("S1" %in% colnames(mat)) { print(cc);print(mat); print(mat2) }
    if(replace.na) { mat[is.na(mat)] <- force.percentage(replace.with[[1]]) }
    if(length(Dim(mat))<2) {
      snp.counts <- mat  # in case only 1 child may also be a row, not a matrix
      snp.counts2 <- mat2 # (affected)
    } else {
      snp.counts <- colSums(mat,na.rm=T) # add the number of times passed to children to the snp count
      snp.counts2 <- colSums(mat2,na.rm=T) # (affected)
    }
    countmat[1,dels] <- countmat[1,dels] + snp.counts
    countmat[10,dels] <- countmat[10,dels] + snp.counts2  # affected
    countmat[3,dels] <- countmat[3,dels] + nrow(mat)
    countmat[12,dels] <- countmat[12,dels] + nrow(mat2)
    if(length(dels)==0) { stop("length of jj was 0") }
    for(jj in 1:length(dels)) {
      # keep track of which affected kids are added to the 'passed on' count, in case doubleup in dads
      more.kids <- rownames(mat[mat[,jj]>0,])
      more.affs <- rownames(mat2[mat2[,jj]>0,])
      if(length(more.affs)>0) { counted.affs[[dels[jj]]] <- c(counted.affs[[dels[jj]]],more.affs) }
      if(length(more.kids)>0) { counted.kids[[dels[jj]]] <- c(counted.kids[[dels[jj]]],more.kids) }
    }
  }
  ## DADS ##
  for(cc in 1:nrow(ff)) {
    loop.tracker(cc+1,2*nrow(ff)+nrow(kk)+1)  # allows loop tracker to span all three loops
    dels <- which(TDT[rownames(ff)[cc],]>thresh) # which DELs does the next dad have
    if(length(dels)<1) { next }  # skip if this father has none
    their.kids <- rownames(kk)[fklink[[cc]]]
    if(length(their.kids)<1) { next }
    mat <- TDT[their.kids,dels,drop=FALSE]  # extract his children for these DEL snps
    mat[mat>1] <- 1
    mat2 <- mat[rownames(mat) %in% affect.samps,,drop=FALSE] # only with affected kids
    #if("S1" %in% colnames(mat)) { print(cc);print(mat); print(mat2) }
    if(replace.na) { mat[is.na(mat)] <- force.percentage(replace.with[[2]]) }
    if(length(Dim(mat))<2) {
      snp.counts <- mat  # in case only 1 child, a row, not a matrix
      snp.counts2 <- mat2 # (affected)
    } else {
      snp.counts <- colSums(mat,na.rm=T) # add the number of times passed to children to the snp count
      snp.counts2 <- colSums(mat2,na.rm=T) # (affected)
    }
    countmat[2,dels] <- countmat[2,dels] + snp.counts
    countmat[11,dels] <- countmat[11,dels] + snp.counts2  # affected
    countmat[4,dels] <- countmat[4,dels] + nrow(mat)
    countmat[13,dels] <- countmat[13,dels] + nrow(mat2)
    if(length(dels)==0) { stop("length of jj was 0") }
    for(jj in 1:length(dels)) {
      again.affs <- which(rownames(mat2[mat2[,jj]>0,]) %in% counted.affs[[dels[jj]]])
      again.kids <- which(rownames(mat[mat[,jj]>0,]) %in% counted.kids[[dels[jj]]])
      countmat[21,dels[jj]] <- countmat[21,dels[jj]] + length(again.kids)
      countmat[20,dels[jj]] <- countmat[20,dels[jj]] + length(again.affs)
      # add the number of aff kids already counted in the mum to this row total
    }
  }
  countmat[5,] <- countmat[1,]+countmat[2,] # add passed-on CNVs from mum and dad to 'total' row
  countmat[6,] <- countmat[3,]+countmat[4,] # add child count from mum and dad to 'total' row
  
  ## now all in reverse, looking up the parents for each child ##
  
  ## KIDS ##
  for(cc in 1:nrow(kk)) {
    loop.tracker(cc+1,2*nrow(ff)+nrow(kk)+1) # allows loop tracker to span all three loops
    dels <- which(TDT[rownames(kk)[cc],]>thresh) # which DELs does the next kid have
    if(length(dels)<1) { next }  # skip if this mother has none
    their.folks <- c(rownames(mm)[kmlink[cc]],rownames(ff)[kflink[cc]])
    #cat("processing child:",rownames(kk)[cc],"with parents",paste(their.folks, collapse=","),"\n")
    if(length(their.folks)<1) { next }
    mat <- TDT[their.folks,dels,drop=FALSE]  # extract mother+father of child for these DEL/DUP snps
    #if("S6" %in% colnames(mat)) { print(mat) }
    if(replace.na) { mat[is.na(mat)] <- force.percentage(replace.with[[3]]) }
   # if("S1" %in% colnames(mat)) { print(cc);print(mat); print(mat2) }
    if(length(Dim(mat))<2) {
      snp.counts <- mat  # in case only 1 parent, a row, not a matrix
    } else {
      snp.counts <- colSums(mat,na.rm=T) 
    }
    in.both <- as.numeric(snp.counts>1)
    countmat[18,dels] <- countmat[18,dels] + in.both # how many times in both parents
    snp.counts[snp.counts>1] <- 1 # add '1' if at least one parent had the same CNV found in the child
    #if(!all(snp.counts %in% c(0,1))) { warning("snp.counts seems invalid: ",paste(snp.counts,collapse=",")) }
    countmat[7,dels] <- countmat[7,dels] + snp.counts # increments snps if at least 1 parent has the cnv
    countmat[8,dels] <- countmat[8,dels] + (1-snp.counts) # adds denovos for each snp
    countmat[9,dels] <- countmat[9,dels] + nrow(mat) # adds to parent count for each snp
    if(rownames(kk)[cc] %in% affect.samps) {
      countmat[19,dels] <- countmat[19,dels] + in.both
      countmat[14,dels] <- countmat[14,dels] + (1-snp.counts) # increments snps for aff if at least 1 parent has the cnv
      countmat[15,dels] <- countmat[15,dels] + 1 # alternate way to calculate # of affected kids for each cnv
      countmat[16,dels] <- countmat[16,dels] + snp.counts # affected kids count, when a parent has cnv
      countmat[17,dels] <- countmat[17,dels] + nrow(mat) # affected kids parent count
    }
    ## run TDT tests ##
    trios.with <- countmat[22,] <- (countmat[12,]+countmat[13,]-countmat[20,]) # number of trios with an affect child where 1 parent has the cnv
    trios.with1 <- countmat[6,] <- (countmat[3,]+countmat[4,]-countmat[21,]) # number of trios with an 
    trios.with2 <- trios.with1 - trios.with
    did.trans <-  (countmat[16,]) # number of transmissions to an affected child where 1 parent has cnv [subtract where 2 parents have]
    did.trans1 <- countmat[7,]-countmat[16,]
    did.not <- trios.with-did.trans
    did.not1 <- trios.with2 - did.trans1
    tdt.chi <- ((did.not-did.trans)^2)/(did.not+did.trans)
    tdt.chi1 <- ((did.not1-did.trans1)^2)/(did.not1+did.trans1)
    tdt.p <- pchisq(tdt.chi,1,lower.tail=F)
    tdt.p1 <- pchisq(tdt.chi1,1,lower.tail=F)
    countmat[23,] <- tdt.chi
    countmat[24,] <- tdt.p
    countmat[25,] <- tdt.chi1
    countmat[26,] <- tdt.p1
    countmat[27,] <- did.trans/(did.not+did.trans)
    countmat[28,] <- did.trans1/(did.not1+did.trans1)
  }
  cat("\n")
  ###
  return(countmat)
}







#play with TDT


if(F & !exists("tdt3")) {
  # only bother recalculating if these vars not already present
  #print(load("/chiswick/data/ncooper/ImmunochipFamilies/RESULTS/TDT_results.RData"))
  #print(load(cat.path(dir$res,"TDT_results",suf=suffix,ext="RData")))
  dir <- make.dir("/chiswick/data/ncooper/ImmunochipFamilies")
  
  dir <- make.dir("/chiswick/data/ncooper/ImmunochipFamilies")
  family.check.and.validate(dir,"famstats_trios",suffix="")
  plot.each.family.for.cnv(dir,reg="S1",chromo=1,cnv.bounds=c(197158752,197170596),suffix=48)
  #tdt3 <- add.all.ids(tdt3,ped,dir)
  #tdt4 <- add.all.ids(tdt4,ped,dir)
  #add in the snp data # note that missing will get coded as 0 (zero) in this step
  #  ww[["S51"]] <- as.numeric(tdt3[,"S51"][match(rownames(ww),rownames(tdt3))])
  
  # key families are those which have a parent with the cnv
  #  keyfams <- names(which(tapply(ww$S51[parz],factor(ww$familyid[parz]),function(X) { any(X==2) })))
  #just the key families
  #  with(ww[ww$familyid %in% keyfams,],table(affected,S51))
  # key families are those which have a child with the cnv
  #  keyfams2 <- names(which(tapply(ww$S51[kidz],factor(ww$familyid[kidz]),function(X) { any(X==2) })))
  #look at just the families with affected kids but not parents
  #  with(ww[ww$familyid %in% keyfams2 & !ww$familyid %in% keyfams,],table(affected,S51)) 
}




# ONLY USE below to make plots of families for the CNV 'S6'
# DN --\
# ---NN, NN, DN, DN
# NN --/ 
#   
#   DN --\:q

# ---DD, DN, DN, NN
# DN --/ 
#   
#   
#   105889 mother has, 1 af child dont, 1 af does, 1 unaf dont
# 157828 father has, 1 af does, 1 af qc-ed
# 162020 mother has, 1 af child dont,  1 af qc-ed, 1 unaf qc-ed
# 231234 father has, 1 af child dont,  1 af qc-ed, 1 unaf qc-ed
# 263224 father has, 2 af child dont, 1 unaf qc-ed
# 282433 father has, 1 unaf has, 1 af doesn't, 1 af qc-ed
# 297972 mother has, 1 af has, 1 af qc-ed
# 298496 mother has, 1 af dont, 1 af qc-ed
# 468821 father has, 1 unaf dont, 2 af dont
# 
# 1 1
#   1
# 1
# 1
# 2
# 1
#   1
# 1
# 2
# 
# 
# 9 dont, 3 do
# 
# 36/12 = 3


# want to know the general rate of transmission for all DELs/DUPs
# want to know the denovo rate of DELs/DUPs
# 
# 8:1 'denovos' are affected
# next step, count within each family how many times each CNV is passed 
# on versus not, sum for all. get transmission rate.

