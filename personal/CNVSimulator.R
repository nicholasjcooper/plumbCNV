sim.dir <- make.dir("/chiswick/data/ncooper/simRunTest/")
#init.dirs.fn(sim.dir)

#play with sim functions
if(F) {
plot(sim.snp.baf(1000,"dup4",noise=.5))
plot(sim.snp.lrr(1000,"dup4",noise=1))
plot(sim.snp.baf(1000,"dup3",noise=1.5))
plot(sim.snp.lrr(1000,"dup3",noise=4))
plot(sim.snp.baf(1000,"del1",noise=.5))
plot(sim.snp.lrr(1000,"del1",noise=1.75))
plot(sim.snp.baf(1000,"del0",noise=.5))
plot(sim.snp.lrr(1000,"del0",noise=.5))
plot(sim.snp.baf(1000,"loh",noise=1.15))
plot(sim.snp.lrr(1000,"loh",noise=.95))
plot(sim.snp.baf(1000,"normal",noise=1))
plot(sim.snp.lrr(1000,"normal",noise=.75))

}


  
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")

dir <- make.dir("/chiswick/data/ncooper/immunochipRunTest/")
DT <- read.data.tracker(dir)
bln <- getSlot(DT,"big.lrr")
bl <- getBigMat(bln[[1]],dir=dir$big)
uu <- get.plate.info.for.big(bl,dir)

# yields matrices of 49 plates x 196K snps
vv <- get.plate.snp.stats.for.big(bl,dir,func=mean,n.cores=22)
ww <- get.plate.snp.stats.for.big(bl,dir,func=sd,n.cores=22)

vv <- vv-vv
ww <- ww-ww

print(mean(vv)); print(mean(ww))

nsub <- nrow(vv)

# do this for each bigMat
# create new bigMat with this data, each plate like a sample
# forms basis to overlay some fake CNVs
# overlay CNVs with different strengths and noise levels
# see which are detected and which not


#setupframe
pfb <- reader("/chiswick/data/ncooper/immunochipRunTest/PENNCNV/BAF.pfb")

select <- match(rownames(pfb),colnames(vv))
nn <- c(0.5,0.75,1,1.25,1.5,1.75,2)
BAF <- LRR <- matrix(numeric(),nrow=nrow(pfb),ncol=nsub*length(nn))
rownames(BAF) <- rownames(LRR) <- rownames(pfb)
colnames(BAF) <- colnames(LRR) <- idz <- paste("n",rep(nn,each=nsub),rep(rownames(vv),length(nn)),sep="_")

idz.list <- tapply(idz,factor(rep(c(1:length(nn)),each=nsub)),list)

print(load("/chiswick/data/ncooper/immunochipRunTest/RESULTS/cnvResultsPCA2424.RData"))
deL <- toGenomeOrder(cnvResults[[4]],strict=T); duP <- toGenomeOrder(cnvResults[[5]],strict=T)
sim.dels <- sample(nrow(deL),nsub*2); sim.dups <- sample(nrow(duP),nsub*2)
sim.dels1 <- sort(sim.dels[1:nsub]); sim.dups3 <- sort(sim.dups[1:nsub])
sim.dels0 <- sort(sim.dels[c(1:nsub)+nsub]); sim.dups4 <- sort(sim.dups[c(1:nsub)+nsub])
snp.info <- read.snp.info(dir)

deL[["id2"]] <- deL$id; duP[["id2"]] <- duP$id
IDVEC <- deL$id; IDVEC2 <- duP$id
IDVEC2[sim.dups3] <- IDVEC2[sim.dups4] <- idz[[1]]
IDVEC[sim.dels1] <- IDVEC[sim.dels0] <- idz[[1]]
deL[["id"]] <- IDVEC; duP[["id"]] <- IDVEC2
new.sim.cnvs <- sim.cnvs <- rbind(deL[sim.dels1,],deL[sim.dels0,],duP[sim.dups3,],duP[sim.dups4,])
for (cc in 2:length(nn)) {
  sim.cnvs.next <- sim.cnvs
  sim.cnvs.next[["id"]] <- idz[[cc]]
  new.sim.cnvs <- rbind(new.sim.cnvs,sim.cnvs.next)
}
save(new.sim.cnvs,file=cat.path(sim.dir$ano,"simulatedCNVsFLAT.RData"))

# convert ranges into snp index in pfb file for each set
st.en.snps <- cbind(start.snp(snp.info,deL[sim.dels1,]),end.snp(snp.info,deL[sim.dels1,]))
st.en.snps[,1] <- match(st.en.snps[,1],rownames(pfb))
st.en.snps[,2] <- match(st.en.snps[,2],rownames(pfb))
sim.dels1 <- st.en.snps
st.en.snps <- cbind(start.snp(snp.info,deL[sim.dels0,]),end.snp(snp.info,deL[sim.dels0,]))
st.en.snps[,1] <- match(st.en.snps[,1],rownames(pfb))
st.en.snps[,2] <- match(st.en.snps[,2],rownames(pfb))
sim.dels0 <- st.en.snps
st.en.snps <- cbind(start.snp(snp.info,duP[sim.dups3,]),end.snp(snp.info,duP[sim.dups3,]))
st.en.snps[,1] <- match(st.en.snps[,1],rownames(pfb))
st.en.snps[,2] <- match(st.en.snps[,2],rownames(pfb))
sim.dups3 <- st.en.snps
st.en.snps <- cbind(start.snp(snp.info,duP[sim.dups4,]),end.snp(snp.info,duP[sim.dups4,]))
st.en.snps[,1] <- match(st.en.snps[,1],rownames(pfb))
st.en.snps[,2] <- match(st.en.snps[,2],rownames(pfb))
sim.dups4 <- st.en.snps






                    
plt <- rep(rownames(vv),length(nn))
plt.info <- cbind(colnames(LRR),plt,rep(0,ncol(LRR)))
colnames(plt.info) <- c("id","plate","phenotype")
write.table(plt.info,file=cat.path(sim.dir$ano,"plate.lookup.txt"),quote=F,row.names=F)

for (dd in 1:length(nn)) {
  baf <- lrr <- matrix(numeric(),nrow=nrow(pfb),ncol=nsub)
  for (cc in 1:nsub) {
    #make1sample:
    lrr[,cc] <- (sim.snp.lrr(nrow(pfb),"normal",noise=nn[dd],baf=pfb$PFB)) + vv[cc,select]
    baf[,cc] <- (sim.snp.baf(nrow(pfb),"normal",noise=nn[dd],baf=pfb$PFB))
    #insert sim.dels1
    pos <- c(sim.dels1[cc,1]:sim.dels1[cc,2])
    lrr[pos,cc] <- (sim.snp.lrr(length(pos),"del",noise=nn[dd],baf=pfb$PFB[pos])) + vv[cc,select[pos]]
    baf[pos,cc] <- (sim.snp.baf(length(pos),"del",noise=nn[dd],baf=pfb$PFB[pos]))
    #insertsim.dels0
    pos <- c(sim.dels0[cc,1]:sim.dels0[cc,2])
    lrr[pos,cc] <- (sim.snp.lrr(length(pos),"del0",noise=nn[dd],baf=pfb$PFB[pos])) + vv[cc,select[pos]]
    baf[pos,cc] <- (sim.snp.baf(length(pos),"del0",noise=nn[dd],baf=pfb$PFB[pos]))
    #insertsim.dups3
    pos <- c(sim.dups3[cc,1]:sim.dups3[cc,2])
    lrr[pos,cc] <- (sim.snp.lrr(length(pos),"dup",noise=nn[dd],baf=pfb$PFB[pos])) + vv[cc,select[pos]]
    baf[pos,cc] <- (sim.snp.baf(length(pos),"dup",noise=nn[dd],baf=pfb$PFB[pos]))
    #insertsim.dups4
    pos <- c(sim.dups4[cc,1]:sim.dups4[cc,2])
    lrr[pos,cc] <- (sim.snp.lrr(length(pos),"dup4",noise=nn[dd],baf=pfb$PFB[pos])) + vv[cc,select[pos]]
    baf[pos,cc] <- (sim.snp.baf(length(pos),"dup4",noise=nn[dd],baf=pfb$PFB[pos]))
  }
  LRR[,(((dd-1)*nsub)+c(1:nsub))] <- lrr
  BAF[,(((dd-1)*nsub)+c(1:nsub))] <- baf
  cat(".")
}

LRRBIG <- as.big.matrix(LRR,backingfile="LRRbckFile",descriptorfile="LRRdescrFile",backingpath=sim.dir$big)
BAFBIG <- as.big.matrix(BAF,backingfile="BAFbckFile",descriptorfile="BAFdescrFile",backingpath=sim.dir$big)

#dir <- make.dir("/chiswick/data/ncooper/simRunTest/")

auxdir <- "/home/ncooper/Documents/necessaryfilesSIM"
map.fn <- cat.path(auxdir,"snpdata.map")

chip.settings <- list(dir.raw=NULL,dir.base=sim.dir$base,snp.support=map.fn,
                      gsf=F,aux.files.dir=auxdir)
#chip.settings <- list(dir.raw=dir_rawM,dir.base=dir_baseM,snp.support=s.supM,
#                     gsf=gsfM,aux.files.dir=auxdirM)


base.settings <- list(dt.name="datatracker",
                      delete.as.we.go=F,
                      grps=c(1),snp.fields=NULL,geno.file=NULL,
                      big.lrr=paste("LRR","descrFile",sep=""),big.baf="BAFdescrFile",
                      run.mode=c("scratch","normal","big")[3],
                      snp.run.mode=c("normal","skip","plink")[2],
                      plink.imp=F,HD.mode=F,
                      n.cores=22,low.ram=F,
                      hide.penn.plink=T,penn.path="/usr/local/bin/penncnv64/",
                      ucsc="hg18",erase.previous=F,verbose=F)

snp.settings <- list(callrate.samp.thr=.94,callrate.snp.thr=.985,hwe.thr=0.0001,
                     snp.grp.miss=T,grp.hwe.z.thr=4,grp.cr.thr=.001,het.lo=.17,het.hi=.25)

samp.settings <- list(nSD=3,mean.thr=c("LB","UB"),dlrs.thr=c(NA,"UB"),gc.thr=c("LB","UB"),
                      badPlateThresh=0.50,skip.chr.ab=T,lob=2,hib=2.5,pctile.bound=0.01,cohort.pc.correct=F,
                      batch="plate",other.batch=list(),
                      lrr.report=T,chr.ab.report=T,plate.report=T)

pca.settings <- list(num.pcs=2,pc.to.keep=.20,assoc=F,n.store=50,correct.sex=F,
                     comparison=T,comp.gc=T,comps="plate")

penn.settings <- list(hmm="hh550.hmm",relative=T,run.manual=F,print.cmds=F,q.cores=50,
                      grid.id="all.q",cluster.fn="q.cmd")

cnv.settings <- list(out.format="Ranges",results="everything",print.summary.overlaps=T,
                     cnv.qc=T,rare.qc=F,plate.qc=T,pval=0.01,del.rate=0.4,dup.rate=0.2,thr.sd=2,plate.thr=3,
                     rmv.low.plates=F,min.sites=10,rare.olp=0.5,rare.pc=0.04,rmv.bad.reg=T)

settings <- c(chip.settings,base.settings,snp.settings,samp.settings,pca.settings,penn.settings,cnv.settings)

###################


dir <- make.dir("/chiswick/data/ncooper/simRunTest/")

cnvResult <- plumbcnv(settings,start.at=3,pause.after=6,restore.mode=T,result.pref="cnvResultsPCAFLAT2")


#FOR A SIMULATION:
# run on all sim, or real plate based data
# don't make cnvs too common or PCA gets screwed up (alternatively skip pca)
# evaluate by overlaps with original set
# evaluate QS's for true overlaps versus false overlaps; true should be high, false near 0
# graph/tabulate missed and try to ascertain why:
#  previously due to reasons including: overlapping dup/del; PCA artifact, HLA region removed

