source('/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R')
load.all.libs()

if(F) {
source('/chiswick/data/ncooper/ImmunochipReplication/Scripts/plotSampleCNVvsLRRBAF.R')
print(load("/chiswick/data/ncooper/immunochipRunTest/RESULTS/cnvResultsPCA2424.RData"))
cnvTest <- cnvResults[[4]]; cnvTest22 <- cnvTest[22]
dir <- make.dir("/chiswick/data/ncooper/immunochipRunTest/")
snp.info <- read.snp.info(dir)
sample.info <- read.sample.info(dir)

cnv.plot(dir="",samples="",LRR=T,BAF=F,PREPOSTPC=F,chr.expand=0,Chr=1:22,Pos=NA,Cnv=Pos,scl=10^6,
         tag.cnvs=F,medSmooth=F,med.rat=10,gcOverlay=F,geneOverlay=F,exons=F,
         bafOverlay=F,pfb.rng=NULL,hzOverlay=F,hz.rng=NULL,
         cnvPlotFileName="",show.fail=F,
         lrr.file="LRRFiltSortdescrFile",baf.file="BAFFiltSortdescrFile",
         snp.info=NULL,scheme=NA,col1="lightblue",col2="lightgreen",
         cnv.col="orange", cnv.lty="dashed", gc.col="blue",ov.lty="dotted",
         hz.col="navy",pfb.col="brown",gene.col="green",
         c.xlim=NULL,c.ylim=NULL,cust.sub="",ucsc="hg18",DT=NULL)



cnv.plot(dir=dir,samples="5299857114_R01C02",BAF=T,Chr=16,
         Pos=c(25,50)*10^6,Cnv=c(31239727,47042416),tag.cnvs=T,lrr.file="big.lrr")

cnv.plot(dir=dir,samples=jj,BAF=T,Chr=19,show.fail=T,
         Pos=c(55224893,55445233),Cnv=c(55327893,55342233),tag.cnvs=T,lrr.file="big.lrr")


#chr19:60007252-60069023   immunoKIR

ii <- reader("/home/ncooper/kirs.txt")
mt.ids <- reader("/home/ncooper/Documents/necessaryfilesICHIP/alt.id.lookup.txt")
jj <-narm(mt.ids[[1]][match(ii,mt.ids[[2]])])
cnv.plot(dir=dir,samples=jj,BAF=T,Chr=19,show.fail=T,
         Pos=c(59000000,60080000),Cnv=c(60007252,60069023),tag.cnvs=T,lrr.file="big.lrr")

# oo <- convert.textpos.to.data(rep("chr19:55327893-55342233",length(jj)))
# oo <- data.frame.to.ranged(as.data.frame(oo))
# oo[["id"]] <- jj
# plot.all.ranges(oo,DT=NULL,file="niko.ranges.pdf",dir="",pc.flank=0.5,snp.info=NULL,
#                 col1="black",scheme="mono",LRR=T,BAF=T,bafOverlay=F,hzOverlay=F,
#                 PREPOSTPC=F,show.fail=T,baf.file="big.baf",lrr.file="big.lrr") 



plot.all.ranges(cnv.ranges,DT=NULL,file="all.ranges.pdf",dir="",pc.flank=5,snp.info=NULL,
                col1="black",scheme="mono",LRR=T,BAF=F,bafOverlay=F,hzOverlay=F,
                PREPOSTPC=F,baf.file="big.baf",lrr.file="big.pcc",...) 


plot.all.ranges(cnvTest22[7:9,],file="all.ranges.pdf",dir=dir,pc.flank=5,snp.info=snp.info[,-2],
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=T,lrr.file="big.lrr",bafOverlay=T) 

plot.all.ranges(x.result[[4]],file="all.ranges.nonroh.gcon.pdf",dir=dir,pc.flank=5,
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=F,lrr.file="big.pcc",bafOverlay=F)

plot.all.ranges(DUP[1:50,],file="all.ranges.nonroh.gcoff.DUP.pdf",dir=dir,pc.flank=5,
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=F,lrr.file="big.pcc",bafOverlay=F)



plot.all.ranges(cnvResults[[5]],file="all.ranges.pdf",dir=dir,pc.flank=5,
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=T,lrr.file="big.lrr",bafOverlay=T)

send.to.pwf("all.ranges.pdf",dir$res)

plot.all.ranges(cnvTest[1][1:100,],file="all.ranges1_1_100pc.pdf",dir=dir,pc.flank=5,snp.info=snp.info[,-2],
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=F,lrr.file="big.pcc",geneOverlay=T)

#bafOverlay=T,gcOverlay=T
}





#dd <- as.numeric(snpRaw[1,]); extractROH(dd)


#dir <- make.dir("/chiswick/data/ncooper/immunochipRunTest/")

#DT <- read.data.tracker(dir)
#paste("chr",chr(dd[starts[1]]),start(dd[starts[1]),"-",end(dd[ends[1]]))
#rm(snpMat)

doROH <- F

if(doROH) {
  overhang.thr <- .20 # if cnv is only e.g, 20% of a region, be suspicious
  cnvResults <- get(paste(load(getSlot(DT,"cnv.result",n.pcs=6))))
  all.samp <- c(cnvResults$rareDEL$id,cnvResults$allDel$id)
  all.samp <- unique(all.samp)
  sml.fn <- getSlot(DT,"snp.stats")
  snpMatLst <- get(paste(load(sml.fn)))
  ROH.list <- lapply(snpMatLst,get.ROH.for.SnpMatrix,dir=dir,min.size=100,sample.list=all.samp)
  if(length(ROH.list)>1) { ROH.list <- do.call("c",args=ROH.list) }  
  save(ROH.list,file=cat.path(dir$ano,"ROH.list.RData"))
  NN <- nrow(cnvResults$rareDEL)
  flag <- rep(F,NN)
  for(cc in 1:NN) {
    loop.tracker(cc=cc,max=NN)
    if(cnvResults$rareDEL$cn[cc]==1) {
        #pc of ROH that CNV is
        subj <- cnvResults$rareDEL[cc,]; qur <- ROH.list[[paste(cnvResults$rareDEL$id[cc])]]
        pc.of.ROH <- overlap.pc(subject=subj,query=qur,name.by.gene=F,rel.query=T,fill.blanks=F)$pc
        #pc of CNV covered by ROH
        pc.cnv.cov <- overlap.pc(subject=cnvResults$rareDEL[cc,],query=ROH.list[[paste(cnvResults$rareDEL$id[cc])]],name.by.gene=F,rel.query=F,fill.blanks=F)$pc
        pc.of.ROH <- as.numeric(pc.of.ROH); pc.cnv.cov <- as.numeric(pc.cnv.cov)
        if(length(pc.cnv.cov)>0) {
          if(pc.cnv.cov==1) {
            if(length(pc.of.ROH)>0) {
              # set marker if CNV is just a small island within a big ROH region
              if(pc.of.ROH<overhang.thr) { flag[cc] <- T }
            }
          }
        }
    }
  }
  print(which(flag))
  save(flag,file="/chiswick/data/ncooper/immunochipRunTest/flag3.RData")
  cnvResults$rareDEL[["roh"]] <- flag
}




#dir <- make.dir("/chiswick/data/ncooper/immunochipRunTest/")
#cnv.plot(sample="5412557027_R02C02",Chr=1,Pos=c(159.7,159.9)*10^6
#         ,Cnv=c(159765328, 159776645),dir=dir,BAF=T,geneOverlay=T,tag.cnvs=T)



#ROH.list2 <- get.ROH.for.SnpMatrix(snpMatLst[[2]],dir=dir,min.size=100)

#save(ROH.list2,file=cat.path(dir$ano,"ROH.list.post2.RData"))

#ROH.list1 <- get.ROH.for.SnpMatrix(snpMatLst[[1]],dir=dir,min.size=100)
#save(ROH.list1,file=cat.path(dir$ano,"ROH.list.post1.RData"))
#ROH.list3 <- get.ROH.for.SnpMatrix(snpMatLst[[3]],dir=dir,min.size=100)
#save(ROH.list3,file=cat.path(dir$ano,"ROH.list.post3.RData"))



## deprecated
# get.flanks.from.big.mat.alt <- function(ranges,bigMat,ratio=5,bp=NA,snp.info=NULL,L=F,R=F,both=T) {
#   fl <- get.ratio.set(ranges,ratio=ratio,bp=bp)
#   idz <- (ranges$id)
#   if(!L & !R & !both) { warning("L, R and both all unselected, returning nothing"); return(NULL) }
#   if(L | both) {
#     left <- RangedData(ranges=IRanges(start=fl[,1],end=fl[,2]),space=chr(ranges))
#     flanking_1 <- range.snp(snp.info,toGenomeOrder(left,strict=T))
#   }
#   if(R | both) {
#     right <- RangedData(ranges=IRanges(start=fl[,3],end=fl[,4]),space=chr(ranges))
#     flanking_2 <- range.snp(snp.info,toGenomeOrder(right,strict=T))
#   }
#   if(nrow(flanking_1)!=length(idz)) { stop("unequal lengths, or id column missing from ranges") }
#   if(both) {
#     flanking <- lapply(list(flanking_1,flanking_2),big.extract.snp.ranges,samples=idz,bigMat=bigMat,snp.info=snp.info)
#     flankingB <- mapply(FUN=c,flanking[[1]],flanking[[2]])
#   } else { flankingB <- NULL }
#   if(L) {
#     if (both) { flankingL <- flanking[[1]] } else {
#       flankingL <- big.extract.snp.ranges(flanking_1,samples=idz,bigMat=bigMat,snp.info=snp.info)
#     }
#   } else { flankingL <- NULL }
#   if(R) {
#     if (both) { flankingR <- flanking[[2]] } else {
#       flankingR <- big.extract.snp.ranges(flanking_2,samples=idz,bigMat=bigMat,snp.info=snp.info)
#     }
#   } else { flankingR <- NULL }
#   outlist <- list(both=flankingB, left=flankingL, right=flankingR)
#   return(flanking)
# }


if(F) {
  cnt <- 2
  non.cnv <- c(bafDat2L[[cnt]],bafDat2R[[cnt]])
  one.cnv <- cnvbaf[[cnt]]
  all.dat <- c(one.cnv,non.cnv)
  labs <- c(rep(0,length(one.cnv)),rep(1,length(non.cnv)))
  ii <- glm(labs~all.dat, family = binomial(logit))
  mysumfun(ii)
}





runQS <- F

if(runQS) {
  dir <- make.dir("/chiswick/data/ncooper/simRunTest/") #make.dir("/chiswick/data/ncooper/immunochipRunTest/")
  DT <- read.data.tracker(dir)
  cnvset <- get(paste(load(getSlot(DT,"cnvresult"))))
  resort <- match(cnvset[[4]]$id,toGenomeOrder(cnvset[[4]])$id)
  
  DEL <- toGenomeOrder(cnvset[[4]])
  DUP <- toGenomeOrder(cnvset[[5]])
  QS.results <- get.quality.scores(DEL,dir)
  QS.results.matrix <- make.qs.table(QS.results)
  # quality score 'pass-rate' at a given threshold
  x <- QS.results.matrix[QS.results.matrix[[4]]>=10,1]; 
  print(length(which(x<=.95))/length(x))
  # quality score by snp count
  ii <- tapply(QS.results.matrix[[1]],factor(QS.results.matrix[[4]]),mean)
  print(ii)
  #DUP <- toGenomeOrder(cnvset[[5]])
  QS.results.matrix2 <- get.quality.scores(DUP,dir)
  
  # quality score 'pass-rate' at a given threshold
  x <- QS.results.matrix2[QS.results.matrix2[[4]]>=10,1]; 
  print(length(which(x<=.95))/length(x))  
  # quality score by snp count
  ii2 <- tapply(QS.results.matrix2[[1]],factor(QS.results.matrix2[[4]]),mean)
  print(ii2)
}
  
  
