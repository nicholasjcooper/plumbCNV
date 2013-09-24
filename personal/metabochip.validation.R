
DT <- read.data.tracker(dir)
DT2 <- read.data.tracker(dir2)

gs <- get.gene.annot(dir)


cat("load metabochip results\n")
cnvResults <- get(paste(load(getSlot(DT2,"cnv.result",n.pcs=n.pcs2))))
metabo <- toGenomeOrder(cnvResults[[4]]) ; metabo <- metabo[metabo$numSnps>=min.snps,]; metabo <- annot.cnv(metabo,gs)
metabo2 <- toGenomeOrder(cnvResults[[5]]); metabo2 <- metabo2[metabo2$numSnps>=min.snps,]; metabo2 <- annot.cnv(metabo2,gs)

cat("load immunochip results\n")
cnvResults <- get(paste(load(getSlot(DT,"cnv.result",n.pcs=n.pcs))))
#print(load("/chiswick/data/ncooper/immunochipRunTest/RESULTS/cnvResultsPCA66.RData"))
#print(load("/chiswick/data/ncooper/immunochipRunTest/RESULTS/cnvResultsPCAAssoc9.RData"))
immuno <- toGenomeOrder(cnvResults[[4]]); immuno <- immuno[immuno$numSnps>=min.snps,]; immuno <- annot.cnv(immuno,gs)
immuno2 <-  toGenomeOrder(cnvResults[[5]]); immuno2 <- immuno2[immuno2$numSnps>=min.snps,]; immuno2 <- annot.cnv(immuno2,gs)

cat("add immunochip sample ids to metabochip dataset\n")
mt.ids <- reader(cat.path(dir$ano,"metabochipsamples.txt"))
metabo[["immuno"]] <- mt.ids[,1][match(metabo$id,mt.ids[,2])]
metabo2[["immuno"]] <- mt.ids[,1][match(metabo2$id,mt.ids[,2])]

#match("130205_C12_1958_BC1040955",mt.ids[[2]])


if(print.common) {
  length(which(immuno$id %in% metabo$immuno))
  im.list <- metabo[(metabo$immuno %in% immuno$id),"immuno"][[1]]
  ui <- get.all.samp.fails(dir)
  cat("[note,",length(which(im.list %in% ui)),"ids were failed samples]\n")
  
  for (cc in 1:length(im.list)) { 
    print(metabo[which(metabo$immuno==im.list[cc]),]); print(immuno[immuno$id==im.list[cc],]);
  }
}
  

ofn <- cat.path(dir$res,"metaboverlap",suf=suffix,ext="txt")
sink(ofn)
cat("find number of overlapping CNVs between immunochip and metabochip\n")
yy <- find.overlaps(immuno,ref=metabo)
cat("number of immunochip DELs that overlap with metabochip",length(yy[yy!=""]),"\n")
yyy <- yy[yy!=""]; uuu <- unlist(yyy)
cat("number of metabochip DEL hits:",length(uuu),"; unique hits:",length(unique(uuu)),"\n")
yy2 <- find.overlaps(immuno2,ref=metabo2)
cat("number of immunochip DUPs that overlap with metabochip",length(yy2[yy2!=""]),"\n")
yyy2 <- yy2[yy2!=""]; uuu2 <- unlist(yyy2)
cat("number of metabochip DUP hits:",length(uuu2),"; unique hits:",length(unique(uuu2)),"\n")
sink()

cat("find possible number of overlapping CNVs between immunochip and metabochip but seeing what snps are there\n")
snp.info <- read.snp.info(dir); snp.info2 <- read.snp.info(dir2)
iii <- find.overlaps(snp.info,ref=snp.info2)
cat("number of metabochip snps that are on immunochip",length(iii[iii!=""]),"\n")

metz <- dgv.subset.cov.by.n.snps(dir=dir,snp.set=snp.info,dgv=metabo,add.col=T,db="metabochip")
metz2 <- dgv.subset.cov.by.n.snps(dir=dir,snp.set=snp.info,dgv=metabo2,add.col=T,db="metabochip")
metz[["immuno"]] <- mt.ids[,1][match(metz$id,mt.ids[,2])]
metz2[["immuno"]] <- mt.ids[,1][match(metz2$id,mt.ids[,2])]


lkly <- 5 # number of snps intersecting that constitutes the claim 'likely'

n1 <- (length(which(paste(1:nrow(metz))[metz$chip.coverage>0] %in% paste(unique(uuu)))))
n2 <-(length(which(paste(1:nrow(metz))[metz$chip.coverage>lkly] %in% paste(unique(uuu)))))
n3 <- (length(which(paste(1:nrow(metz2))[metz2$chip.coverage>0] %in% paste(unique(uuu2)))))
n4 <- (length(which(paste(1:nrow(metz2))[metz2$chip.coverage>lkly] %in% paste(unique(uuu2)))))

sink(ofn,append=T)
cat("number of metabochip DELs that could be detectable on immunochip",nrow(metz[metz$chip.coverage>0,]),"; found:",n1,"\n")
cat("number of metabochip DELs that are likely detectable on immunochip",nrow(metz[metz$chip.coverage>lkly,]),"; found:",n2,"\n")
cat("number of metabochip DUPs that could be detectable on immunochip",nrow(metz2[metz2$chip.coverage>0,]),"; found:",n3,"\n")
cat("number of metabochip DUPS that are likely detectable on immunochip",nrow(metz2[metz2$chip.coverage>lkly,]),"; found:",n4,"\n")
sink()

dup.mode <- F
if(dup.mode) {
  #swap metz and metz2 around so extra processing is done on the DUPs
  metz3 <- metz; metz <- metz2; metz2 <- metz3
  thr <- .95 #max(get.pctile(QS.results.matrix[[1]],pc.thr))
  # all dups in common were validated, but some clearer than others
  gd.id <- reader("/chiswick/data/ncooper/immunochipRunTest/RESULTS/CLEAR_METABO.txt")
  bd.id <- reader("/chiswick/data/ncooper/immunochipRunTest/RESULTS/ROUGH_METABO.txt")
} else {
  thr <- .95 #max(get.pctile(QS.results.matrix[[1]],pc.thr))
  # only about half the common dels looked real
  gd.id <- reader("/chiswick/data/ncooper/immunochipRunTest/RESULTS/GOOD_METABO.txt")
  bd.id <- reader("/chiswick/data/ncooper/immunochipRunTest/RESULTS/BAD_METABO.txt")
}

if(do.qs) {
  mtab <- metz[metz$chip.coverage>lkly,]
  mtab <- (toGenomeOrder(mtab,strict=T))
  
  # use metabo ranges, but change IDs to immuno equivalents, then extraction
  # will use the immuno big.matrices in the same positions
  itab <- toGenomeOrder(mtab)  #[good.or.bad==1,])
  itab[["metabo"]] <- itab[["id"]]
  itab[["id"]] <- itab[["immuno"]]
  itab <- itab[!is.na(itab$id),]
  
  # quality scores for same locations in immunochip dataset
  QS.results2 <- get.quality.scores(itab,dir,n.pcs=n.pcs)
  QS.results.matrix2 <- make.qs.table(QS.results2)
  
  good.or.bad <- rep(NA,times=nrow(mtab))
  good.or.bad[mtab$id %in% gd.id] <- 1
  good.or.bad[mtab$id %in% bd.id] <- 0
  good.or.badI <- rep(NA,times=nrow(itab))
  good.or.badI[itab$metabo %in% gd.id] <- 1
  good.or.badI[itab$metabo %in% bd.id] <- 0
  # bonus or exclusion from badness for >40 SNPs?
  immuno.result <- cbind(itab$id,good.or.badI,QS.results.matrix2)
  immuno.result[[1]]

  # most wrong dups just seem like inexact window, but there still is a CNV
  res.mat2 <- (cbind(itab$id,good.or.badI,QS.results.matrix2))
  
  ofn <- cat.path(dir$res,"qs.immuno.results",suf=suffix,ext="txt")
  write.table(immuno.result,file=ofn,quote=F)
  cat("wrote:",ofn,"\n")
  
  if(do.metabo) {
    # original quality scores in the metabochip dataset
    QS.results <- get.quality.scores(mtab,dir2,n.pcs=n.pcs2) #,nsnp=31)
    QS.results.matrix <- make.qs.table(QS.results)
    metabo.result <- cbind(mtab$id,good.or.bad,QS.results.matrix)
    ofn <- cat.path(dir$res,"qs.metabo.results",suf=suffix,ext="txt")
    write.table(metabo.result,file=ofn,quote=F)
    cat("wrote:",ofn,"\n")
    res.mat <- (cbind(mtab$id,good.or.bad,QS.results.matrix))
    res.mat[order(res.mat[[2]]),]
    FN <- res.mat[which(res.mat[[2]]==1 & res.mat[[3]]<thr),]
    FP <- res.mat[which(res.mat[[2]]==0 & res.mat[[3]]>thr),]
 #   pdf(cat.path(dir$res,"plotQSvsNsnp_by_valid.pdf"))
    plot(res.mat$cnvqs,res.mat$n.snps,col=res.mat$good.or.bad+1,log="y",xlim=c(0,1))
    dev.off()
    
    must.use.package("Epi")
    pc.thr <- do.roc(y=good.or.bad,x=QS.results.matrix$cnvqs)
    
    grp <- QS.results.matrix$cnvqs>thr
    print(table(grp,good.or.bad))
    
    fps <- which(mtab$id %in% FP[[1]])
    fns <- which(mtab$id %in% FN[[1]])
    ##good.ones <- which(QS.results.matrix[[1]] > .95)
    mtab1 <- mtab[fps,] # 13/15, ,15.5/72 bad using QS
    mtab0 <- mtab[fns,] # 15/17 good using QS
  }
}


explore.all <- F
## EXPLORATORY PLOTS OF ALL ##
if(explore.all) {
  plotabo <- mtab

  #olp <- find.overlaps(immuno,ref=mtab[good.ones,])
  plot.range <- 1:nrow(plotabo)
  # use metabo ranges, but change IDs to immuno equivalents
  itabo <- plotabo
  itabo[["id"]] <- itabo[["immuno"]]
  itabo <- itabo[!is.na(itabo$id),]
  plot.all.ranges(itabo,file="immunoDUPs_in_metabo_spots.pdf",dir=dir,pc.flank=5,
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=F,lrr.file="big.pcc",bafOverlay=F)

  plot.all.ranges(plotabo,file="metaboDUPs_in_metabo_spots.pdf",dir=dir2,pc.flank=5,
                col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                show.fail=F,lrr.file="big.pcc",bafOverlay=F)

}
####
if(do.plots & do.qs) {
  plotabo <- mtab0

  #olp <- find.overlaps(immuno,ref=mtab[good.ones,])
  plot.range <- 1:nrow(plotabo)
  # use metabo ranges, but change IDs to immuno equivalents
  itabo <- plotabo
  itabo[["id"]] <- itabo[["immuno"]]
  itabo <- itabo[!is.na(itabo$id),]
  plot.all.ranges(itabo,file="FNimmunoDELs_in_metabo_spots.pdf",dir=dir,pc.flank=5,
                  col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                  show.fail=F,lrr.file="big.pcc",bafOverlay=F)
  
  plot.all.ranges(plotabo,file="FNmetaboDELs_in_metabo_spots.pdf",dir=dir2,pc.flank=5,
                  col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                  show.fail=F,lrr.file="big.pcc",bafOverlay=F)
  
  plotabo <- mtab1
  
  #olp <- find.overlaps(immuno,ref=mtab[good.ones,])
  plot.range <- 1:nrow(plotabo)
  # use metabo ranges, but change IDs to immuno equivalents
  itabo <- plotabo
  itabo[["id"]] <- itabo[["immuno"]]
  itabo <- itabo[!is.na(itabo$id),]
  plot.all.ranges(itabo,file="FPimmunoDELs_in_metabo_spots.pdf",dir=dir,pc.flank=5,
                  col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                  show.fail=F,lrr.file="big.pcc",bafOverlay=F)
  
  plot.all.ranges(plotabo,file="FPmetaboDELs_in_metabo_spots.pdf",dir=dir2,pc.flank=5,
                  col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=F,tag.cnv=T,
                  show.fail=F,lrr.file="big.pcc",bafOverlay=F)
  
#  plot.all.ranges(mtab[good.or.bad==0,][1:10,],file="whybiased.pdf",dir=dir2,pc.flank=5,
#                  col1="black",scheme="mono",LRR=T,BAF=T,PREPOSTPC=T,tag.cnv=T,
#                  show.fail=T,lrr.file="big.lrr")
  
}

# testing CNV set for overlaps with custom ranges ...
# number of metabochip DELs that could be detectable on immunochip 236 ; found: 83 
# number of metabochip DELs that are likely detectable on immunochip 142 ; found: 79 
# number of metabochip DUPs that could be detectable on immunochip 146 ; found: 68 
# number of metabochip DUPS that are likely detectable on immunochip 71 ; found: 46 

if(ind.examples) {
  #DUP genes
  dgg1 <- c("CHL1","CNTN6","RPL23AP38","CNTN4","IL5RA","TRNT1","CRBN","MPHOSPH9","C12orf65","CDK2AP1","SBNO1","SETD8","RILPL2","DENND1B","PLEKHG6","PFKP","PITRM1")
  #Dups T1D have but no controls have: (but both similar freq in metabochip - 15ish)
  #FCGR3B ZNF804B 
  #12      12 
  metabo2[grep("FCGR3B",metabo2$gene),]
  metabo2[grep("ZNF804B ",metabo2$gene),]
  for (dd in 1:length(dgg1)) {
    print(metabo2[grep(dgg1[dd],metabo2$gene),])
  }
  
  #DEL genes
  dgg2 <- c("LIME1","SLC2A4RG","PKIA","LRRC56","C11orf35","RASSF7","ARHGAP11B","MTMR15","MTMR10","TRPM1","KLF13","OTUD","XRCC6BP1")
  
    
  # Dels T1D have but no controls have:
  #   GJB6   PKIA       SULT1A1 , SULT1A2
  # 12     9 vs 2     7 vs  3
  metabo[grep("GJB6",metabo$gene),]
  metabo[grep("PKIA",metabo$gene),]
  metabo[grep("SULT1A",metabo$gene),]
  for (dd in 1:length(dgg2)) {
    print(metabo[grep(dgg2[dd],metabo$gene),])
  }
  
}
