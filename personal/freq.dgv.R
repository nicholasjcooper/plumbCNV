### this does the DGV validation analysis all in one hit
### filtering by rare % on the DGV, then optionally on QS thr, etc for CNV set.
## most of the functions slightly modified from the original dgv.validation script
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/validation.functions.R")


n.cores <- 1
mf <- .02 # maximum frequency of SNPs
dir <- make.dir("/chiswick/data/ncooper/immunochipRunTwo/")
from.scratch <- F
thr <- .5

exclude.keyword.from.dgv.col <- function(dgv,keywords,col="SampleSize",verbose=F,dir=NULL,inv=F,sel=F)
{
  badones <- rep(F,nrow(dgv))
  #keywords <- get.vec.multi.type(keywords,dir)
  if(length(keywords)<1) { warning("keywords list empty"); return(dgv) }
  if(!is.character(keywords)) { warning("keywords list not character()"); return(dgv) }
  if(all(keywords=="")) { return(dgv) }
  for (cc in 1:length(keywords)) {
    badz <- grep(paste(keywords[cc]),x=paste(dgv[[col]]),ignore.case=T)
    if(verbose) { cat(length(badz),"found with",keywords[cc],"\n") }
    badones[badz] <- T
  }
  if(inv) { badones <- !badones }
  if(sel) { return(badones) }
  if(verbose) { cat("Excluded ",length(which(badones)),"/",length(badones)," rows\n",sep="") }
  return(dgv[!badones,])
}


create.validation.set <- function(dir,ukw,n.cores=1,mode="cghseq",verbose=F) {
  dgv.methods <- c("MCD_analysis","MassSpec","Optical_mapping","Paired-end_mapping",
                   "OEA_assembly","Composite_approach","FISH","MLPA","PCR","Southern",
                   "qPCR","Read-depth_analysis","ROMA","BAC_aCGH","Oligo_aCGH",
                   "Sequence_alignment","Sequencing","SNP_array","SNP_genotyping_analysis")
  ## create a set of DGV cnv's that are likely valid and we have good chance of finding
  dgv <- dgv.subset.cov.by.n.snps(dir=dir,snp.range=1:20,max.freq=mf,
                                  do.plot=F,n.cores=n.cores,add.col=T)
  dgv2 <- dgv[dgv$chip.coverage>0,] # only parts of the DGV with some coverage by ichip markers
  
  #dgv2 <- exclude.keyword.from.dgv.col(dgv2,ukw,col="SampleSize",verbose=verbose) # remove non euro
  kw <- c("snp","sequenc","cgh")
  ## create logical selection vectors for studies with sequencing, CGH, bead and other
  meth.col <- c("Method.platform","method")
  meth.col <- meth.col[meth.col %in% colnames(dgv2)]
  SEQ <- exclude.keyword.from.dgv.col(dgv2,"sequenc",col=meth.col,verbose=verbose,inv=F,sel=T)
  CGH <- exclude.keyword.from.dgv.col(dgv2,"cgh",col=meth.col,verbose=verbose,inv=F,sel=T)
  BEAD <- exclude.keyword.from.dgv.col(dgv2,"snp",col=meth.col,verbose=verbose,inv=F,sel=T)
  OTHER <- exclude.keyword.from.dgv.col(dgv2,kw,col=meth.col,verbose=verbose,inv=T,sel=T)

  if(tolower(mode)=="cgh")  dgv3 <- dgv2[CGH,] # select only CGH datasets
  if(tolower(mode)=="cghseq")  dgv3 <- dgv2[CGH | SEQ,] # select only CGH and SEQUENCE datasets
  if(tolower(mode)=="bead")  dgv3 <- dgv2[BEAD,] # select only Beadchip datasets
  if(tolower(mode)=="all")  dgv3 <- dgv2 # select all relevant DGV datasets

  dgv3$Loss[is.na(dgv3$Loss)] <- 0  # set NA values for DELs/DUPs to zeros
  dgv3$Gain[is.na(dgv3$Gain)] <- 0
  
  # create DGV subsets for dels/dups separately across levels of n-SNP coverage
  dgv.dels.10 <- dgv3[dgv3$chip.coverage>9 & dgv3$Loss>0,]
  dgv.dups.10 <- dgv3[dgv3$chip.coverage>9 & dgv3$Gain>0,]
  dgv.dels.8 <- dgv3[dgv3$chip.coverage>7 & dgv3$Loss>0,]
  dgv.dups.8 <- dgv3[dgv3$chip.coverage>7& dgv3$Gain>0,]
  dgv.dels.6 <- dgv3[dgv3$chip.coverage>5 & dgv3$Loss>0,]
  dgv.dups.6 <- dgv3[dgv3$chip.coverage>5 & dgv3$Gain>0,]
  dgv.dels.4 <- dgv3[dgv3$chip.coverage>3 & dgv3$Loss>0,]
  dgv.dups.4 <- dgv3[dgv3$chip.coverage>3 & dgv3$Gain>0,]

  # merge into a list
  del.validation.set <- list(d10=dgv.dels.10,d8=dgv.dels.8,d6=dgv.dels.6,d4=dgv.dels.4)
  dup.validation.set <- list(d10=dgv.dups.10,d8=dgv.dups.8,d6=dgv.dups.6,d4=dgv.dups.4)
  dgv.validation.set <- list(DEL=del.validation.set,DUP=dup.validation.set)
  return(dgv.validation.set)
}


do.dgv.n.overlaps <- function(cnvResults,dir,comps=c(1:2),dbs=c(1:4),DVS,...) {
  # pass the results list from plumbCNV() and include any args of 'find.overlaps()' for filtering
  fail <- F
  if(!is.list(cnvResults)) { fail <- T }
  if(length(cnvResults)!=2) { fail <- T }
  if(!all(sapply(cnvResults,is)[1,]=="RangedData")) { fail <- T }
  if(fail) { 
    warning("cnvResults must be a list produced by plumbCNV() containing:\n",
            " 5 RangedData objects: allCNV, allDel, allDup, rareDel, rareDup")
    return(NULL) 
  }
  # initialise values for each run of the loop
  if(!all(comps %in% 1:2)) { comps <- 1:2 }; n.c <- length(comps)
  if(!all(dbs %in% 1:4)) { dbs <- 1:4 }; n.d <- length(dbs)
  DUPs <- c(F,T,F,T)[comps]
  DELs <- c(T,F,T,F)[comps]
  dupdel <- c("DEL","DUP","DEL","DUP")[comps]
  compz <- c("Rare DELs","Rare DUPs")[comps]
  set.n <- c(1:2)[comps]
  db <- c("d10","d8","d6","d4")[dbs]
  mns <- c(10,8,6,4)[dbs]
  lab <- c("DGV-10snp","DGV-8snp","DGV-6snp","DGV-4snp")[dbs]
  dbz <- vector("list",n.d); names(dbz) <- db
  # loop for finding overlaps and doing counts
  for (dd in 1:n.d) {
    # setup list substructures and names
    dbz[[dd]] <- vector("list",2); names(dbz[[dd]]) <- c("overlaps","counts")
    dbz[[dd]]$overlaps <- vector("list",n.c); names(dbz[[dd]]$overlaps) <- compz
    dbz[[dd]]$counts <- dbz[[dd]]$overlaps
    for (cc in 1:n.c) {
      tl <- paste("Running overlap analysis for:",compz[cc],"with",lab[dd])
      cat("\n",tl,"\n")  #,paste(rep("=",nchar(tl)),collapse=""),"\n\n",sep="")
      #require("proftools")
      #tf <- "Rproftemp.out"
      #Rprof("Rproftemp.out")
      
      dbz[[dd]]$overlaps[[cc]] <- oo <- (find.overlaps(cnvResults[[set.n[cc]]],
                                        ref=DVS[[dupdel[cc]]][[db[dd]]],min.sites=mns[dd],
                                        DEL=DELs[cc],DUP=DUPs[cc],dir=dir,fill.blanks=F,...))
      #Rprof()
      #rd <- readProfileData("Rproftemp.out")
      #tab <- flatProfile(rd, F)
      #print(tab)
      dbz[[dd]]$counts[[cc]] <- count.cnvs(oo[!is.na(as.numeric(names(oo)))],content.txt=lab[dd],
                       cnv.sample.map=cnvResults[[set.n[cc]]],by.phenotype=T,tot.cnvs=nrow(cnvResults[[set.n[cc]]]))
    }
  }
  return(dbz)
}






ukw.fn <- cat.path(dir$ano,"unwantedKeywords.txt") # list of keywords indicating non-european samples
ukw <- reader(ukw.fn)

sav.fn <- character(3)
sav.fn[1] <- cat.path(dir$ano,"dgv.validation.allCGHSEQ.RData")
sav.fn[2] <- cat.path(dir$ano,"dgv.validation.allBEAD.RData")
sav.fn[3] <- cat.path(dir$ano,"dgv.validation.allALL.RData")

if(from.scratch) {
  # don't exclude non europeans
  dgv.validation.set <- create.validation.set(dir,"",n.cores=n.cores,mode="cghseq",verbose=F)
  save(dgv.validation.set,file=sav.fn[1])
  dgv.validation.set <- create.validation.set(dir,"",n.cores=n.cores,mode="bead",verbose=F)
  save(dgv.validation.set,file=sav.fn[2])
  dgv.validation.set <- create.validation.set(dir,"",n.cores=n.cores,mode="all",verbose=F)
  save(dgv.validation.set,file=sav.fn[3])
}



####################
# START MAIN PART ##
validate <- T
# here, takes the combined set of all CNVs from all runs from 'parse.results.files'
# and then filters on QS
# then generates the full set of overlaps
# these overlaps can then be used to split results for each RUN

if(validate) {
  SCT.list <- vector("list",length(sav.fn)); names(SCT.list) <- basename(sav.fn)
  scrs <-  all.runs.del$score
  cat(length(scrs[scrs>=thr]),"/",length(scrs)," DELs pass ",thr," threshold\n",sep="")
  print(summary(scrs))
  all.del <- all.runs.del[scrs>=thr,]
  scrs2 <- all.runs.dup$score
  cat(length(scrs2[scrs2>=thr]),"/",length(scrs2)," DUPs pass ",thr," threshold\n",sep="")
  print(summary(scrs2))
  all.dup <- all.runs.dup[scrs>=thr,]
  # save QS to file
  del.denom <- rev(cumsum(rev(table(all.del$numSnps))))[c("4","6","8","10")]
  dup.denom <- rev(cumsum(rev(table(all.dup$numSnps))))[c("4","6","8","10")]
  x.result <- list(all.del,all.dup)
  names(x.result) <- c("Rare DELs","Rare DUPs")
  # do bulk of processing
  size.dgv <- list();
  count.per.cnv.del <- count.per.cnv.dup <- list.per.cnv.del <- list.per.cnv.dup <- list()
  for (cc in 1:length(sav.fn)) {
    dgv.validation.set <- get(paste(load(sav.fn[cc])))
    #c
    size.dgv[[cc]] <- cbind(sapply(dgv.validation.set[[1]],nrow),sapply(dgv.validation.set[[2]],nrow))
    ovlp.sum <- do.dgv.n.overlaps(x.result,dir=dir,DVS=dgv.validation.set,n.cores=n.cores,comps=c(1,2))
    count.per.cnv.del[[cc]] <- lapply(ovlp.sum,function(X) { X$counts[["Rare DELs"]]$overall.counts[[1]] } )
    count.per.cnv.dup[[cc]] <- lapply(ovlp.sum,function(X) { X$counts[["Rare DUPs"]]$overall.counts[[1]] } )
    list.per.cnv.del[[cc]] <- lapply(ovlp.sum,function(X) { X$overlaps[["Rare DELs"]] })
    list.per.cnv.dup[[cc]] <- lapply(ovlp.sum,function(X) { X$overlaps[["Rare DUPs"]] })
    SCT.list[[cc]] <- summary.counts.table(ovlp.sum,print.result=F)
    #  prt <- print(pheno.ratios.table(dir,sum.table=SCT.list[[cc]]))
    for (ee in 1:2)  { 
      prt <- pheno.ratios.table(dir,sum.table=SCT.list[[cc]])[[ee]]
      if("pheno.ratios" %in% names(prt)) { prt <- prt[["pheno.ratios"]][["1/0"]] 
        print(prt)
      } else {
        warning("phenotype ratio table for group ",cc,",",ee," could not be produced, perhaps not enough CNVs in some group(s)") 
      }
    }
  }
  
  ##### START HERE (after previous bit has run)
  # need to separate the sensitivity results by run
  # also rerun based on a lower frequency cutoff
  ########; dd <- ss <- 1
  
  ## this should extract the sensitivity for each run, listed
  ## by comparison source and number of snps
  db <- c(1,2,3,4); src <- c(1,2,3); ns <- c(10,8,6,4); srcs <- c("cghseq","bead","all")
  #count.per.cnv.del[[src]][[db]]
  #count.per.cnv.dup[[src]][[db]]
  res.list <- list(d10=list(),d8=list(),d6=list(),d4=list())
  res.list <- list(cghseq=res.list,bead=res.list,all=res.list)
  n.snp.del <- all.del$numSnps
  n.snp.dup <- all.dup$numSnps
  for(ss in src){
    for(dd in db) {
      ## get the full cnv list for a given source(cgh, etc) and given nsnp=10,8,6,4
      XL <- list.per.cnv.del[[src[ss]]][[db[dd]]] #DEL
      XP <- list.per.cnv.dup[[src[ss]]][[db[dd]]] #DUP
      # allow selections of only >nsnp
      right.size.l <- which(n.snp.del>=ns[db[dd]]) 
      right.size.p <- which(n.snp.dup>=ns[db[dd]])
      #head(rownames(all.del))
      cnv.nums <- names(XL); cnv.nums2 <- names(XP)
      # get just the ones with some DGV overlap and ensure select by
      #  right size (to make sense of original index)
      # RUN #'s
      xl.run.list <- factor(all.del$RUN[right.size.l[as.numeric(cnv.nums)]])
      xp.run.list <- factor(all.dup$RUN[right.size.p[as.numeric(cnv.nums2)]])
      # list of all cnvs in each run
      list.per.run.del <- Unlist(tapply(XL,xl.run.list,c),1)
      list.per.run.dup <- Unlist(tapply(XP,xp.run.list,c),1)
      # sensitivity, how many CNVs found in DGV (count unique DGV ids in each)
      unq.del <- sapply(list.per.run.del,function(X) { length(unique(X)) })
      unq.dup <- sapply(list.per.run.dup,function(X) { length(unique(X)) })
      pc.del <- unq.del/size.dgv[[src[ss]]][db[dd],1]
      pc.dup <- unq.dup/size.dgv[[src[ss]]][db[dd],2]
      # specificity, how many CNVs found in DGV (count total number in each run)
      n.del <- table(xl.run.list)
      n.dup <- table(xp.run.list)
      pc.del2 <- n.del/table(all.del$RUN[right.size.l])
      pc.dup2 <- n.dup/table(all.dup$RUN[right.size.p])
      res.list[[ss]][[dd]] <- cbind(unq.del,unq.dup,pc.del,pc.dup,n.del,n.dup,pc.del2,pc.dup2)
    }
  }
  
  ## setup the column indexes for source, nsnp and RUN
  for(cc in 1:3) { for(dd in 1:4) {
      df <- as.data.frame(res.list[[cc]][[dd]])
      df[["src"]] <- srcs[cc]
      df[["db"]] <- ns[dd]
      df[["run"]] <- rownames(df)
      res.list[[cc]][[dd]] <- df
    }
  } 
  res.list <- unlist(res.list,recursive=FALSE)
  res.list <- do.call("rbind",args=res.list)
  X <- res.list
#  X <- X[(as.numeric(X$run) %% 2) == 0,] # optionally select for odd/even runs; (cnv-qc)
  # OR use one prepared earlier, e.g, run on stats3, etc
  #save(res.list,file="/chiswick/data/ncooper/ImmunochipReplication/Scripts/myreslist.rda")
  #res.list <- reader("/chiswick/data/ncooper/ImmunochipReplication/Scripts/myreslist.rda")
  pdf("~/Documents/goodmaterials/overallSensSpecNsnpDbAll50.pdf")
  #SENSITIVTY VS SPECIFITY PLOTS with annotation, might need updated numbers
  #DEL
  plot(c(0.2,.65,0.2,.65),c(0.2,0.2,.8,.8),bty="l",
       xlab="sensitivity",ylab="specificity",main="rDEL",col="white")
  abline(a=0,b=1,lty="dashed",col="grey")
  abline(a=1,b=-1,lty="dashed",col="grey")
  points(X[,3],X[,7],pch=c(3,1,10)[match(X[,9],srcs)],
         col=c("red","purple","lightblue","blue")[match(X[,10],ns)])
  ##DUP
  plot(c(0.2,1,0.2,1),c(0,0,1,1),bty="l",
       xlab="sensitivity",ylab="specificity",main="rDUP",col="white")
  abline(a=0,b=1,lty="dashed",col="grey")
  abline(a=1,b=-1,lty="dashed",col="grey")
  points(X[,4],X[,8],pch=c(3,1,10)[match(X[,9],srcs)],
         col=c("red","purple","lightblue","blue")[match(X[,10],ns)])
  
  
  ## make a RINGs plot of the overall sensitivity and specificity against
  # all nsnps and CGHSEQ vs BEAD
  ## DUP
  library(car)
  plot(c(0.2,.8,0.2,.8),c(0.1,0.1,1,1),
       xlab="sensitivity",ylab="specificity",main="rDUP",col="white",bty="l")
  for(dd in 1:2) {
    for(cc in 1:4) {
      sel <- X[,9]==srcs[dd] & X[,10]==ns[cc] & X[,4]<.9
      dataEllipse(X[sel,4],X[sel,8], plot.points=F,levels=0.5,
                  col=c("red","purple","lightblue","blue")[cc],add=T,lty=c("solid","dotted")[dd])
    }
  }
  legend("topleft",legend=c("CGH/Sequence datasets","Bead chip datasets",
                             paste(c(10,8,6,4),"-SNP 50% confidence ellipse",sep="")),
         col=c("black","black","red","purple","lightblue","blue"),
         lty=c("solid","dotted","solid","solid","solid","solid"),bty="n",lwd=1.5)
  ## DEL
  plot(c(0.2,.6,0.2,.6),c(0.35,0.35,.8,.8),bty="l",
       xlab="sensitivity",ylab="specificity",main="rDEL",col="white")
  for(dd in 1:2) {
    for(cc in 1:4) {
      sel <- X[,9]==srcs[dd] & X[,10]==ns[cc]
      dataEllipse(X[sel,3],X[sel,7], plot.points=F,levels=0.5,
                  col=c("red","purple","lightblue","blue")[cc],add=T,lty=c("solid","dotted")[dd])
    }
  }
  legend("topleft",legend=c("CGH/Sequence datasets","Bead chip datasets",
                             paste(c(10,8,6,4),"-SNP 50% confidence ellipse",sep="")),
         col=c("black","black","red","purple","lightblue","blue"),
         lty=c("solid","dotted","solid","solid","solid","solid"),bty="n",lwd=1.5)
  
  ## DEL counts
  plot(c(0,1000,0,1000),c(0,1,0,1),bty="l",
       xlab="#rDELs",ylab="specificity",main="rDEL",col="white")
  for(dd in 1:2) {
    for(cc in 1:4) {
      sel <- X[,9]==srcs[dd] & X[,10]==ns[cc]
      dataEllipse(X[sel,1],X[sel,7], plot.points=F,levels=0.5,
                  col=c("red","purple","lightblue","blue")[cc],add=T,lty=c("solid","dotted")[dd])
    }
  }
  legend("topleft",legend=c("CGH/Sequence datasets","Bead chip datasets",
                            paste(c(10,8,6,4),"-SNP 50% confidence ellipse",sep="")),
         col=c("black","black","red","purple","lightblue","blue"),
         lty=c("solid","dotted","solid","solid","solid","solid"),bty="n",lwd=1.5)
  ## DUP counts
  library(car)
  plot(c(0,1000,0,1000),c(0.1,0.1,1,1),
       xlab="#rDUPs",ylab="specificity",main="rDUP",col="white",bty="l")
  for(dd in 1:2) {
    for(cc in 1:4) {
      sel <- X[,9]==srcs[dd] & X[,10]==ns[cc] & X[,4]<.9
      dataEllipse(X[sel,2],X[sel,8], plot.points=F,levels=0.5,
                  col=c("red","purple","lightblue","blue")[cc],add=T,lty=c("solid","dotted")[dd])
    }
  }
  legend("topleft",legend=c("CGH/Sequence datasets","Bead chip datasets",
                            paste(c(10,8,6,4),"-SNP 50% confidence ellipse",sep="")),
         col=c("black","black","red","purple","lightblue","blue"),
         lty=c("solid","dotted","solid","solid","solid","solid"),bty="n",lwd=1.5)
  dev.off()
  #regression results#
  summary(lm(X[,4][X[,9]==srcs[1]]~X[,10][X[,9]==srcs[1]])) # * DUP sensitivity versus nsnp [CGHSEQ]
  summary(lm(X[,3][X[,9]==srcs[1]]~X[,10][X[,9]==srcs[1]])) # *** DEL sensitivity versus nsnp [CGHSEQ]
  summary(lm(X[,8][X[,9]==srcs[1]]~X[,10][X[,9]==srcs[1]])) # .06 DUP specificity versus nsnp [CGHSEQ]
  summary(lm(X[,7][X[,9]==srcs[1]]~X[,10][X[,9]==srcs[1]])) # .91 DEL specificity versus nsnp [CGHSEQ]
  
  ##########################################################
  
  
  ######### FIGS FROM PARSE ###############
  # change to 1 or 2 to plot graphs for DELs/DUPs and run from here down to # ENDPLOT #
  dow <- c("DEL","DUP")[1]
  
  ### setup the vector of odds ratios for each run del,dup
  # vary '#' as to whether we look at the ratio considering quality score or not
  vec <- numeric()
  for (cc in as.numeric(names(mnz))) {
    if(dow=="DEL") {
      vec[cc] <- RUNS[[cc]][["divisors"]][7,10] # DELs with QS
    } else {
      vec[cc] <- RUNS[[cc]][["divisors"]][8,11] # DUPs with QS
    }
    #vec[cc] <- RUNS[[cc]][["divisors"]][5,8] # DELs without QS
    #vec[cc] <- RUNS[[cc]][["divisors"]][6,9] # DUPs without QS
  }
  
  lwdz <- c(0.75,1.00,2.00)
  cnv <- c("DEL","DUP")
  db <- c("CGHSEQ","BEAD","ALL")
  an <- c("SENS","SPEC")
  st <- c("COUNT","PC")
  sdbs <- c(10,8,6,4)
  
  pdf("~/Documents/goodmaterials/QS.OR50.pdf") #,height=7,width=5)
  par(mfrow=c(2,1))
  for(aa in 1:length(cnv)) {
    for(bb in 1:length(db)) {
      for(cc in 1:length(sdbs)) {
        for(dd in 1:length(an)) {
          # plot QS and OR side by side
          plot.dgv.result(qs=T,cnv=cnv[aa],an=an[dd],db=db[bb],nsnp=sdbs[cc])
          plot.dgv.result(qs=F,cnv=cnv[aa],an=an[dd],db=db[bb],nsnp=sdbs[cc])
        }
      }
    }
  }
  dev.off()
  
  ###############################
  ## this will do results for the overall - i.e, every freaking CNV we found
  # print key summary of each
  big.result <- as.data.frame(matrix(ncol=6,nrow=(length(sav.fn)*8)))
  colnames(big.result) <- c("measure","DELs","DUPs","DELsPc","DUPsPc","dB")
  for (cc in 1:length(sav.fn)) {
    SCT <- SCT.list[[cc]]
    #cat(sav.fn[cc])
    ij <- 3; ik <- c(3,4,5)
    TAB <- rbind(SCT[[1]][ij,ik],SCT[[2]][ij,ik],SCT[[3]][ij,ik],SCT[[4]][ij,ik])
    ij <- 2; ik <- c(3,4,5)
    TAB2 <- rbind(SCT[[1]][ij,ik],SCT[[2]][ij,ik],SCT[[3]][ij,ik],SCT[[4]][ij,ik])
    TAB <- as.data.frame(TAB); TAB2 <- as.data.frame(TAB2)
    TAB[[2]] <- as.numeric(TAB[[2]]); TAB2[[2]] <- as.numeric(TAB2[[2]])
    TAB[[3]] <- as.numeric(TAB[[3]]); TAB2[[3]] <- as.numeric(TAB2[[3]])
    big.result[((cc-1)*8)+1:4,1:3] <- (TAB)
    big.result[((cc-1)*8)+1:4,4:5] <- round(TAB[,2:3]/cbind(rev(del.denom),rev(dup.denom)),3)
    big.result[((cc-1)*8)+5:8,1:3] <- (TAB2)
    big.result[((cc-1)*8)+5:8,4:5] <- round(TAB2[,2:3]/size.dgv[[cc]],3)
    big.result[((cc-1)*8)+1:8,6] <- rmv.ext(basename(sav.fn[cc]))
  }
  big.result[[6]] <- gsub("dgv.validation.","",big.result[[6]])
  ofn <- cat.path(dir$res,"dgv.f",round(runif(1),2)*100,".validation.results",suf=suffix,ext="txt")
  write.table(big.result,file=ofn,quote=F)
  cat("wrote:",ofn,"\n")
}

save(vec, RUNS, runconfig, res, mnz, all.runs.del, all.runs.dup,
     lwdz, cnv, db, an, st, sdbs, srcs, res.list, src,
     ns, thr, mf, file="/chiswick/data/ncooper/ImmunochipReplication/Scripts/SaveThr00.RData")

#load("/chiswick/data/ncooper/ImmunochipReplication/Scripts/SaveThr00.RData")
#load("/chiswick/data/ncooper/ImmunochipReplication/Scripts/SaveThr50.RData")
#load("/chiswick/data/ncooper/ImmunochipReplication/Scripts/SaveThr90.RData")

## read in CNV indiv file
#cnv.ind <- "/chiswick/data/ncooper/immunochipRunTest/CNVQC/mytest.cnv.indiv"
#cnv.ind <- "/chiswick/data/ncooper/immunochipRunTest/CNVQC/rareDEL.cnv.indiv"

#myf <- read.table(cnv.ind,header=T)
#pheno.comp1 <- tapply(myf$KB,factor(myf$PHE),mean)
#pheno.comp2 <- tapply(myf$KBAVG,factor(myf$PHE),mean)
#print(pheno.comp1); print(pheno.comp2)

gnz <- c("LIME1","SLC2A4RG","PKIA","LRRC56","C11orf35","RASSF7","ARHGAP11B","MTMR15","MTMR10","TRPM1","KLF13","OTUD","XRCC6BP1")
#get.GO.for.genes(gnz)
gnz <- c("CHL1","CNTN6","RPL23AP38","CNTN4","IL5RA","TRNT1","CRBN","MPHOSPH9","C12orf65","CDK2AP1","SBNO1","SETD8","RILPL2","DENND1B","PLEKHG6","PFKP","PITRM1")
#get.GO.for.genes(gnz)


# oo <- extract.cnv.regions(dir,type="del",by.cnv=F,lwr=0.5,upr=1.5)
cbind(oo[[1]],oo[[2]],substr(oo[[3]],1,10))[order(oo[[1]]/oo[[2]]),]



fisher.test(x=matrix(c(4,6292,0,8332),nrow=2), y = NULL, workspace = 200000, hybrid = FALSE,
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

c(6,6292,1,8332)

mm <- matrix(c(5,6292,0,8332)
aylmer.test(mm)
       

# DELs Enriched in affected:
# POOL   IID PHE   CHR  BP1         BP2         KB     
# [1,] "S120" "7" "6:1" "20" "61838983"  "61844833"  "5.85"    LIME1; 'T cells' SLC2A4RG"  immune,  huntingtons                      
# [2,] "S129" "7" "6:1" "8"  "79647809"  "79670750"  "22.941"  "PKIA"                                    
# [3,] "S144" "6" "5:1" "11" "544299"    "561316"    "17.017"  "LRRC56;C11orf35;RASSF7 'rheum athritis, TNF'"                  
# [4,] "S167" "5" "4:1" "15" "28737821"  "30198807"  "1460.99" "ARHGAP11B;MTMR15;MTMR10;TRPM1 'melanoma' ;KLF13 'lymphocyte survival';OTUD"
# [5,] "S170" "5" "4:1" "12" "56619684"  "56623938"  "4.254"   "XRCC6BP1"     

# NB: none of these found in metabochip ctrls, also those more common in ctrls mostly non-genic, and genic also found in metabo
# (same can't be said for DUPs)
#6:1 p-value = 0.0005942 ; 5:1 p-value = 0.00366  ;  4:1 p-value = 0.01899

#DEL genes
#c("LIME1","SLC2A4RG","PKIA","LRRC56","C11orf35","RASSF7","ARHGAP11B","MTMR15","MTMR10","TRPM1","KLF13","OTUD","XRCC6BP1")
#DUP genes
#c("CHL1","CNTN6","RPL23AP38","CNTN4","IL5RA","TRNT1","CRBN","MPHOSPH9","C12orf65","CDK2AP1","SBNO1","SETD8","RILPL2","DENND1B","PLEKHG6","PFKP","PITRM1")

# DUPs Enriched in affected:
#  POOL   IID PHE   CHR  BP1         BP2         KB                                                  
#[1,] "S357" "6" "5:1" "3"  "259356"    "3174635"   "2915.28" "CHL1 'kidney cancer';CNTN6 [metabo];RPL23AP38;CNTN4 'autism';IL5RA;TRNT1;C"
#[2,] "S380" "5" "4:1" "12" "122170153" "122488359" "318.206" "MPHOSPH9 'multiple scl' ;C12orf65 'Encephalomyopathy' ;CDK2AP1;SBNO1;SETD8;RILPL2 'ciliary membrane content"
#[3,] "S416" "5" "5:0" "1"  "195836641" "195855698" "19.057"  "DENND1B"   - asthma                              
#[4,] "S440" "4" "4:0" "12" "6285172"   "6294322"   "9.15"    "PLEKHG6"                                 
#[5,] "S444" "4" "4:0" "10" "2001233"   "3196256"   "1195.02" "PFKP 'glycolosis, fructose';PITRM1 ' degrades mitochondrial peptides'                            
#[6,] "S462" "4" "4:0" "7"  "50247689"  "50257229"  "9.54"    ""        





# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.euroCGHSEQ.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit       570       629
# 2  CNVs with at least one DGV-8snp hit       672       724
# 3  CNVs with at least one DGV-6snp hit       908       979
# 4  CNVs with at least one DGV-4snp hit      1431      1356
# Rare DELs Rare DUPs
# 1 0.3147432 0.3307045
# 2 0.2839037 0.3213493
# 3 0.2915864 0.3447183
# 4 0.3006934 0.3205674
# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.euroBEAD.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit       472       836
# 2  CNVs with at least one DGV-8snp hit       627       950
# 3  CNVs with at least one DGV-6snp hit       829      1263
# 4  CNVs with at least one DGV-4snp hit      1555      2063
# Rare DELs Rare DUPs
# 1 0.2606295 0.4395373
# 2 0.2648923 0.4216600
# 3 0.2662171 0.4447183
# 4 0.3267493 0.4877069
# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.euroALL.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit       923      1086
# 2  CNVs with at least one DGV-8snp hit      1166      1260
# 3  CNVs with at least one DGV-6snp hit      1470      1634
# 4  CNVs with at least one DGV-4snp hit      2454      2553
# Rare DELs Rare DUPs
# 1 0.5096632 0.5709779
# 2 0.4926067 0.5592543
# 3 0.4720617 0.5753521
# 4 0.5156545 0.6035461
# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.allCGHSEQ.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit       830       643
# 2  CNVs with at least one DGV-8snp hit       998       747
# 3  CNVs with at least one DGV-6snp hit      1121      1020
# 4  CNVs with at least one DGV-4snp hit      1684      1481
# Rare DELs Rare DUPs
# 1 0.4583103 0.3380652
# 2 0.4216308 0.3315579
# 3 0.3599872 0.3591549
# 4 0.3538559 0.3501182
# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.allBEAD.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit       499       851
# 2  CNVs with at least one DGV-8snp hit       664       967
# 3  CNVs with at least one DGV-6snp hit       880      1282
# 4  CNVs with at least one DGV-4snp hit      1641      2085
# Rare DELs Rare DUPs
# 1 0.2755384 0.4474238
# 2 0.2805239 0.4292055
# 3 0.2825947 0.4514085
# 4 0.3448203 0.4929078
# /chiswick/data/ncooper/immunochipRunTest/ANNOTATION/dgv.validation.allALL.RData                                 Count Rare DELs Rare DUPs
# 1 CNVs with at least one DGV-10snp hit      1061      1130
# 2  CNVs with at least one DGV-8snp hit      1353      1313
# 3  CNVs with at least one DGV-6snp hit      1679      1695
# 4  CNVs with at least one DGV-4snp hit      2706      2634
# Rare DELs Rare DUPs
# 1 0.5858642 0.5941115
# 2 0.5716096 0.5827785
# 3 0.5391779 0.5968310
# 4 0.5686069 0.6226950
