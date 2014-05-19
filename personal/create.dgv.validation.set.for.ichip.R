n.cores <- 1
mf <- NA # maximum frequency of SNPs

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


do.dgv.n.overlaps <- function(cnvResults,dir,comps=c(1:4),dbs=c(1:4),DVS,...) {
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
  if(!all(comps %in% 1:4)) { comps <- 1:4 }; n.c <- length(comps)
  if(!all(dbs %in% 1:4)) { dbs <- 1:4 }; n.d <- length(dbs)
  DUPs <- c(F,T,F,T)[comps]
  DELs <- c(T,F,T,F)[comps]
  dupdel <- c("DEL","DUP","DEL","DUP")[comps]
  compz <- c("All Deletions","All Duplications","Rare DELs","Rare DUPs")[comps]
  set.n <- c(2:5)[comps]
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



if(validate) {
  
  # a lot of this now duplicated in process.quality.scores() in validation.functions.R
  SCT.list <- vector("list",length(sav.fn)); names(SCT.list) <- basename(sav.fn)
  res.fn <- getSlot(DT,"cnvresult",n.pcs=pca.set)
  x.result <- get(paste(load(cat.path(dir$res,res.fn))))  #cnvResultsPCAAssoc9.RData"))))
  #x.result$rareDEL <- x.result$rareDEL[!x.result$rareDEL$roh,]
  ofn <- cat.path(dir$res,"qs.del.backup",suf=suffix,ext="RData")
  if(file.exists(ofn)) { print(load(ofn)) } else {
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
  cnvResults$rareDEL <- x.result$rareDEL; 
  cnvResults$rareDEL[["score"]] <- scrs;
  save(cnvResults,file=res.fn)
  x.result$rareDEL <- x.result$rareDEL[scrs>=thr,]
  qs.sc2 <- get.quality.scores(x.result$rareDUP,dir,n.pcs=pca.set)
  ofn <- cat.path(dir$res,"qs.dup.backup",suf=suffix,ext="RData")
  if(file.exists(ofn)) { print(load(ofn)) } else {
    qs.sc2 <- get.quality.scores(x.result$rareDUP,dir,n.pcs=pca.set)
    save(qs.sc2,file=ofn)
  }
  qs.mat2 <- make.qs.table(qs.sc2)
  ofn <- cat.path(dir$res,"qs.dup.results",suf=suffix,ext="txt")
  write.table(qs.mat2,file=ofn,quote=F);   cat("wrote:",ofn,"\n")
  scrs2 <- qs.mat2[[2]]; scrs2[is.na(scrs2)] <- qs.mat2[[3]][is.na(scrs2)]
  cat(length(scrs2[scrs2>=.7]),"/",length(scrs2)," DUPs pass .7 threshold\n",sep="")
  print(summary(scrs2))
  # save QS to file
  cnvResults$rareDUP <- x.result$rareDUP; 
  cnvResults$rareDUP[["score"]] <- scrs2;
  save(cnvResults,file=res.fn)
  x.result$rareDUP <- x.result$rareDUP[scrs2>=thr,]
  del.denom <- rev(cumsum(rev(table(x.result[[4]]$numSnps))))[c("4","6","8","10")]
  dup.denom <- rev(cumsum(rev(table(x.result[[5]]$numSnps))))[c("4","6","8","10")]
  # do bulk of processing
  size.dgv <- list(); 
  for (cc in 1:length(sav.fn)) {
    dgv.validation.set <- get(paste(load(sav.fn[cc])))
    #c
    size.dgv[[cc]] <- cbind(sapply(dgv.validation.set[[1]],nrow),sapply(dgv.validation.set[[1]],nrow))
    ovlp.sum <- do.dgv.n.overlaps(x.result,dir=dir,DVS=dgv.validation.set,n.cores=n.cores)
    SCT.list[[cc]] <- summary.counts.table(ovlp.sum,print.result=F)
    #  prt <- print(pheno.ratios.table(dir,sum.table=SCT.list[[cc]]))
    for (ee in 1:4)  { 
      prt <- pheno.ratios.table(dir,sum.table=SCT.list[[cc]])[[ee]]
      if("pheno.ratios" %in% names(prt)) { prt <- prt[["pheno.ratios"]][["1/0"]] 
        print(prt)
      } else {
        warning("phenotype ratio table for group ",cc,",",ee," could not be produced, perhaps not enough CNVs in some group(s)") 
      }
    }
  }
  # print key summary of each
  big.result <- as.data.frame(matrix(ncol=6,nrow=(length(sav.fn)*8)))
  colnames(big.result) <- c("measure","DELs","DUPs","DELsPc","DUPsPc","dB")
  for (cc in 1:length(sav.fn)) {
    SCT <- SCT.list[[cc]]
    #cat(sav.fn[cc])
    ij <- 3; ik <- c(3,6,7)
    TAB <- rbind(SCT[[1]][ij,ik],SCT[[2]][ij,ik],SCT[[3]][ij,ik],SCT[[4]][ij,ik])
    ij <- 2; ik <- c(3,6,7)
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
  ofn <- cat.path(dir$res,"dgv.validation.results",suf=suffix,ext="txt")
  write.table(big.result,file=ofn,quote=F)
  cat("wrote:",ofn,"\n")
}




## read in CNV indiv file
#cnv.ind <- "/chiswick/data/ncooper/immunochipRunTest/CNVQC/mytest.cnv.indiv"
#cnv.ind <- "/chiswick/data/ncooper/immunochipRunTest/CNVQC/rareDEL.cnv.indiv"

#myf <- read.table(cnv.ind,header=T)
#pheno.comp1 <- tapply(myf$KB,factor(myf$PHE),mean)
#pheno.comp2 <- tapply(myf$KBAVG,factor(myf$PHE),mean)
#print(pheno.comp1); print(pheno.comp2)

get.GO.for.genes <- function(gene.list,bio=T,cel=F,mol=F) {
  must.use.package(c("biomaRt","genoset","gage"),T)
  ens <- useMart("ENSEMBL_MART_ENSEMBL",
               dataset="hsapiens_gene_ensembl",
               host="may2009.archive.ensembl.org",
               path="/biomart/martservice",
               archive=FALSE)
  ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  data(egSymb)
  base.attr <- c("hgnc_symbol", "chromosome_name")
  if(bio) { base.attr <- c(base.attr,"go_biological_process_description") }
  if(cel) { base.attr <- c(base.attr,"go_cellular_component_description") }
  if(mol) { base.attr <- c(base.attr,"go_molecular_function_description") }
  dat <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
                            "start_position", "end_position", "band"), filters = "hgnc_symbol",
               values = egSymb[,2], mart = ens)
  results <- getBM(attributes = base.attr, filters = "hgnc_symbol",
                   values = c(gene.list), mart = ens)
  return(results)
}

gnz <- c("LIME1","SLC2A4RG","PKIA","LRRC56","C11orf35","RASSF7","ARHGAP11B","MTMR15","MTMR10","TRPM1","KLF13","OTUD","XRCC6BP1")
#get.GO.for.genes(gnz)
gnz <- c("CHL1","CNTN6","RPL23AP38","CNTN4","IL5RA","TRNT1","CRBN","MPHOSPH9","C12orf65","CDK2AP1","SBNO1","SETD8","RILPL2","DENND1B","PLEKHG6","PFKP","PITRM1")
#get.GO.for.genes(gnz)

extract.cnv.regions <- function(dir, type="DEL", by.cnv=F, enriched=T, genes=T) {
  #probably only work when 2 phenotypes are present
  dir <- validate.dir.for(dir,"cnv.qc")
  must.use.package("genoset",T)
  type <- toupper(type); if(type!="DUP") { type <- "DEL" }
  cnv.ov <- cat.path(dir$cnv.qc,pref="rare",fn=type,ext="cnv.overlap")
  ## read in cnv overlap file
  tt <- reader(cnv.ov)

  inter <- tt[tt[,2]=="CON",]
  un <- tt[tt[,2]=="UNION",] # CNVRs
  pp <- strsplit(un[,4],":",fixed=T)
  case <- as.numeric(sapply(pp,"[",1)); ctrl <- as.numeric(sapply(pp,"[",2))
  ii <- ((case/ctrl)[case!=0 & ctrl!=0]) #[1:100]
  cat("Ratio summary:\n"); print(summary(ii))
  case0 <- case; case0[case0==0] <- 0.5; ctrl0 <- ctrl; ctrl0[ctrl0==0] <- 0.5;
  ii0 <- (case0/ctrl0)
  main <- tt[tt[,2]!="UNION" & tt[,2]!="CON",]  # results for individuals

  CNVRs <- RangedData(IRanges(start=as.numeric(un[,6]),end=as.numeric(un[,7]),
                              names=un[,1]),space=as.numeric(un[,5]),cases=case,ctrls=ctrl)
  if(enriched) {
    many.ctrl <- un[which(ii0<=.25),]
    many.case <- un[which(ii0>=4),]
    if(genes) {
      UNAFFECTED.REG <- many.ctrl[,"POOL"]; indx2 <- match(UNAFFECTED.REG,rownames(CNVRs))
      AFFECTED.REG <- many.case[,"POOL"]; indx <- match(AFFECTED.REG,rownames(CNVRs))
      CNVRs[["genes"]] <- (find.overlaps(CNVRs,db="gene",vec.out=T,delim=";"))
      many.ctrl <- cbind(many.ctrl,substr(CNVRs$genes[indx2],1,40))[,-match(c("FID","TYPE","SCORE"),colnames(many.ctrl))]
      many.case <- cbind(many.case,substr(CNVRs$genes[indx],1,40))[,-match(c("FID","TYPE","SCORE"),colnames(many.case))]
    } 
    cat("\nEnriched in unaffected summary:\n"); print(many.ctrl); 
    cat("\nEnriched in affected:\n"); print(many.case) # examples where case/ctrl outnumber each other >=4:1
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

# oo <- extract.cnv.regions(dir,type="del",by.cnv=F)
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
