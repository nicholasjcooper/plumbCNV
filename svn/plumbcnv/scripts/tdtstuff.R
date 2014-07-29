#play with TDT

if(!exists("tdt3")) {
  # only bother recalculating if these vars not already present
  print(load("/chiswick/data/ncooper/ImmunochipFamilies/RESULTS/TDT_results.RData"))
  dir <- make.dir("/chiswick/data/ncooper/ImmunochipFamilies")
  tdt3 <- add.all.ids(tdt3,ped,dir)
  tdt4 <- add.all.ids(tdt4,ped,dir)
}
# DN --\
# ---NN, NN, DN, DN
# NN --/ 
#   
#   DN --\
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


trio.analysis <- function(dir=NULL, cnvResults, ped.file) {
  dir <- validate.dir.for(dir,c("res","cnv.qc"))
  if(!is.list(cnvResults) | length(cnvResults)!=5 | !is(cnvResults[[1]])[1]=="RangedData") { 
    warning("cnvResults object should be a list of RangedData objects, length 5"); return(NULL) }
  cat("\nRunning TDT analysis for family trios\n")
  oo1 <- extract.cnv.regions(dir,type="del",by.cnv=FALSE,lwr=0.25,upr=4,FET=T,prt=F) # regionlist
  oo2 <- extract.cnv.regions(dir,type="dup",by.cnv=FALSE,lwr=0.25,upr=4,FET=T,prt=F) # regionlist
  oo3 <- extract.cnv.regions(dir,type="del",by.cnv=TRUE,lwr=0.25,upr=4,FET=T,prt=F) # cnv list
  oo4 <- extract.cnv.regions(dir,type="dup",by.cnv=TRUE,lwr=0.25,upr=4,FET=T,prt=F) # cnv list
  oo1[["tdt"]] <- NA; oo2[["tdt"]] <- NA
  
  double.del.table <- cnvResults[[4]][which(cnvResults[[4]]$cn==0),] # rareDELs
  double.dup.table <- cnvResults[[5]][which(cnvResults[[5]]$cn==4),] # rareDUPs
  
  oo3 <- update.cnvrs.with.cn(oo3,double.del.table)
  oo4 <- update.cnvrs.with.cn(oo4,double.dup.table,DEL=FALSE)
  
  tdt3 <- make.cnv.reg.snp.matrix(oo3)
  tdt4 <- make.cnv.reg.snp.matrix(oo4)
  
  ped <- get.pedData(ped.file)
  
  tdt3 <- add.all.ids(tdt3, ped, dir)
  tdt4 <- add.all.ids(tdt4, ped, dir)
  
  ii <- tdt.snp(data=ped,snp.data=tdt3)
  tt <- p.value(ii,1)
  oo1[["tdt"]][match(names(tt),rownames(oo1))] <- tt
  ii <- tdt.snp(data=ped,snp.data=tdt4)
  tt <- p.value(ii,1)
  oo2[["tdt"]][match(names(tt),rownames(oo2))] <- tt
  CNVR <- list(deletions=oo1,duplications=oo2)
  sum.del <- ranged.to.data.frame(oo1,T)[order(oo1[[5]]),]
  sum.dup <- ranged.to.data.frame(oo2,T)[order(oo2[[5]]),]
  colnames(sum.del)[c(7,9)] <- colnames(sum.dup)[c(7,9)] <- c("p.fet","p.tdt") 
  ## remove NAs..
  sum.del2 <- sum.del[!is.na(sum.del$p.tdt),]; sum.dup2 <- sum.dup[!is.na(sum.dup$p.tdt),]
  sum.del2[["genes"]] <- substr(sum.del2[["genes"]],1,10); 
  sum.dup2[["genes"]] <- substr(sum.dup2[["genes"]],1,10)
  cat("\nDeletions TDT Results\n") ; print(sum.del2[sum.del2[,"p.tdt"]<.05,])
  cat("\nDuplications TDT Results\n") ; print(sum.dup2[sum.dup2[,"p.tdt"]<.05,])
  ofn=cat.path(dir$res,"TDT_results.RData")
  save(oo1,oo2,sum.del,sum.dup,tdt3,tdt4,ped,file=ofn)
  cat(" wrote TDT family analysis results to",ofn,"\n")
  return(CNVR)
}



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
  new.ids.fail <- all.ids[!all.ids %in% cur.ids & sample.info$QCfail==1]
  CNVRs <- colnames(tdt.cnv)
  out <- as.data.frame(matrix(0,nrow=length(all.ids),ncol=length(CNVRs)))
  rownames(out) <- paste(all.ids); colnames(out) <- paste(CNVRs) 
  out[new.ids.fail,] <- NA
  out <- data.frame.to.SnpMatrix(out)
  out[rownames(tdt.cnv),] <- tdt.cnv[rownames(tdt.cnv),]
  #assume all new.ids.pass == 1 [no need to do it explicitly though]
  return(out)
}


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
        samps <- narm(match(id.lists[[cc]],Xd$id[rel.rows]))
        out[samps,cnvz] <- 2 # homozygous = 2*dup/2*del
      }
      #loop.tracker(cc,length(id.lists))
    }
  }
  return(data.frame.to.SnpMatrix(out))
}


list.to.env <- function(list) {
  if(!is.list(list)) { stop("this function's sole parameter must be a list object")}
  if(is.null(names(list))) { stop("list elements must be named") }
  if(length(list)>1000) { warning("list contains over 1000 elements, this operation will crowd the workspace") }
  for(cc in 1:length(list)) {
    assign(x=names(list)[cc],value=list[[cc]],pos=parent.frame())
  }
  return(NULL)
}


get.ped.linked.sets <- function(ped) {
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
  #add in the snp data # note that missing will get coded as 0 (zero) in this step
  ww[["S51"]] <- as.numeric(tdt3[,"S51"][match(rownames(ww),rownames(tdt3))])
  #choose just the mothers and fathers
  parz <- which(is.na(ww$father) & is.na(ww$mother))
  #choose just the children
  kidz <- which(!is.na(ww$father) & !is.na(ww$mother))
  # key families are those which have a parent with the cnv
  keyfams <- names(which(tapply(ww$S51[parz],factor(ww$familyid[parz]),function(X) { any(X==2) })))
  #just the key families
  with(ww[ww$familyid %in% keyfams,],table(affected,S51))
  # key families are those which have a child with the cnv
  keyfams2 <- names(which(tapply(ww$S51[kidz],factor(ww$familyid[kidz]),function(X) { any(X==2) })))
  #look at just the families with affected kids but not parents
  with(ww[ww$familyid %in% keyfams2 & !ww$familyid %in% keyfams,],table(affected,S51))
  # fix errors in ped file coding
  badz <- with(ped,which(father==1 & mother==1))
  ped[badz,"mother"] <- 2
  badz <- with(ped,which(father==2 & mother==2))
  ped[badz,"father"] <- 1
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

# sumof{for each mum, sum of how many times she passes a del to her children} / numkids
# sumof{for each dad, sum of how many times  he passes a del to his children} / numkids

get.CNV.wise.inheritance.counts <- function(tdt.snp,ped=NULL,only.doubles=FALSE) {
  if(is(tdt.snp)[1]=="SnpMatrix") {  TDT <- SnpMatrix.to.data.frame(tdt.snp) } else { TDT <- tdt.snp }
  if(!is.data.frame(TDT)) { stop("tdt.snp must be a SnpMatrix or data.frame with snp-ids as rownames") }
  if(is.null(ped)) { stop("a valid 'pedData' object must be inputted (see snpStats:tdt.snp documentation)")}
  ped.list <- get.ped.linked.sets(ped)
  list.to.env(ped.list) # take the variables from the list and assign them into the local environment
  if(!exists("kk") | !exists("fklink")) { stop("failure to import variables from get.ped.linked.sets() function") }
  countmat <- matrix(0,nrow=5,ncol=ncol(TDT))
  rownames(countmat)  <- c("mother-pass","father-pass","mother-kidcount","father-kidcount","total")
  colnames(countmat) <- colnames(TDT)
  del.count <- 0
  if(only.doubles) { thresh <- 1 } else { thresh <- 0 }
  for(cc in 1:nrow(mm)) {
    loop.tracker(cc,2*nrow(mm))
    dels <- which(TDT[rownames(mm)[cc],]>thresh) # which DELs does the next mum have
    del.count <- del.count + length(dels) # keep a total count of dels
    if(length(dels)<1) { next }  # skip if this mother has none
    mat <- TDT[rownames(kk)[mklink[[cc]]],dels,drop=FALSE]  # extract her children for these DEL snps
    if(length(Dim(mat))<2) {
      snp.counts <- mat  # in case only 1 child, a row, not a matrix
    } else {
      snp.counts <- colSums(mat,na.rm=T) # add the number of times passed to children to the snp count
    }
    countmat[1,dels] <- countmat[1,dels] + snp.counts
    countmat[3,dels] <- countmat[3,dels] + nrow(mat)
  }
  
  for(cc in 1:nrow(ff)) {
    loop.tracker(cc+1,2*nrow(ff)+1)
    dels <- which(TDT[rownames(ff)[cc],]>thresh) # which DELs does the next dad have
    del.count <- del.count + length(dels) # keep a total count of dels
    if(length(dels)<1) { next }  # skip if this father has none
    mat <- TDT[rownames(kk)[fklink[[cc]]],dels,drop=FALSE]  # extract his children for these DEL snps
    if(length(Dim(mat))<2) {
      snp.counts <- mat  # in case only 1 child, a row, not a matrix
    } else {
      snp.counts <- colSums(mat,na.rm=T) # add the number of times passed to children to the snp count
    }
    countmat[2,dels] <- countmat[2,dels] + snp.counts
    countmat[4,dels] <- countmat[4,dels] + nrow(mat)
  }
  return(countmat)
}






