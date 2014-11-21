
##### FUNCTION INDEX #####
# validate.dir.for - in case the 'dir' input list object is not the standardised list form, convert allows flexible use of list or regular directory specifications in plumbCNV functions
# jlapply
# table2d
# remove.duplicated.id.ranges
# moved to snpMatLst: == > get.SnpMatrix.in.file
# logiti
# FET
# find.overlaps
# get.overlap.stats.chr - prv(next.chr, x, y, xx, yy)
# reduce.list.to.scalars.old
# nmoo
# reduce.list.to.scalars.2
# rlts
# reduce.list.to.scalars
# fill.blank.matches - ensure a blank entry for cnvs with no overlaps full.len is an optional desired out-length for the list (which if trailing cells are blank might differ from the guess) head(paste(names(newlist)))
# overlap.pc
# import.marker.data
# chip.coverage - calculate % of chromosomes/genome covered by a set of microarray markers in terms of possible CNVs detectable targ.int 100000 size of windows to find min.snps minimum number of snps per window can return the gaps - they'll be unmerged...
# calc.cov
# get.dgv.ranges - download or use local version of DGV (database for Genomic Variants)
# draw.cnv.bounds
# read.penn.cnv.file
# rmv.dir.plink.file
# remove.samp.from.plink - remove specific sample from a plink file - eg because had too many cnvs
# plink.to.Ranges
# Ranges.to.cnvgsa
# add.genes.to.GSA
# get.geneData.obj
# full.cnvGSA
# plink.to.cnvGSA
# rmv.dir.penn.cnv.file
# stats.on.CNV.file
# convert.penn.to.plink
# read.plink.file
# extractROH
# convert.snp.indx.to.pos 
# get.ROH.for.SnpMatrix - either use 'dir' to automatically get snp.info and failures
# tdt.snp2
# either.side - retrieve the nearest genes either side of a set of locations, SNPs, CNVs or intervals
# super.annotate.cnvs - add genes, bands, annotate unnamed genes and intergenic ranges,  optionally add haplosufficiency prediction categories
################


### VERY GENERAL ###




# internal function
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



# add genes, bands, annotate unnamed genes and intergenic ranges, 
# optionally add haplosufficiency prediction categories
super.annotate.cnvs <- function(oo1,hap=FALSE) {
  oo1 <- Gene.pos(ranges=oo1,bioC=T,build=36)
  oo1 <- oo1[,-which(colnames(oo1) %in% "index")]
  oo1 <- Band.pos(ranges=oo1,bioC=T,build=36)
  oo1[["gene"]][oo1[["gene"]]== ""] <- "unnamed-gene"
  oo1[["gene"]][oo1[["gene"]] %in% "intergenic"] <- paste0("intergenic [",either.side(oo1[oo1[["gene"]] %in% "intergenic",],build=36,tracker=FALSE),"]")
  oo1[["gene"]][oo1[["gene"]] %in% "unnamed-gene"] <- paste0("unnamed-gene [",either.side(oo1[oo1[["gene"]] %in% "unnamed-gene",],build=36,tracker=FALSE),"]")
  if(hap) {
    JS <- make.hap()
    oo1[["Hap55"]] <- suppressWarnings(hap.mean(oo1,JS,FUN=num.more.than.55,n=.55))
    oo1[["Hap65"]] <- suppressWarnings(hap.mean(oo1,JS,FUN=num.more.than.55,n=.65))
  }
  oo2 <- oo1
  oo2[["gene"]] <- compact.gene.list(oo2[["gene"]])
  return(oo2)
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


# retrieve the nearest genes either side of a set of locations, SNPs, CNVs or intervals
# stored as a ranged data object
# slow so don't run for duplicate ranges, and only run for 'intergenic' ranges or much
# time will be wasted!
either.side <- function(ranges,limit=NULL,build=NULL,collapse=TRUE,sep=",",tracker=TRUE) {
  typ <- is(ranges)[1]
  if(!typ %in% c("GRanges","RangedData","ChipInfo")) { stop("invalid 'ranges' object") }
  chrz <- chr2(ranges)
  stz <- start(ranges)
  enz <- end(ranges)
  left <- right <- character(length(stz))
  midz <- round(rowMeans(cbind(stz,enz)))
  ## preload gene annotation to save time ## 
  ga <- get.gene.annot(build=build,GRanges=F) 
  if(!exists("ga")) { stop("couldn't find gene annotation") }  ## load object: ga [gene database]
  ga <- ga[ga$gene!="",]
  ##
  for (cc in 1:length(stz)) {
    left[cc] <- nearest.gene(chrz[cc], midz[cc], n=1, side="left",ids=TRUE,limit=NULL,build=build, ga=ga) 
    right[cc] <- nearest.gene(chrz[cc], midz[cc], n=1, side="right",ids=TRUE,limit=NULL,build=build, ga=ga) 
    if(tracker) { loop.tracker(cc,length(stz)) }
  }
  if(collapse) {
    return(paste(left,right,sep=sep))
  } else {
    return(data.frame(left=left,right=right))
  }
}


## forces a 2d table with every possible cell (allow zero counts)
#col and row are the categories of which to force tabulation
# cn and rn are optional col and rownames for the output table
table2d <- function(...,col,row,rn=NULL,cn=NULL,remove.na=TRUE) {
  clmn <- paste(col); rowe <- paste(row)
  if(remove.na) { clmn <- narm(clmn); rowe <- narm(rowe); rn <- narm(rn); cn <- narm(cn) }
  inp <- list(...); bad1 <- F 
  if(length(inp)<2) { bad1 <- T } else { 
    if(length(inp[[1]])!=length(inp[[2]])) { bad1 <- T 
    } else { if(length(inp[[1]])<1) { bad1 <- T } }
  }
  if(bad1) {
    warning("at least 2 arguments of equal length required for table2d, passing to table")
    return(table(...)) 
  }
  rrr <- inp[[1]]; ccc <- inp[[2]]
  # add 1 of each possible cell to the 2 vectors, to ensure each table cell is represented
  rrr <- c(rep(sortna(rowe),each=length(clmn)),rrr)
  ccc <- c(rep(sortna(clmn),length(rowe)),ccc)
  #prv(rrr,ccc)
  nmz <- sapply(match.call(expand.dots=TRUE)[-1], deparse) # get ... arg names
  TT <- table(rrr,ccc,dnn=nmz[1:2])-1
  #prv(TT,rn) ; print(cn)
  if(is.character(rn)) { if(length(rn)==length(rowe)) { rownames(TT) <- rn[1:nrow(TT)] } }
  if(is.character(cn)) { if(length(cn)==length(clmn)) { colnames(TT) <- cn[1:ncol(TT)] } }
  return(TT)    
}




# removes all but the first of any overlapping ranges within the same sample
# should be a list of cnvs with a column 'id', e.g, sample id
remove.duplicated.id.ranges <- function(X,column="id") {
  if(column %in% colnames(X)) { 
    coln <- which(colnames(X)==column)
  } else {
    return(X) ; warning("column '",column,"' not found in X, duplicates may still exist") 
  }
  if(!is(X)[1]=="RangedData") { stop("X wasn't RangedData") }
  prn <- nrow(X)
  XX <- vector("list",length(chrNames2(X)))
  do.one.chr <- function(X,coln) {
    idz <- X[[coln]]; 
    idz <- idz[duplicated(idz)]
    udz <- unique(idz)
    n.id <- length(udz)
    #prv(udz)
    for (cc in 1:n.id) {
      select <- which(X[[coln]] %in% udz[cc])
      U <- X[select,]
      ov <- findOverlaps(U)
      qq <- queryHits(ov)
      ss <- subjectHits(ov)  
      tt <- table(qq,ss)
      diag(tt) <- 0
      tt[lower.tri(tt)] <- 0
      to.kill <- colnames(tt)[colSums(tt)>0]
      if(length(to.kill)>0) {
        X <- X[-c(select[as.numeric(to.kill)]),]
      }
      # loop.tracker(cc,n.id); #cat("x")
    }
    #print(X)
    return(X)
  }
  ct <- 1
  for(chrz in chrNames2(X)) {
    XX[[ct]] <- do.one.chr(X[chrz],coln)
    ct <- ct + 1
  }
  X <- do.call("rbind",args=XX)
  ND <- prn - nrow(X)
  if(ND>0) { cat("NOTE: removed",ND,"ids that were duplicates\n") }
  return(X)
}


# specific to me? MOVED TO SNPMATRIXLIST.R
# reads snp matrices from RData binary file. if it finds multiple, will attempt to join them together
# can also handle an XSnpMatrix
#get.SnpMatrix.in.file <- function(file,warn=TRUE){
#  obj.nm <- paste(load(file))
#  ## NOW MAKE SURE WE HAVE EXACTLY ONE SNPMATRIX OBJECT FROM THIS FILE ##
#  if(length(obj.nm)>0) {
#    typz <- sapply(obj.nm,function(X) { is(get(X))[1] })
#    vld <- which(typz %in% c("XSnpMatrix","SnpMatrix"))
#    if(length(vld)<1) { stop("no SnpMatrix objects found in file")}
#    if(length(vld)>1) {
#      if(warn) { warning("found multiple SnpMatrix objects in datafile, attempting to rbind() them") }
#      concat.snp.matrix <- NULL; 
#      vld <- vld[order(names(typz)[vld])] # alphabetical order should ensure consistency across multiple data files
#      try(concat.snp.matrix <- do.call("rbind",args=lapply(obj.nm[vld],function(X) get(X))))
#      if(is.null(concat.snp.matrix)) { stop("SnpMatrix objects had different numbers of SNPs [cols], could not rbind")}
#      obj.nm <- "concat.snp.matrix"
#    } else {
#      obj.nm <- obj.nm[vld]
#    }
#  }
#  ret.out <- get(obj.nm[1])
#  if(is(ret.out)[1]=="XSnpMatrix") { if(warn) { warning("read XSnpMatrix from ",file) } }
#  return(ret.out)
#}



### MORE GENERAL ###


# specific
# basically same format as the fishers test but uses logistic regression (more power)
logiti <- function(case,ctrl,case.d=NA,cont.d=NA,inclusive=T,...,stat=c("p.value","conf.int","estimate","all")) { 
  stat <- tolower(paste(stat[1]))
  if(inclusive) {
    case.d <- case.d-case
    cont.d <- cont.d-ctrl
  }
  if(length(case)>1 | length(ctrl)>1) { case <- case[1]; ctrl <- ctrl[1]; warning("only first elements of count vectors will be used") }
  if(!stat %in% c("p.value","conf.int","estimate","all")) { 
    warning("invalid stat value, setting to p.value"); stat <- "p.value" }
  case.vec <- c(rep(0,(case.d)),rep(1,case))
  ctrl.vec <- c(rep(0,(cont.d)),rep(1,ctrl))
  ph <- c(rep(0,length(ctrl.vec)),rep(1,length(case.vec)))
  cnv <- c(ctrl.vec,case.vec)
 # print(table(cnv,ph))
  res <- glm(ph~cnv,family=binomial("logit"))
 # return(res) ##### HERE!!!!
  out <- (logistic.summary(res,ci=T))
  #OR OR-low OR-hi  p-value
  if(stat=="conf.int") { 
    Out <- paste(round(as.numeric(out[[1]][,"OR-low"]),3),
     round(as.numeric(out[[1]][,"OR-hi"]),3),sep=",") 
    #return(out)
  } else {
    if(stat=="estimate") {
      Out <- out[[1]][,"OR"]
    } else { 
      if(stat=="all") {
        Out <- out
      } else {
        Out <- (out[[1]][,"p-value"])
      }
    }
  }
  return(Out)
}

# specific
# fishers exact test given case and control obs counts, then total group counts from info object or manual counts.
FET <- function(case,ctrl,dir=NULL,sample.info=NULL,case.d=NA,cont.d=NA,stat=c("p.value","conf.int","estimate","all"),inclusive=T) {
  skip.c <- F
  DIGITS <- 3
  stat <- tolower(paste(stat[1]))
  if(!stat %in% c("p.value","conf.int","estimate","all")) { 
    warning("invalid stat value, setting to p.value"); stat <- "p.value" }
  if(is.null(sample.info)) {
    if(!is.null(dir)) {
      sample.info <- read.sample.info(dir)
    } else {
      if(!all(!is.na(c(case.d,cont.d)))) {
        warning("must provide case.d, cont.d, or 'dir' or sample.info to ",
                "calculate the case/control denominators"); return(NULL)
      } else {
        skip.c <- T
      }
    }
  }
  if(!skip.c) {
    p.c <- table(sample.info$phenotype,sample.info$QCfail)[,1]
    case.d <- p.c[2]; cont.d <- p.c[1]
  }
  if(stat=="all") { fet <- vector("list",length(case)) } else { fet <- numeric(length(case)) }
  for (cc in 1:length(case)) {
    if(length(case[cc])>0 & length(ctrl[cc])==0) { ctrl[cc] <- 0 }
    if(length(case[cc])==0 & length(ctrl[cc])>0) { case[cc] <- 0 }
    if(length(case[cc])>0 & length(ctrl[cc])>0) { 
      if(inclusive) {
        case.d2 <- case.d-case[cc]
        cont.d2 <- cont.d-ctrl[cc]
      } else { case.d2 <- case.d; cont.d2 <- cont.d }
      ftr <- fisher.test(x=matrix(c(case[cc],case.d2,ctrl[cc],cont.d2),nrow=2))
      FTR <- ftr[[stat]]
      if(stat=="conf.int") { fet[cc] <- paste(round(FTR,digits=DIGITS),collapse=",") } else { fet[cc] <- FTR }
      if(stat=="all") { fet[[cc]] <- ftr }
    } else {
      if(stat=="all") { fet[[cc]] <- NA } else { fet[cc] <- NA }
    }
    #loop.tracker(cc,length(case))
  }
  return(fet)
}



# specific to plumb cnv
# like 'findOverlaps' but according to a specific percentage
# specify own reference comparison or use standard 'gene', 'exon' or 'dgv'
# can filter on several criteria prior to matching (length/snps/dels/dups), or after matching (overlap %)
# can get a list of overlap hits, or just a filtered RangedData object
find.overlaps <- function(cnv.ranges, thresh=0, geq=T, rel.ref=T, pc=T, ranges.out=F, vals.out=F, vec.out=F, delim=",",
                          ref=NULL, none.val=0, fill.blanks=T, DEL=T, DUP=T, min.sites=0, len.lo=NA, len.hi=NA,
                          autosomes.only=T, alt.name=NULL, testy=F,
                          db=c("gene","exon","dgv"), txid=F, build="hg18", n.cores=1, dir="", quiet=F) {
  if(is(cnv.ranges)[1]!="RangedData") {
    warning("'cnv.ranges' must be 'RangedData' type; returning null")
    return(NULL)
  } else {
    if(nrow(cnv.ranges)<1) { warning("cnv.ranges contains zero rows, returning NULL"); return(NULL) }
  }
  build <- ucsc.sanitizer(build)
  uv <- tolower(universe(cnv.ranges)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { build <- uv } }
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
  if(is(ref)[1]=="GRanges") { ref <- as(ref,"RangedData") }
  if(is(ref)[1]=="RangedData" & all(db==c("gene","exon","dgv"))) {
    db <- "custom ranges"
  }
  db <- tolower(paste(db)[1])
  if(is(ref)[1]!="RangedData") {
    ; ano <- "gene" 
    if(length(grep("exon",db))>0) { ano <- "exon" }
    if(length(grep("dgv",db))>0) { ano <- "dgv" }
    if(!quiet) { cat(" loading",ano,"annotation...\n") }
    ref <- switch(ano,exon=get.exon.annot(dir=dir,build=build,GRanges=FALSE),gene=get.gene.annot(dir=dir,build=build,GRanges=FALSE),
                  dgv=get.dgv.ranges(dir=dir,build=build,compact=T,GRanges=FALSE)) # choose which reference to use
    
  } else {
    # ref is already a ranges object for comparison passed in by parameter
    if(any(tolower(db) %in% c("exon","gene","dgv")))  { ano <- db } else { ano <- "custom ranges" ; txid <- F }
  }
  do.check <- testy
  if(ano=="gene" | ano=="exon") { nbg <- T } else { nbg <- F } # use gene names if gene or exon
  if(ano=="dgv") { if(!(DUP & DEL)) {
    if(!DUP) { if("Loss" %in% colnames(ref)) { ref <- ref[!is.na(ref$Loss),] } else { cat("Loss column not found\n")} }
    if(!DEL) { if("Gain" %in% colnames(ref)) { ref <- ref[!is.na(ref$Gain),] } else { cat("Gain column not found\n")} }
  } }
  # rel.query=T gives %'s relative to the query (e.g gene), false gives %'s relative to subject (ie, cnv)
  if(!quiet) { cat(" testing CNV set for overlaps with",ano,"...\n") }
 # if(do.check) { prv(ref,cnv.ranges) }
  mm <- overlap.pc(query=ref,subj=cnv.ranges,name.by.gene=nbg,
                   rel.query=rel.ref,fill.blanks=fill.blanks, txid=txid,alt.name=alt.name,
                   n.cores=n.cores,none.val=none.val,autosomes.only=autosomes.only)

#  if(do.check) { prv(mm) }
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
    pass.thresh <- jlapply(list1=mm[[3]],list2=thresh,pc=pc,FUN=fun) # get list of passing threshold
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

# internal specific to plumb
## mainly to tidy up the function below - does the processing for 1 chromosome
get.overlap.stats.chr <- function (next.chr, x, y, xx, yy, name.by.gene=T, rel.query=T, txid=F, prog=F, alt.name=NULL) {
  #prv(next.chr, x, y, xx, yy)
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
     # if(next.chr %in% c(1,11)) { prv(xx[next.chr]) }
      ind <- rownames(xx[next.chr])
      genes.per.cnv <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
    } else {
      if(!is.null(alt.name)) {
        #prv.large(xx[next.chr])
        ind <- xx[next.chr][,paste(alt.name)]
        genes.per.cnv <- sapply(genes.per.cnv.n,match.num.to.names,simplify=F,ind=ind)
      } else {
        # use row-numbers for matchs unless using DGV ids
        genes.per.cnv <- genes.per.cnv.n 
      }
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


# internal plumb
# flattens list - old version
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

# make null vector
nmoo <- function(x) { x <- rep(NULL,1) }

# flattens list 
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

# gated list flattening
rlts <- function(X) {
  if(is(X)[1]=="list") { 
    X <- reduce.list.to.scalars(X) 
  } else {
    if(length(X)>1 | !is.numeric(X)) { X <- NULL }
  }
  return(X)
}

# flattens list 
reduce.list.to.scalars <- function(ll) {
  if(is.list(ll)) {
    ll <- rev(lapply(rev(ll),rlts))
  }
  return(ll)
}


# for all non overlaps enter some missing value to preserve length of output list
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


# uses 'findOverlaps' and reports the percentage overlap for each region/etc
# most sense when: x is reference (e.g, genes); y is regions of interest, eg.  CNVs
overlap.pc <- function(query,subj,name.by.gene=F,rel.query=T,fill.blanks=T,autosomes.only=T,
                       text.out=F,delim=",",none.val=0,n.cores=1, txid=F, prog=F, alt.name=NULL) {
  must.use.package(c("genoset","IRanges"),T)
  if(is(query)[1]!="RangedData" | is(subj)[1]!="RangedData") {
    warning("'query' and 'subj' must both be 'RangedData' type; returning null")
    return(NULL)
  }  
  blnk.out <- list(ids=NULL,width=NULL,pc=NULL)
  if(nrow(query)<1 | nrow(subj)<1) { return(blnk.out) }
  # force this so order can be constructed with respect to the initial order
  subj.nrow <- nrow(subj)
  #  prv(query,subj)
  rownames(subj) <- paste(1:subj.nrow)
  # query names only matter for result text, so only add names if no text.
  # genes/exons have a separate column treated specially as otherwise there are duplicate entries not suited to rownames
  #print(length(which(is.na(rownames(query)))))
  if(is.null(rownames(query))) { rownames(query) <- paste(1:nrow(query)) }
  ## reduce dimension of subj; 
  to.cut <- which(tolower(colnames(subj)) %in% c("cn","numsnps","fid","score"))
  if(length(to.cut)>0) { subj <- subj[,-to.cut] } # removing unnecessary cols speeds up
  to.cut2 <- which(!tolower(colnames(query)) %in% c("gene","txid","txids","txname","txnames",alt.name))
  if(length(to.cut2)>0) { query <- query[,-to.cut2] } # removing unnecessary cols speeds up
  #print(head(query))
  xlist <- set.chr.to.numeric(query,table.out=T)
 # print(head(xlist[[1]]))
  #  prv(xlist)
  xx <- xlist[[1]]
  xx=toGenomeOrder2(xx,strict=T)
  #print(head(xx))
  ili <- set.chr.to.numeric(subj,table.in=xlist[[2]],table.out=F)
  yy=toGenomeOrder2(ili,strict=T)
  if(autosomes.only) {
    xx=select.autosomes(xx)
    yy=select.autosomes(yy)
  }
  #  xx=toGenomeOrder(select.autosomes(query),strict=T)
  #  yy=toGenomeOrder(select.autosomes(subj),strict=T)
  chr.set.x <- chrNums(xx)
  chr.set.y <- chrNums(yy)
  #  prv(xx,yy,chr.set.x,chr.set.y)
  common <- which(sortna(chr.set.x) %in% sortna(chr.set.y))
  if(length(common)<1 | nrow(xx)<1 | nrow(yy)<1) { return(blnk.out) }
  chr.set <- paste(sortna(chrNums(xx[paste(sortna(chr.set.x)[common])])))
  #print(head(rownames(xx))); print(head(rownames(yy)),7)
  x <- ranges(xx[chr.set]); y <- ranges(yy[chr.set]) # keep only common + convert to IRangeslist
  xx <- xx[chr.set]; yy <- yy[chr.set] # keep only common in originals (these store rownames, unfiltered positions)
  #print(head(rownames(xx))); print(head(rownames(yy)))
  n.chr <- length(chr.set)
  by.chr.list <- vector("list",n.chr); names(by.chr.list) <- chr.set
  if(prog) { cat("|") }
  if(n.cores>1) {
    must.use.package("parallel")
    by.chr.list <- parallel::mclapply(X=1:n.chr, FUN=get.overlap.stats.chr, x=x, y=y, xx=xx, yy=yy, alt.name=alt.name,
                                       name.by.gene=name.by.gene, rel.query=rel.query, txid=txid, mc.cores=n.cores,prog=prog)
  } else {
    by.chr.list <- lapply(X=1:n.chr, FUN=get.overlap.stats.chr, x=x, y=y, xx=xx, yy=yy,  alt.name=alt.name,
                          name.by.gene=name.by.gene, rel.query=rel.query, txid=txid, prog=prog)
  }
  if(prog) { cat("|\n") }
  # prv(by.chr.list)
  #  return(by.chr.list) ; # save(by.chr.list,file="by.chr.list.RData")
  namez <- as.character(unlist(sapply(lapply(by.chr.list,"[[",1),names)))
  #  namez <- as.character(unlist(sapply(by.chr.list,names)))
  idz <- unlist(sapply(by.chr.list,"[[",1),recursive=F)
  if(all(is.null(idz))) { if(name.by.gene) { warning("blank names: nb: name.by.gene should be false if gene names not present") } }
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


# iFunctions? probably not, more plumbCNV specific
#create a snp.info object
# calculate chromosome-wise position and ID for each SNP in list()
# dir should be contain the annotation directory assumed to contain the snp list (dir.ano)
import.marker.data <- function(dir, markerinfo.fn="snpdata.map",snp.fn="snpNames.txt", anot=c("bim","vcf","map","map3")[4],
                               snp.col=NA, pos.col=NA, chr.col=NA, verbose=F )
{
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
  snpData <- toGenomeOrder2(snpData,strict=T)
  return(snpData)
}





### LESS GENERAL ###

# calculate relative chip coverage vs gwas
chip.coverage <- function(snp.info,targ.int=100000,by.chr=F,min.snps=10,add.missing.chr=TRUE,
                          full.chr.lens=T,verbose=T,dir="",build="hg18",ranges=F) {
  # calculate % of chromosomes/genome covered by a set of microarray markers in terms of
  # possible CNVs detectable
  # targ.int #100000  ## size of windows to find
  # min.snps ## minimum number of snps per window
  # can return the gaps - they'll be unmerged...
  must.use.package(c("genoset","IRanges"),T)
  build <- ucsc.sanitizer(build)
  if(!is(snp.info)[1]=="RangedData") { warning("snp.info wasn't RangedData") ; return(NULL) }
  snp.info <- select.autosomes(snp.info)
  snp.info <- toGenomeOrder2(snp.info,strict=T)
  chr.set <- chrNums(snp.info); n.chr <- length(snp.info)
  # single core -  multicore doesn't really help much (like 5 seconds saving for 22 cores for a large dataset)
  rd <- calc.cov(snp.info, targ.int, min.snps)
  rd <- reduce(rd) # combine adjacent ranges
  # calculate total length as total length of all chromosomes
  if(is.null(dir)) { dir <- "" }
  chr.lens.full <- (get.chr.lens(dir,build=build)[chr.set])-(targ.int)
  # calculate total  length using the first and last pos on each chromosome [allows calc for subranges]
  chr.lens.range <- as.numeric(sapply(snp.info,function(x) { diff(range(start(x))) })) #- (targ.int)
  if (full.chr.lens) { 
    # add gaps before the first snp and after the last snp for each chr
    chr.lens <- chr.lens.full 
    stz <- sapply(snp.info,function(x) { min(start(x)) }) #- (targ.int)
    enz  <- (chr.lens.full - sapply(snp.info,function(x) { rev(sortna(start(x)))[min.snps-1] })) - (targ.int)
    rd2 <- RangedData(IRanges(start=rep(1,n.chr)[stz>0],end=stz[stz>0]),space=chr.set[stz>0])
    rd3 <- RangedData(IRanges(start=(chr.lens.full-enz-targ.int)[enz>0],end=(chr.lens.full-targ.int)[enz>0]),space=chr.set[enz>0])
    rd  <- rbind(toGenomeOrder2(rd,strict=T),toGenomeOrder2(rd2,strict=T),toGenomeOrder2(rd3,strict=T))
  } else { 
    cat("using range restricted by start and end of chip coverage for each chromosome")
    chr.lens <- chr.lens.range  
  }
  rd <- toGenomeOrder2(rd,strict=T)
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


# calculate relative chip coverage vs gwas
chip.coverage2 <- function(snp.info,targ.int=100000,min.snps=10,by.chr=F,pad.missing.autosomes=TRUE,
                           verbose=T,dir="",build="hg18",ranges=F,alt.ref=NULL) {
  # calculate % of chromosomes/genome covered by a set of microarray markers in terms of
  # possible CNVs detectable
  # targ.int #100000  ## size of windows to find
  # min.snps ## minimum number of snps per window
  # can return the gaps - they'll be unmerged...
  if(!is.null(alt.ref) & !by.chr) {
    if(is(alt.ref)[1]!=is(snp.info)[1]) { stop("alt.ref should be NULL, or same type as snp.info") }
    if(length(grep("chr",chrNames2(snp.info)))>0) { alt.ref <- set.chr.to.char(alt.ref) } else { alt.ref <- set.chr.to.numeric(alt.ref) }
    snp.info.sub <- subsetByOverlaps(snp.info,alt.ref)
    #prv(snp.info.sub)
    covered <- calc.cov2(snp.info.sub, targ.int, min.snps)
    cov <- sum(as.numeric(width((covered))))
    total.ref <- sum(as.numeric(width(alt.ref)))
    cat("Base pairs covered: ",cov,"\n")
    cat("Base pairs of gap: ",total.ref-cov,"\n")
    cat("Coverage % of alternative reference, bp=",(total.ref),": ",cov/(total.ref),"\n")
    return(cov/(total.ref))
  } else {
    covered <- calc.cov2(snp.info, targ.int, min.snps)
    not.covered <- invert.granges(as(covered,"GRanges"),pad.missing.autosomes=pad.missing.autosomes,build=build)
  }
  if(ranges) { return(covered) }
  if(!by.chr) {
    cov <- sum(as.numeric(width((covered))))
    not.cov <- sum(as.numeric(width((not.covered))))
    cat("Base pairs covered: ",cov,"\n")
    cat("Base pairs of gap: ",not.cov,"\n")
    pc <- cov/(cov+not.cov)
    cat("Coverage % of genome bp=",(cov+not.cov),": ",pc,"\n")
    return(pc)
  } else {
    covered <- as(covered,"RangedData")
    not.covered <- as(not.covered,"RangedData")
    chr.cov <- chr.gap <- chr.pc <- vector("list",length(covered))
    names(chr.pc) <- names(chr.gap) <- names(chr.cov) <- chrNames(covered)
    for (cc in 1:length(chr.cov)) {
      chr.cov[[cc]] <- sum(as.numeric(width((covered[cc]))))
      chr.gap[[cc]] <- sum(as.numeric(width((not.covered[cc]))))
      chr.pc[[cc]] <- chr.cov[[cc]]/(chr.cov[[cc]]+chr.gap[[cc]])
    }
    return(list(pc=chr.pc,cov=chr.cov,gap=chr.gap))
  }
}
  
calc.cov2 <- function (snp.info, targ.int, min.snps) {
  typ <- is(snp.info)[1]
  if(!typ %in% c("RangedData","GRanges","ChipInfo")) { stop("snp.info must be GRanges, RangedData, or ChipInfo type") }
  snp.info <- toGenomeOrder2(snp.info)
  st <- start(snp.info)
  ch <- chr2(snp.info)
  tch <- table(ch)
  if(sum(as.numeric(tch))<1) { stop("no SNPs found in snp.info") }
  chr.list <- names(tch)[as.numeric(tch)>0]
  nchr <- length(chr.list)
  chr.cov <- vector("list",nchr)
  names(chr.cov) <- chr.list
  if(length(nchr)<1) { prv(snp.info); stop("no chromosomes found in snp.info") }
  for (cc in 1:nchr) {
    #print(cc)
    nch <- chr.list[cc]
    sel <- ch==nch
    nsnp <- length(st[sel])  
    if(nsnp<(min.snps+1)) { 
      warning("chr",nch," had too few snps for any coverage")
      chr.cov[[cc]] <- RangedData(IRanges(start=1,end=1),space=nch)
    } else {
      chz <- ch[sel][1:(nsnp-(min.snps-1))]
      stz <- st[sel][1:(nsnp-(min.snps-1))]
      enz <- st[sel][min.snps:nsnp]
      if(any(stz>enz)) { stop("start should never be > end") }
      range6 <- RangedData(IRanges(start=stz,end=enz),space=chz)   
      toolong <- width(range6)>targ.int
      chr.cov[[cc]] <- reduce(range6[!toolong,])
    }
    #loop.tracker(cc,nchr)
  }
 # print(sapply(chr.cov,is))
  #return(chr.cov)
  out.obj <- toGenomeOrder2(do.call("rbind",args=chr.cov))
  out.obj <- as(out.obj,typ)
  return(out.obj)
}


# internal function used by chip.coverage() to calculate coverage gaps
calc.cov <- function (snp.info, targ.int, min.snps) {
  cur.chr <- chr2(snp.info)
  if(nrow(snp.info)<(min.snps+1)) { warning("not enough snps to calculate coverage"); return(NA) }
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
        skipnum <- 0
        while(dif[cc]<gap) {
          # see how many snps we need to skip ahead until the gap to 10th [min.snp] snp is <= targ.int
          # record this distance as a gap where a CNV can't be detected, then start again from the 'skipnum' snp
          skipnum <- skipnum+1
          dif[cc] <- startz[cc+skipnum] - startz[cc]   ## old idea, don't think it works : dif2 <- startz[cc+min.snps+skipnum+1] - startz[min.snps+cc]
        }
        beg[cc] <- startz[cc]; fin[cc] <- startz[cc+skipnum]; Chr[cc] <- cur.chr[cc]
      } else {
        # otherwise 'min.snps' were in an interval <=targ.int, so all is kosher
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
  rd <- toGenomeOrder2(rd,strict=T)
  return(rd)
}


get.genome <- function(build=NULL,num.names=TRUE) {
  ccc <- get.chr.lens(build=build,names=T)
  ch <- names(ccc); if(num.names) { ch <- gsub("chr","",ch) }
  ii <- RangedData(IRanges(start=rep(1,length(ccc)),end=as.numeric(ccc)),space=ch)
  return(toGenomeOrder2(ii))
}


#dir <- make.dir("/chiswick/data/ncooper/ImmunochipFamilies")
#hh550 <- reader("~/Downloads/hh550.map")
#colnames(hh550) <- c("chr","pos")
#hh550 <- data.frame.to.ranged(hh550)
#hh550 <- toGenomeOrder2(subsetByOverlaps(hh550, get.genome(build=36))) # make sure no illegal SNP positions!
#ichip.regions <- reader("/chiswick/data/ncooper/imputation/COMMON/iChipFineMappingRegionsB36.RData")
#snp.info <- read.snp.info(dir)
#i12 <- set.chr.to.numeric(ichip.regions)
#t1d <- get.t1d.subset(i12,build=36)
#sum <- make.summary.of.ichip.vs.550.coverage(hh550,snp.info,alt.ref2,alt.ref3)
  
make.summary.of.ichip.vs.550.coverage <- function(hh550,snp.info,t1d,i12,build="hg18") {
  #source("~/github/plumbCNV/generalCNVFunctions.R")
  cnv.sizes <- c(5000,10000,20000,50000,100000,200000,400000,1000000,3000000)
  regs <- c("T1D","Immune12","Genome")
  chips <- c("ImmunoChip","HumanHapMap550K")
  cols <- paste(rep(regs,length(chips)),rep(chips,each=length(regs)),sep="-")
  result.mat <- matrix(nrow=length(cnv.sizes),ncol=length(cols))
  colnames(result.mat) <- cols
  rownames(result.mat) <- cnv.sizes
  snp.infoA <- select.autosomes(snp.info)
  hh550A <- select.autosomes(hh550)
  for(cc in 1:length(cnv.sizes)) {
    Header(paste0(cnv.sizes[cc]/1000,"K minimum CNV size"))
    result.mat[cc,1] <-chip.coverage2(snp.infoA,cnv.sizes[cc],6,alt.ref=t1d,build=build) # T1D regions coverage for iChip
    result.mat[cc,2] <-chip.coverage2(snp.infoA,cnv.sizes[cc],6,alt.ref=i12,build=build) # immune regions coverage for iChip
    result.mat[cc,3] <-chip.coverage2(snp.infoA,cnv.sizes[cc],6,build=build)                  # genome coverage for iChip
    
    result.mat[cc,4] <-chip.coverage2(hh550A,cnv.sizes[cc],6,alt.ref=t1d,build=build) # T1D regions coverage for 550k
    result.mat[cc,5] <-chip.coverage2(hh550A,cnv.sizes[cc],6,alt.ref=i12,build=build) # immune regions coverage for 550k
    result.mat[cc,6] <-chip.coverage2(hh550A,cnv.sizes[cc],6,build=build)                  # genome coverage for 550k
  }
  
  print(result.mat,quote=F)
  return(result.mat)
}


# import the DGV into a rangedData object
get.dgv.ranges <- function(dir=NULL,build="hg18",bioC=TRUE,text=FALSE,shortenID=TRUE, compact=FALSE, alt.url=NULL, GRanges=TRUE)
{
  ## download or use local version of DGV (database for Genomic Variants)
  from.scr <- TRUE
  build <- ucsc.sanitizer(build)
  dir <- validate.dir.for(dir,c("ano"),warn=FALSE)
  old.colnm.core <- c("VariationID","Start","End","Chr","Gain","Loss") # order: ids,st,en,chr,...
  colnm.core <- c("variantaccession","start","end","chr","observedgains","observedlosses") # order: ids,st,en,chr,...
  if(!is.null(dir)) {
    dg.fn <- cat.path(dir$ano,"dgvAnnot.RData")
    if(file.exists(dg.fn)) {
      tt <- get(paste(load(dg.fn)))
      from.scr <- FALSE
    }
  }
  if(from.scr) {
    #dgv.path <- paste("http://projects.tcag.ca/variation/downloads/variation.",build[1],".v10.nov.2010.txt",sep="")
    dgv.path <- switch(build,hg19="http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2013-05-31.txt",
                       hg17="http://dgv.tcag.ca/dgv/docs/NCBI35_hg17_variants_2013-05-31.txt",
                       hg18="http://dgv.tcag.ca/dgv/docs/NCBI36_hg18_variants_2013-05-31.txt")
    if(!is.null(alt.url)) { dgv.path <- alt.url }
    downloc <- cat.path(dir$ano,"DGV.variants.txt")
    download.file(url=dgv.path,downloc,quiet=TRUE)
    tt <- read.delim(downloc,header=TRUE,stringsAsFactors=FALSE)
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
                                    width=NULL,chr=colnm.core[4],build=build, exclude=to.cut)
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
  if(bioC & GRanges) {
    outData <- as(outData,"GRanges")
  }
  return(outData)
}



# do abline at start and end of all CNV positions
# ie to make a graph of a CNV easier to read
draw.cnv.bounds <- function(cnv,chr.offset=NA,pos=NA,cnv.lty="dashed",cnv.col="orange",plot.scl=10^6) {
  done <- F
  if(is(cnv)[1]=="RangedData") {
    # highlight CNVs in rangedData object
    if(all(is.na(chr)) & all(is.na(pos))) {
      cnv <- in.window(cnv,chr,pos,full.overlap=F, rmv.dup=T)
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

# read in a penn cnv file (with standard suffixes)
read.penn.cnv.file <- function(filename,readtable=TRUE) {
  cN <- c("coords","n.snps","length","copy.number","file","first.snp","last.snp")
  if(!file.exists(filename)) { stop("specified penn-cnv file did not exist") }
  if(readtable) {
    file.out <- read.table(filename,header=FALSE,stringsAsFactors=FALSE)
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


# remove directories from filenames in a plink file
rmv.dir.plink.file <- function(...)
{
  rmv.dir.penn.cnv.file(...,plink=T)
}

# remove specifix sample list from a plink file, e.g, after applying Qc
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

# read in a plink cnv file and return Ranged data object
plink.to.Ranges <- function(plink.cnv) {
  if(get.ext(plink.cnv)!="cnv") {
    warning(paste("this function is meant for importing plink *.cnv files so using a file with",
                  "extension",get.ext(plink.cnv),"may cause unpredictable results\n"))
  }
  dat <- read.table(plink.cnv,header=TRUE,stringsAsFactors=FALSE)
  # Plink format:
  #FID  IID  CHR  BP1  BP2  TYPE	SCORE	SITES
  outData <- RangedData(ranges=IRanges(start=dat$BP1,end=dat$BP2),id=dat$IID,space=dat$CHR,
                        cn=dat$TYPE,numSnps=dat$SITES,fid=dat$FID,score=dat$SCORE)
  outData <- toGenomeOrder2(outData,strict=T)
  # note the IDs can't be rownames due to duplicates
  return(outData)
}


# convert an Ranged data object to a cnv-gsa object
# for use with the package cnvGSA
Ranges.to.cnvgsa <- function(dat) {
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


## use cnvGSA function 'getCnvGenes' to add gene annotation to a cnv slot object for cnvGSA
add.genes.to.GSA <- function(cnv.frame, delim=";", build="hg18", genemap=NULL) {
  gsa.cols <- c("SampleID","Chr","Coord_i","Coord_f","Type","Genes","CnvID")
  build <- ucsc.sanitizer(build)
  if(!all(colnames(cnv.frame)==gsa.cols))
  { cat(" invalid frame, genes not added\n") ; return(cnv.frame) }
  must.use.package("cnvGSA",T)
  if(is.null(genemap)) {
    genemap <- get.gene.annot(build=build,bioC=F,GRanges=FALSE)
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

# get object needed for cnvGSA package
# contains lookup of SYMBOL vs Entrez ID, and descriptions of genes
get.geneData.obj <- function() {
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

# create a cnvGSA object using cnv.frame or RangedData object generated by plumb cnv
full.cnvGSA <- function(do.assoc=F, cnv.frame=NULL, gmt.file=NULL, grandtotals.mode=c("all","cnv","cnvGen"),
                        sample.classes=c("case","ctrl"),fdr.iter=2,ext.report=200,boxplots=F,
                        limit.type=c("DEL","DUP"),limits.size=NULL,rem.genes=NULL,
                        delim=";", build="hg18", genemap=NULL) {
  must.use.package("cnvGSA",T)
  build <- ucsc.sanitizer(build)
  if(is(cnv.frame)[1]=="RangedData") { 
    uv <- tolower(universe(cnv.frane)); if(length(uv)>0) { if(uv %in% paste("hg",16:20,sep="")) { build <- uv } }
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
  cnv.frame <- add.genes.to.GSA(cnv.frame, delim=";", build=build, genemap=NULL)
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

# read in a plink cnv file and return cnvGSA object
plink.to.cnvGSA <- function(cnv.data) {
  if(get.ext(cnv.data)!="cnv") {
    warning(paste("this function is meant for importing plink *.cnv files so using a file with",
                  "extension",get.ext(cnv.data),"may cause unpredictable results"))
  }
  dat <- plink.to.Ranges(cnv.data)
  dat <- Ranges.to.cnvGSA(dat)
  return(dat)
}

# whether to read files using 'read.table()' or manual method.. speed varies
rmv.dir.penn.cnv.file <- function(filename,append=".nodir",ext=F,verbose=F,
                                  plink=F,readtable=T) {
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


## return detailed summary of CNVs in a plink .cnv file
stats.on.CNV.file <- function(fn,use.dat=F) 
{
  dat <- read.table(fn,header=TRUE,stringsAsFactors=FALSE)  
  tt <- as.numeric(table(dat$IID))[order(as.numeric(table(dat$IID)))]
  # headers in file
  #FID  IID  CHR  BP1  BP2	TYPE	SCORE	SITES
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

## convert from penn cnv format to plink format
## plink has tab delimited and fixed width formats. can do either
convert.penn.to.plink <- function(penn.in,plink.out=NULL,fixed.width=F) {
  penn <- read.penn.cnv.file(penn.in)
  if(ncol(penn)<1 | nrow(penn)<1) { stop("file ",penn.in," was empty, returning NULL") }
  colnames(penn) <- c("location","numsnp","length","cn","id","startsnp","endsnp")
  #prv(penn$location)
  loc <- as.data.frame(convert.textpos.to.data(penn$location),stringsAsFactors=FALSE)
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




# read in a plink format CNV file
read.plink.file <- function(filename,readtable=TRUE) {
  cN <- c("FID","IID","CHR","BP1","BP2","TYPE","SCORE","SITES")
  if(readtable) {
    file.out <- read.table(filename,header=T,stringsAsFactors=FALSE)
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
      file.out <- read.table(filename,stringsAsFactors=FALSE)
    }
  }
  return(file.out)
}



#### ROH functions

# core function to find runs of ROH
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

# for an ROHlist, converts raw snp index information back into a position
convert.snp.indx.to.pos <- function(ROHlist,snpMat,snp.info) {
  #must.use.package("genoset")
  snp.info <- toGenomeOrder2(select.autosomes(snp.info),strict=T)
  valid.snps <- colnames(snpMat)[colnames(snpMat) %in% rownames(snp.info)]
  select <- which(rownames(snp.info) %in% valid.snps)
  if(length(valid.snps)<ncol(snpMat)) { warning("some snp names from snpMat were not found in snp.info - set to NA") }
  snp.info <- snp.info[select,]  
  st <- start(snp.info); en <- end(snp.info); ch <- chr2(snp.info)
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
  snp.info <- toGenomeOrder2(select.autosomes(snp.info),strict=T)
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


## copyof snpStats::tdt.snp = don't think it's used anymore
tdt.snp2 <- function (ped, id, father, mother, affected, data = sys.parent(), 
                      snp.data, rules = NULL, snp.subset = NULL, check.inheritance = TRUE, 
                      robust = FALSE, uncertain = FALSE, score = FALSE) 
{
  if (!is.null(rules) || !check.inheritance) 
    robust <- TRUE
  mcall <- match.call()
  if (!is(snp.data, "SnpMatrix")) 
    stop("snp.data argument must be of class SnpMatrix")
  if (missing(data)) {
    ped <- as.character(ped)
    nped <- length(ped)
    if (nped != nrow(snp.data)) 
      stop("length of `ped' argument incompatible with `snp.data'")
    id <- as.character(id)
    if (length(id) != nped) 
      stop("incompatible length for `id' and `ped' arguments")
    father <- as.character(father)
    if (length(father) != nped) 
      stop("incompatible length for `father' and `ped' arguments")
    mother <- as.character(mother)
    if (length(mother) != nped) 
      stop("incompatible length for `mother' and `ped' arguments")
    affected <- as.logical(affected)
    if (length(affected) != nped) 
      stop("incompatible length for `affected' and `ped' arguments")
    subject.names <- rownames(snp.data)
    in.snp <- 1:nped
    have.snps <- rep(TRUE, nped)
    print("here")
  }
  else {
    data <- as.data.frame(data)
    nped <- nrow(data)
    subject.names <- rownames(data)
    if (missing(ped)) 
      ped <- as.character(data[, 1])
    else ped <- as.character(eval(mcall$ped, envir = data))
    if (is.null(ped)) 
      stop("pedigree identifiers not found in data frame")
    if (missing(id)) 
      id <- as.character(data[, 2])
    else id <- as.character(eval(mcall$id, envir = data))
    if (is.null(id)) 
      stop("subject identifiers not found in data frame")
    if (missing(father)) 
      father <- as.character(data[, 3])
    else father <- as.character(eval(mcall$father, envir = data))
    if (is.null(father)) 
      stop("father identifiers not found in data frame")
    if (missing(mother)) 
      mother <- as.character(data[, 4])
    else mother <- as.character(eval(mcall$mother, envir = data))
    if (is.null(mother)) 
      stop("mother identifiers not found in data frame")
    if (missing(affected)) 
      affected <- (data[, 6] == 2)
    else affected <- as.logical(eval(mcall$affected, envir = data))
    if (is.null(affected)) 
      stop("disease status not found in data frame")
    in.snp <- match(subject.names, rownames(snp.data))
    have.snps <- !is.na(in.snp)
    print("there")
  }
  affected[is.na(affected)] <- FALSE
  s.unique <- paste(ped, id, sep = ":")
  if (any(duplicated(s.unique))) 
    warning("Combination of pedigree ID and ID within pedigree does not generate unique IDs")
  f.unique <- paste(ped, father, sep = ":")
  fpos <- match(f.unique, s.unique)
  m.unique <- paste(ped, mother, sep = ":")
  mpos <- match(m.unique, s.unique)
  #prv(have.snps,affected,fpos,mpos)
  
  trio <- have.snps & affected & (!is.na(fpos)) & (have.snps[fpos]) & 
    (!is.na(mpos)) & (have.snps[mpos])
  ntrio <- sum(trio, na.rm = TRUE)
  if (ntrio == 0) {
    cat("No potentially complete trios to analyse\n")
    return(NULL)
  }
  pd.snps <- in.snp[trio]
  fr.snps <- in.snp[fpos[trio]]
  mr.snps <- in.snp[mpos[trio]]
  clust <- as.integer(factor(ped[trio]))
  cord <- order(clust)
  cat("Analysing", ntrio, "potentially complete trios in", 
      max(clust), "different pedigrees\n")
  scores <- .Call("score_tdt", pd.snps[cord], fr.snps[cord], 
                  mr.snps[cord], clust[cord], snp.data, rules, snp.subset, 
                  check.inheritance, robust, uncertain, PACKAGE = "snpStats")
  chisq <- .Call("chisq_single", scores, PACKAGE = "snpStats")
  if (is.null(rules)) {
    if (is.null(snp.subset)) 
      tested <- colnames(snp.data)
    else tested <- colnames(snp.data)[snp.subset]
  }
  else {
    if (is.null(snp.subset)) 
      tested <- names(rules)
    else tested <- names(rules)[snp.subset]
  }
  if (score) 
    res <- new("SingleSnpTestsScore", snp.names = tested, 
               chisq = chisq, N = scores$N, N.r2 = scores$N.r2, 
               U = scores$U, V = scores$V)
  else res <- new("SingleSnpTests", snp.names = tested, chisq = chisq, 
                  N = scores$N, N.r2 = scores$N.r2)
  res
}
