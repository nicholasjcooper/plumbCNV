
### VERY GENERAL ###

# to make plumbcnv work, need to make this tiny fix to reader's get.delim (should fix in april)
# 
# get.delim <- function(fn,n=10,comment="#",skip=0,
#                       delims=c("\t"," ","\t| +",";",","),large=10,one.byte=TRUE)  
# {
#   # test top 'n' lines to determine what delimeter the file uses
#   if(!file.exists(fn)) { stop(paste("cannot derive delimiter as file",fn,"was not found"))}
#   test.bit <- n.readLines(fn=fn,n=n,comment=comment,skip=skip)
#   #print(test.bit)
#   num.del <- list()
#   if(any(nchar(delims)>1) & one.byte) { delims <- delims[-which(nchar(delims)>1)] }
#   for (cc in 1:length(delims)) {
#     fff <- nchar(delims[[cc]])==1
#     num.del[[cc]] <- sapply(strsplit(test.bit,delims[[cc]],fixed=fff),length)
#   }
#   #prv(num.del)
#   if(all(unlist(num.del)==1)) { 
#     warning("not a delimited file, probably a vector file")
#     return(NA)
#   }
#   # are there some delimiters that produce consistent ncol between rows?
#   need.0 <- sapply(num.del,function(X) { sum(diff(X)) })
#   num.del <- sapply(num.del,"[",1)
#   if(any(!need.0)) {
#     #rng <- range(num.del)
#     candidates <- which(num.del>1 & num.del<=large & !need.0)
#     #print(candidates)
#     if(length(candidates)>0) { out <- candidates[1] 
#     } else {
#       candidates <- which(num.del>large & !need.0)
#       if(length(candidates)>0) { out <- candidates[1]
#       } else {
#         candidates <- which(num.del==1 & !need.0)
#         if(length(candidates)>0) { out <- candidates[1]
#         } else {
#           warning("no delimiters tried were able to produce a valid file spec")
#           out <- NULL
#         }
#       }
#     }
#   } else {
#     warning("no delimiters tried were able to produce a valid file spec")
#     out <- NULL
#   }
#   #print(delims)
#   return(delims[out])
# }


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


# internal function
# normal sort function excludes NAs without asking, this keeps them in
sortna <- function(...) {
  sort(..., na.last=TRUE)
}

# iFunctions
# allows an sapply style function to only work on valid values
clean.fn <- function(x,fail=NA,fn=function(x) { x }) {
  if(!is.null(x)) { 
    x <- x[!is.na(x)]; x <- x[(x!="")]; return(fn(x)) 
  } else {  return(fail) } 
}

# convert text locations like chr5:12344-5253553 into a matrix with cols: chr, start, end
# iFunctions
convert.textpos.to.data <- function(text) {
  do.one <- function(X) {
    chr.pos <- strsplit(X,":",fixed=T)[[1]]  
    chr.txt <- gsub("chr","",chr.pos[[1]],ignore.case=T)
    chr <- chr.txt #allow for X , as.integer(chr.txt); if(is.na(chr)) { print(chr.txt) }
    pos.txt <- strsplit(chr.pos[[2]],"-",fixed=T)[[1]]
    pos12 <- as.integer(pos.txt)
    out <- c(chr,pos12[1],pos12[2])
    names(out) <- c("chr","start","end")
    return(out)
  }
  return(t(sapply(text,do.one)))
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


# iFunctions
## takes a RangedData object from some annotation lookup functions and converts to standard text positions
Ranges.to.txt <- function(chr.list) {
  if(!is(chr.list)[1] %in% c("RangedData","GRanges")) { stop("Not a GRanges or RangedData object") }
  text.out <- paste("chr",chr2(chr.list),":",format(start(chr.list),scientific=F,trim=T),"-",
                    format(end(chr.list),scientific=F,trim=T),sep="")
  return(text.out)
}

# iFunctions
# select autosomes only in RangedData object
select.autosomes <- function(snp.info,deselect=FALSE) {
  #must.use.package("genoset",bioC=T)
  typ <- is(snp.info)[1]
  if(!typ %in% c("RangedData","GRanges")) { warning("not a RangedData or GRanges object"); return(snp.info) }
  if(length(unique(chr2(snp.info))) < length(levels(chr2(snp.info)))) {
    # this fixes the problem when a subset of a ranges object with less 
    #  chromosomes still has empty chr slots from previous object
    if(typ=="RangedData") { snp.info <- snp.info[as.numeric(unique(chr2(snp.info)))] }
  } #else { cat("ok\n") }
  Chrz <- (rownames(chrInfo2(snp.info)))
  chrz <- tolower(paste(Chrz))
  if(length(grep("chr",chrz))>0) {
    select1 <- which(chrz %in% paste("chr",1:22,sep=""))
    select2 <- which(chrz %in% paste(1:22))
    if(length(select2)>length(select1)) { select <- select2 } else { select <- select1 }
  } else {
    select <- which(chrz %in% paste(1:22))
  }
  if(deselect) {
    ok.chrs <- Chrz[!select]
  } else {
    ok.chrs <- Chrz[select]
  }
  if(typ=="RangedData") {
    return(snp.info[ok.chrs])
  } else {
    return(snp.info[chr2(snp.info) %in% ok.chrs,])
  }
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


# iFunctions
# note that this will preserve only chr,start,end, nothing else including rownames
ranged.to.data.frame <- function(ranged,include.cols=FALSE,use.names=TRUE) {
  if(!include.cols) {
    u <- Ranges.to.txt(ranged)
    v <- convert.textpos.to.data(u)
    if(!is.null(rownames(ranged)) & nrow(ranged)==nrow(v) & use.names) { rownames(v) <- rownames(ranged) }
    return(v)
  } else {
    u <- as(ranged,"data.frame")
    cn <- tolower(colnames(u))
    if(is(ranged)[1]=="RangedData") {
      if("names" %in% cn) { 
        rownames(u) <- u[["names"]]
        u <- u[,-which(cn=="names")]
      } 
      if("space" %in% cn) { colnames(u)[which(cn=="space")] <- "chr" }
    } else {
      if(is(ranged)[1]=="GRanges") {
        if("seqnames" %in% cn) { colnames(u)[which(cn=="seqnames")] <- "chr" }
      } else {
        warning("'ranged' should be RangedData or GRanges, coercion could fail")
      }
    }
    return(u)
  }
}


# internal# iFunctions
chrNames2 <- function(X) {
  XX <- chrIndices2(X)
  return(rownames(XX))
}

# internal # iFunctions
# version of toGenomeOrder() that is guaranteed to work for IRanges or GRanges
toGenomeOrder2 <- function(X,...) {
  if(has.method("toGenomeOrder",X)) {
    return(toGenomeOrder(X))
  } else {
    alreadyThere <-("strand" %in% colnames(X))
    out <- toGenomeOrder(as(X,"GRanges"),strict=T)
    X <- as(out,"RangedData")
    if(("strand" %in% colnames(X)) & !alreadyThere) {
      X <- X[,-which(colnames(X) %in% "strand")]
    }
    return(X)
  }
}

#internal # iFunctions
# version of chrInfo() that is guaranteed to work for IRanges or GRanges
chrInfo2 <- function(X) {
  if(has.method("chrInfo",X)) {
    return(chrInfo(X))
  } else {
    out <- chrInfo(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chrIndices() that is guaranteed to work for IRanges or GRanges
chrIndices2 <- function(X,...) {
  if(has.method("chrIndices",X)) {
    return(chrIndices(X,...))
  } else {
    out <- chrIndices(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chr() that is guaranteed to work for IRanges or GRanges
chr2 <- function(X) {
  if(has.method("chr",X)) {
    return(chr(X))
  } else {
    if(is(X)[1]=="RangedData") {
      return(space(X))
    } else {
      warning("chr2() function applies only to RangedData objects, attempting to pass to chr()")
      return(chr(X))
    }
  }
}

#internal # iFunctions
# enter a function as a character or function
# if class is a string, will look for that name, else if an object, will search the class() of that object
has.method <- function(FUN,CLASS) {
  if(!is.character(CLASS)) { CLASS <- class(CLASS) }
  if(!is.character(FUN) & !is.function(FUN)) { stop("FUN must be an R function, as a string or function") }
  test <- showMethods(FUN,classes=CLASS,printTo=F)
  if(length(grep("not an S4 generic function",test))>0) {
    stop(FUN," was not an S4 generic function or required package not loaded")
  }
  return(!(length(grep("No methods",test))>0))
}


# iFunctions internal?
# convenience function to use GRanges
data.frame.to.granges <- function(dat,...) {
  return(data.frame.to.ranged(dat=dat,...,GRanges=TRUE))
}
 
# iFunctions
## convert any data frame with chr,start,end, or pos data into a RangedData object
# not case sensitive
data.frame.to.ranged <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,build=NULL,GRanges=FALSE) 
{
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  if((!chr %in% colnames(dat)) & ("seqnames" %in% colnames(dat)) & GRanges) { ch <- "seqnames" }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("genoset","IRanges"),T)
  g.illegal <- tolower(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))
  if(is.matrix(dat)) { dat <- as.data.frame(dat,stringsAsFactors=FALSE) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  #print(key.nms); print(colnames(dat))
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>2) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { if(!all(c(st,en,ch) %in% colnames(dat))) {
      warning("chromosome and position columns not found") } ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  ## not sure why here are adding 'chr' to X and Y?
  #this was here before? :  if(length(ch)>0) { ch1 <- gsub("Y","chrY",gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T))) } else { ch1 <- NULL }
  if(length(ch)>0) { ch1 <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(length(en1)); print(length(id)); print(length(ch1))
  outData <- GRanges(ranges=IRanges(start=st1,end=en1,names=id),seqnames=ch1); genome(outData) <- build[1]
  #outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  ###  ###  ###  outData <- toGenomeOrder2(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to re-sort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(is(outData)[1]=="GRanges") { more.cols <- more.cols[!more.cols %in% g.illegal] }
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      u <- dat[[more.cols[cc]]][reorder]; #prv(u)
      if(is(outData)[1]=="GRanges") {
        mcols(outData)[[more.cols[cc]]] <- u
      } else {
        outData[[more.cols[cc]]] <- u
      }
    }
  }
  if(GRanges) {
    return(as(outData,"GRanges"))
  } else {
    cncn <- colnames(mcols(outData))
    outData <- as(outData,"RangedData")
    if(any(cncn %in% "strand")) {
      outData <- outData[,-which(cncn=="strand")]
    }
    #outData <- toGenomeOrder2(outData,strict=T)
    return(outData)
  }
}

# iFunctions
# convert snpStats:SnpMatrix object nicely to a dataframe where coding becomes 0,1,2,NA
SnpMatrix.to.data.frame <- function(X) {
  if(!any(c("SnpMatrix","snp.matrix") %in% is(X))) { 
    warning("Invalid input type for 'X', this function is only designed ",
            "for snpStats::SnpMatrix or chopsticks::snp.matrix objects") }
  cov.data <- as.data.frame(X)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
  return(cov.data)
}

# iFunctions
# convert from a dataframe to a SnpMatrix
data.frame.to.SnpMatrix <- function(X){
  must.use.package("snpStats",T)
  if(!any(c("data.frame","matrix") %in% is(X))) { 
    warning("Invalid input type for 'X', this function is only designed ",
            "for data.frame or matrix objects") }
  if(is.matrix(X)) { X <- as.data.frame(X) }
  if(is.data.frame(X)) {
    if(any(sapply(lapply(X,is),"[",1) %in% c("character","factor"))) {
      for(cc in 1:ncol(X)) {
        X[[cc]] <- as.numeric(X[[cc]])
      }
    }
  }
  ## Test whether there are any valid values in X
  if(!all(is.na(X))) {
    mxx <- max(X,na.rm=TRUE)
    if(mxx>3) { warning("X does not appear to contain allele codes") }
    X <- round(X)
    if(mxx==3) { X <- X-1 ; X[X<0] <- NA }
  }
  NN <- as.matrix(X)
  #NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}


# specific to me
# reads snp matrices from RData binary file. if it finds multiple, will attempt to join them together
# can also handle an XSnpMatrix
get.SnpMatrix.in.file <- function(file,warn=TRUE){
  obj.nm <- paste(load(file))
  ## NOW MAKE SURE WE HAVE EXACTLY ONE SNPMATRIX OBJECT FROM THIS FILE ##
  if(length(obj.nm)>0) {
    typz <- sapply(obj.nm,function(X) { is(get(X))[1] })
    vld <- which(typz %in% c("XSnpMatrix","SnpMatrix"))
    if(length(vld)<1) { stop("no SnpMatrix objects found in file")}
    if(length(vld)>1) {
      if(warn) { warning("found multiple SnpMatrix objects in datafile, attempting to rbind() them") }
      concat.snp.matrix <- NULL; 
      vld <- vld[order(names(typz)[vld])] # alphabetical order should ensure consistency across multiple data files
      try(concat.snp.matrix <- do.call("rbind",args=lapply(obj.nm[vld],function(X) get(X))))
      if(is.null(concat.snp.matrix)) { stop("SnpMatrix objects had different numbers of SNPs [cols], could not rbind")}
      obj.nm <- "concat.snp.matrix"
    } else {
      obj.nm <- obj.nm[vld]
    }
  }
  ret.out <- get(obj.nm[1])
  if(is(ret.out)[1]=="XSnpMatrix") { if(warn) { warning("read XSnpMatrix from ",file) } }
  return(ret.out)
}

#internal # iFunctions?
chr.sel <- function(X,chr) {
  # One of the main differences between RangedData and GRanges is the way
  # of selecting the subset for a chromosome. RangedData just uses [n] where
  # 'n' is the chromosome name or number. Whereas GRanges, does not have a
  # method like this, so need to select using [chr(X)==chr.num,]
  # This wrapper allows selection of a chromosome or chromosomes regardless of
  # whether the object is RangedData or GRanges type
  typ <- is(X)[1]
  if(!typ %in% c("RangedData","GRanges")) { stop("not a GRanges or RangedData object") }
  if(!(is.character(chr) | is.numeric(chr))) { stop("chr must be character or numeric type") }
  if(is.numeric(chr)) { if(!all(chr %in% 1:99)) { 
    stop("illegal chromosome index, valid range 1-99 [although 1-28 typical for human]") } }
  if(typ=="RangedData") { return(X[chr])}
  all.chr <- chr2(X)
  if(!all(chr %in% unique(all.chr))) { 
    if(!any(chr %in% unique(all.chr))) { 
      warning("none of the specified chromosome indices were present in the GRanges object, returning NULL")
      return(NULL)
    } else { 
      warning("some of the specified chromosome indices were not present in the GRanges object") 
    }
  }
  return(X[all.chr %in% chr,])
}


### MORE GENERAL ###

# iFunctions
#' @param table.out logical, whether to return a lookup table of how names matched to integers
#' @param table.in data.frame/matrix, col 1 is the raw text names, col 2 is the integer that should be assigned,
#'  col 3 is the cleaned text (of col 1) with 'chr' removed. the required form is outputted by this function if
#'  you set 'table.out=TRUE', so the idea is that to standardize coding amongst several RangedData objects you
#'  can save the table each time and ensure future coding is consistent with this. Note that chromosomes 1-22, X,
#'  Y, XY, and MT are always allocated the same integer, so table is only useful where there are extra NT, COX, HLA
#'  regions, etc.
chrNums <- function(ranged,warn=F,table.out=F,table.in=NULL) {
  #must.use.package("genoset",bioC=T)
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) { warning("not a GRanges or RangedData object"); return(NULL) }
  lookup <- c("X","Y","XY","MT")
  txt1 <- chrNames2(ranged)
  txt <- gsub("chr","",txt1,fixed=T)
  nums <- suppressWarnings(as.numeric(txt))
  num.na <- length(nums[is.na(nums)])
  if(num.na>0) { 
    if(warn) { warning(paste("chromosome numbers requested for non-autosomes, will assign numbers >=23 to letters",
                             paste(txt[is.na(nums)],collapse=","))) }
    aux.ind <- match(txt,lookup)
    nums[!is.na(aux.ind)] <- 22+aux.ind[!is.na(aux.ind)]
    unmatched <- txt[is.na(nums)]
    if(!is.null(table.in)) {
      if((all(table.in[,1] %in% unmatched)) | (all(unmatched %in% table.in[,1]))) {
        if(all(unmatched %in% table.in[,1])) {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)] <- as.numeric(out)
        } else {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
          st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
          nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
        }
      } else {
        out <- table.in[,2][match(unmatched,table.in[,1])]
        nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
        st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
        nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
      }
    } else {
      nums[is.na(nums)] <- 27:(27+length(nums[is.na(nums)])-1)
    }
  }
  if(table.out) {
    out <- cbind(txt1,nums,txt)
    return(out)
  } else {
    return(sortna(as.numeric(nums)))
  }
}


# iFunctions
# for given genome ranges will try to find the closest snps to the end of the ranges
end.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(snp.info=snp.info,chr=chr,pos=pos,ranged=ranged,start=F,end=T,nearest=nearest))
}

# iFunctions
# for given genome ranges will try to find the closest snps to the start and end of the ranges
range.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(snp.info=snp.info,chr=chr,pos=pos,ranged=ranged,start=T,end=T,nearest=nearest))
}

# iFunctions - already in?
# get nearby number of snps?? also defined again below????
get.adj.nsnp <- function(snp.info,ranged,nsnp=10) {
  snp.info <- toGenomeOrder2(snp.info,strict=TRUE); rw.cnt <- 1
  all.fl <- matrix(ncol=4, nrow=0)
  for(cc in chrNums(ranged)) {
    nxt.nm <- rownames(snp.info[paste(cc)]); pos <- start(snp.info[paste(cc)])
    rng <- chr.sel(ranged,paste(cc)) # ranged[paste(cc)]
    st.en.snp <- range.snp(snp.info,ranged=rng)
    fl <- matrix(ncol=4, nrow=nrow(st.en.snp))
    fl[,2] <- start(rng); fl[,3] <- end(rng);
    for(dd in 1:nrow(st.en.snp)) {
      x1 <- pos[max(1,match(st.en.snp[dd,1],nxt.nm)-nsnp)]
      x2 <- pos[min(length(nxt.nm),match(st.en.snp[dd,2],nxt.nm)+nsnp)]
      #print(x1); print(x2); print(length(x1)); print(length(x2))
      fl[dd,1] <- x1[1]
      fl[dd,4] <- x2[1]
    }
    all.fl <- rbind(all.fl,fl)
  }
  fl <- (all.fl)
  fl[fl[,1]>fl[,2],1] <- fl[fl[,1]>fl[,2],2]
  fl[fl[,3]>fl[,4],1] <- fl[fl[,3]>fl[,4],4]
  return(fl)
}


# iFunctions - already in? internal?
# for given genome ranges will try to find the closest snps to the start of the ranges
start.snp <- function(snp.info,ranged=NULL,chr=NULL,pos=NULL,start=T,end=F,nearest=T) {
  # will preferably find an exact match but if nearest=T, will fall-back on nearest match
  #must.use.package("genoset",T)
  nmz <- NULL
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    if(!is.null(chr) & !is.null(pos)) {
      if(is.null(dim(pos))) { st <- pos[1]; en <- pos[2] } else {
        st <- pos[,1]; en <- pos[,2]
      }
      if(length(st)>length(chr)) { chr <- rep(chr[1],length(st)) } else { chr <- chr[1:length(st)] }
    } else {
      stop("if not using 'ranged' input, then chr and pos must be valid")
    }
  } else {
    st <- start(ranged); en <- end(ranged); chr <- chr2(ranged)
    nmz <- rownames(ranged)
  }
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    stop("snp.info must be of type RangedData or GRanges")
  } else {
    if(is.null(rownames(snp.info))) {  rownames(snp.info) <- paste(1:nrow(snp.info)) }
  }
  st.snps <- en.snps <- character(length(chr)) ; prch <- 0
  for (cc in 1:length(chr)) {
    if(chr[cc]!=prch) { 
      ref <- chr.sel(snp.info,paste(chr[cc])) # snp.info[paste(chr[cc])]
      st.ref <- start(ref); rnref <- rownames(ref)
      if(is.null(ref)) { stop(paste("snp.info did not contain chr",chr[cc])) }
    }
    #exact
    if(start) { 
      ind <- match(st[cc],st.ref)
      if(any(is.na(ind)) & nearest) {
        difs <- abs(st[cc]-st.ref)
        ind <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind)==0) { st.snps[cc] <- NA } else {
        st.snps[cc] <- rnref[ind]
      }
    }
    if(end){
      ind2 <- match(en[cc],st.ref)
      if(any(is.na(ind2)) & nearest) {
        difs <- abs(en[cc]-st.ref)
        ind2 <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind2)==0) { en.snps[cc] <- NA } else {
        en.snps[cc] <- rnref[ind2]
      }
    }
    prch <- chr[cc]
  }
  if(start & !end) {
    return(st.snps)
  }
  if(!start & end) {
    return(en.snps)
  }
  #otherwise looks like want both
  out <- cbind(st.snps,en.snps)
  if(!is.null(nmz)) { if(length(nmz)==nrow(out)) { rownames(out) <- nmz } }
  return(out)
}


# iFunctions? internal ?
# force a valid pair of positions, given the inputted chromosome
force.chr.pos <- function(Pos,Chr,snp.info=NULL,build="hg18",dir=NULL) {
  build <- ucsc.sanitizer(build)
  # convert any non autosomes to numbers:
  Chr[grep("c6",Chr,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
  Chr[grep("X",Chr,ignore.case=T)] <- 23
  Chr[grep("Y",Chr,ignore.case=T)] <- 24
  Chr[grep("M",Chr,ignore.case=T)] <- 25
  Chr[grep("NT",Chr,ignore.case=T)] <- 26  # prevent issues with NT_11387, etc
  Chr <- as.numeric(Chr)
  if(any(!paste(Chr) %in% paste(c(1:26)))) { stop("invalid chromosome(s) entered") }
  if(any(paste(Chr) == paste(26))) { warning("'NT' chromosome(s) entered, not supported, NAs produced") }
  if(length(Pos)==2 & is.numeric(Pos)) {
    if(is(snp.info)[1]!="RangedData" & is(snp.info)[1]!="GRanges") { 
      maxln <- get.chr.lens(dir=dir,mito=T,autosomes=FALSE,build=build)[Chr] 
    } else { 
      maxln <- end(tail(snp.info[paste(Chr)],1)) # force start and end to be within 1:chr.len
    }
    mbs <- min(max(1,Pos[1]),(maxln-1)); mbe <- min(max(2,Pos[2]),maxln)
    return(c(mbs,mbe))
  } else {
    Pos <- NA; warning("Pos needs to be numeric length 2, min, max")
  }
  return(round(Pos))
}



# iFunctions
# retreive GO terms from biomart for a given gene list
# can retrieve biological function, cellular component, or molecular description
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


# iFunctions
# select everything in a ranges format that is in a window of a chromosome; 
unique.in.range <- function(ranged,chr,pos,full.overlap=F, unit=c("b","kb","mb","gb"), rmv.dup=T) {
  if(length(pos)>2 | !is.numeric(pos)) { warning("pos should be a start and end numeric range"); return(NULL) }
  if(length(pos)==1) { pos <- rep(pos,2) }
  if(length(chr)>1 | is.na(as.numeric(chr))) { warning("chr should be a single number"); return(NULL) }
  if(!any(is(ranged)[1] %in% c("RangedData","IRanges","GRanges","RangesList"))) { 
    warning("'ranged' should be a RangedData type or similar"); return(NULL) }
  unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  # get set of genes in a position range for a chromosome
  chr.genez <- chr.sel(ranged,paste(chr)) 
  if(full.overlap) {
    ranged <- chr.genez[which(start(chr.genez)>min(pos) & end(chr.genez)<max(pos)),]
  } else {
    # any overlap
    ranged <- chr.genez[which(end(chr.genez)>min(pos) & start(chr.genez)<max(pos)),]
  }
  if(rmv.dup) {
    # remove duplicate genes/exons
    ranged <- ranged[!(duplicated(start(ranged)) & duplicated(end(ranged)) & duplicated(width(ranged))),]
  }
  return(ranged)
}


# iFunctions?
# add gene annotation to existing plot with y range 'ylim',
#  in range 'sect' using 'parse.annot.file' formmated
#  annotation data read in from dataframe 'DF'. step is a scaling
#  factor that will adjust the gene position units based on the 
#  the step size of the rest of the plot.
# more args to 'rect' ...
plot.gene.annot <- function(gs=NULL, chr=1, pos=NA, x.scl=10^6, y.ofs=0, width=1, txt=T, chr.pos.offset=0,
                            build="hg18", dir="", box.col="green", txt.col="black", join.col="red", ...)
{
  dir <- validate.dir.for(dir,"ano")
  build <- ucsc.sanitizer(build)
  if(is(gs)[1]!="RangedData") { gs <- get.gene.annot(dir=dir,build=build,GRanges=FALSE) }
  if(!"gene" %in% colnames(gs)) { warning("didn't find 'gene' column in annotation") ; return(NULL) }
  Col <- c("green", "darkgreen")
  if(all(is.na(pos))) { pos <- c(1,Inf) } 
  # get set of genes in range of the graph section + remove duplicate genes/exons
  rng.genez <- unique.in.range(gs,chr,pos,full.overlap=F, rmv.dup=T)
  if(nrow(rng.genez)<1) { warning("no genes found in range") ; return(NULL) }
  old.cc <- 1 ; tp <- 2
  # set vertical alignments for annotation
  old.auto <- F
  y.cent <- y.ofs
  y.bot <- y.ofs-(width/2)
  y.top <- y.ofs+(width/2)
  # text alignment
  tps <- y.bot + c(.18,.29,.46,.64,.82)[c(1,3,5)]*width
  # x position with scaling (e.g, Mb units = 10^6)
  cnrlo <- start(rng.genez)/x.scl+chr.pos.offset
  cnrhi <- end(rng.genez)/x.scl+chr.pos.offset
  gnnm <- (rng.genez$gene)
  n.genes <- length(cnrlo)
  txt.cex <- .75; if(n.genes>10) { txt.cex <- .5 } ; if(n.genes>100) { txt.cex <- .35 } # more genes = smaller labels
  for (cc in 1:n.genes) {
    # draw rectangle and label for each gene
    rect(cnrlo[cc],y.top,cnrhi[cc], y.bot,border=box.col,...)
  }
  for (cc in 1:n.genes) {
    if (gnnm[old.cc]==gnnm[cc] & cc!=1)
    {
      link <- c(cnrhi[old.cc],cnrlo[cc])
      if(link[1]<link[2]) {
        lines(link,y=rep(y.cent,2),lwd=1,col=join.col,lty="dotted")
      }
    } else {
      if(txt) {
        if(cnrlo[cc] < min(pos/x.scl)) { 
          txt.x <- mean(c(min(pos/x.scl),min(cnrhi[cc],max(pos/x.scl))),na.rm=T)  
        } else { txt.x <- cnrlo[cc] }
        text(txt.x,tps[tp],gnnm[cc],col=txt.col,cex=txt.cex,pos=4,las=2,offset=0)
      }
    }
    old.cc <- cc  
    tp <- tp+1; if(tp==4) {tp <- 1}; #if(n.genes <=5) { tp <- 2 }
  }
  return(rng.genez)
}


# iFunctions?
# should enter a ranged object with data for only 1 chromosome
# yadj controls offset from the default y-axis location of the plotted ranges (default is integers)
# adjust scl if the plot is in Mb units or Kb units (ranged object should always be in base-pairs)
# full vertical plots the ranges with abline() instead of as horizontal lines
# plot a set of ranges from a RangedData object
plot.ranges <- function(rangedat,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid",
                        full.vertical=FALSE,ylim=NULL,scl=c("b","Kb","Mb","Gb"),...) {
  #if(is(ranges)[1]!="RangedData") { warning("need RangedData object") ; return(NULL) }
  chk <- chrNums(rangedat)
  if(length(chk)>1) { warning(length(chk)," chromosomes in 'rangedat', only using the first, chr",chk[1]) ; rangedat <- rangedat[1] }
  if(is.character(scl[1])) { scltxt <- tolower(scl[1]) } else { scltxt <- "b" }
  if(scltxt %in% c("b","kb","mb","gb")) {
    scl <- 1
    if(scltxt=="kb") { scl <- 10^3 }
    if(scltxt=="mb") { scl <- 10^6 }
    if(scltxt=="gb") { scl <- 10^9 }
  } else {
    scl <- 1; warning("scale entered with illegal value, reverting to a base-pair scale")
  }
  xl <- range(c(start(rangedat),end(rangedat)))
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  nr <- nrow(rangedat); if(is.null(nr)) { nr <- length(rangedat) }
  yl <- c(0,(nr+2))
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim)
    ydif <- diff(ylim)
    yl <- ylim
  }
  YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  #print(YY)
  colz <- get.distinct.cols(nr); if(nr>22) { colz <- rep("black",nr) }
  if(is.null(labels)) { lab <- rownames(rangedat) } else { lab <- rangedat[[labels]] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) }
  if(!skip.plot.new) {
    if(full.vertical) {
      plot(x=c(start(rangedat[1,]),end(rangedat[1,]))/scl,y=YY[c(1,1)],
           xlim=xl,ylim=yl,type="l",col="white",lty=lty,...)
      abline(v=c(start(rangedat[1,]),end(rangedat[1,]))/scl,col=colz[1])
    } else {
      plot(x=c(start(rangedat[1,]),end(rangedat[1,]))/scl,y=YY[c(1,1)],
           xlim=xl,ylim=yl,type="l",col=colz[1],lty=lty,...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(full.vertical) {
        abline(v=c(start(rangedat[cc,]),end(rangedat[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        lines(x=c(start(rangedat[cc,]),end(rangedat[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc],lty=lty)
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(full.vertical) { YY <- rep(tail(YY,1),length(YY)) }
      text(x=start(rangedat[cc,])/scl,y=YY[cc]+(diff(YY[1:2])*0.5),labels=lab[cc],cex=0.6,pos=4,offset=0)
    }
  }
}



# iFunctions - internal?
#' @param keep whether to keep as object or just return the chr
set.chr.to.char <- function(ranged,do.x.y=T,keep=T) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(length(grep("chr",chrNames2(ranged)))<length(chrNames(ranged))) {
    ranged <- toGenomeOrder2(ranged,strict=TRUE)
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    RN <- rownames(ranged)
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    if(length(grep("23",paste(mychr2)))>0) { 
      warning("use of arbitrary chromosome numbers for non-autosomes (i.e, >=23)",
              "can lead to annotation issues, try to use labels, X, Y, MT, and XY where possible") }
    sar <- select.autosomes(ranged)
    if(nrow(sar)>0) {
      all.nums.t <- chrNums(sar,table.in=NULL,table.out=T) 
      all.nams <- all.nums.t[,1]
      all.nums <- all.nums.t[,2]
      #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
      for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- paste("chr",all.nums[cc],sep="") }
    } else {
      # no autosomes
    }
    if(do.x.y) {
      mychr2 <- gsub("X","chrX",mychr2)
      mychr2 <- gsub("Y","chrY",mychr2)
      mychr2 <- gsub("23","chrX",mychr2)
      mychr2 <- gsub("24","chrY",mychr2)
      mychr2 <- gsub("chrXchrY","XY",mychr2)
      mychr2 <- gsub("chrYchrX","YX",mychr2) 
      mychr2 <- gsub("MT","chrM",mychr2)
      mychr2 <- gsub("XY","chrXY",mychr2)
      mychr2 <- gsub("chrchr","chr",mychr2)
    }
    #print(tail(mychr2)); print((all.nums))
    #prv(mychr2)
    if(any(is.na(mychr2))) { prv(mychr2[which(is.na(mychr2))]) }
    if(is.null(RN) | length(RN)!=nrow(ranged)) { RN <- 1:nrow(ranged) } #make sure RN's are valid
    if(is(ranged)[1]=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),seqnames=mychr2)
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),space=mychr2)
    }
    out <- toGenomeOrder2(out,strict=TRUE)
    # prv(out)
    ###return(out)
    # need to allow for different indexing of dataframe part for GRanges
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    return(out)
  } else {
    #cat("no change\n")
    return(ranged)      # change not needed
  }
}


# iFunctions - internal?
#' @param table.in, table.out extra parameters for chrNums (e.g, how to convert weird regions)
set.chr.to.numeric <- function(ranged,keep=T,table.in=NULL,table.out=FALSE) {
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(table.out | suppressWarnings(any(is.na(as.numeric(paste(chr2(ranged))))))) {
    silly.name <- "adf89734t5b"
    ranged <- toGenomeOrder2(ranged,strict=T)
    if(typ=="GRanges") {
      mcols(ranged)[[silly.name]] <- paste(1:nrow(ranged))
    } else {
      ranged[[silly.name]] <- paste(1:nrow(ranged))
    }
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    all.nums.t <- chrNums(ranged,table.in=table.in,table.out=T) 
    all.nams <- all.nums.t[,1]
    all.nums <- all.nums.t[,2]
    #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
    for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- all.nums[cc] }
    #print(tail(mychr2)); print((all.nums))
    if(typ=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged)),seqnames=mychr2,silly.name=mcols(ranged)[[silly.name]])
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged)),space=mychr2,silly.name=ranged[[silly.name]])
    }
    out <- toGenomeOrder2(out,strict=T)
    if(typ=="GRanges") {
      oo <- mcols(out)[["silly.name"]]
      rr <- mcols(ranged)[[silly.name]]
    } else {
      oo <- out[["silly.name"]]
      rr <- ranged[[silly.name]]
    }
    if(all(!is.na(oo))) {
      if(is.null(rownames(ranged))) { rownames(ranged) <- paste(1:nrow(ranged)) ; rmv.rn <- TRUE } else { rmv.rn <- FALSE }
      #iioo <- rownames(ranged)[match(oo,rr)]; print(iioo); print(is(iioo))
      rn <- narm(rownames(ranged)[match(oo,rr)])
      if(nrow(out)==length(rn) ) { rownames(out) <- rn } else { warning("rownames did not match number of rows") }
      if(rmv.rn) { rownames(out) <- NULL }
    } else {
      warning("index column was corrupted")
    }
    # prv(out)
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% "silly.name")) { out <- out[,-which(cno %in% "silly.name")] }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% silly.name)) { out <- out[,-which(cno %in% silly.name)] }
    if(table.out) {
      return(list(ranged=out,table.out=all.nums.t))
    } else {
      return(out)
    }
  } else {
    #cat("no change\n")
    cnr <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
    if(any(cnr %in% "silly.name")) { ranged <- ranged[,-which(cnr %in% "silly.name")] }
    return(ranged)      # change not needed
  }
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
  if(is(ref)[1]=="RangedData" & all(db==c("gene","exon","dgv"))) {
    db <- "custom ranges"
  }
  db <- tolower(paste(db)[1])
  if(is(ref)[1]!="RangedData") {
    ; ano <- "gene" 
    if(length(grep("exon",db))>0) { ano <- "exon" }
    if(length(grep("dgv",db))>0) { ano <- "dgv" }
    if(!quiet) { cat(" loading",ano,"annotation...\n") }
    ref <- switch(ano,exon=get.exon.annot(dir=dir,build=build),gene=get.gene.annot(dir=dir,build=build,GRanges=FALSE),
                  dgv=get.dgv.ranges(dir=dir,build=build,compact=T)) # choose which reference to use
    
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
    pass.thresh <- jlapply(list1=mm[[3]],list2=thresh,pc=pc,op=fun) # get list of passing threshold
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

# iFunctions
## get the locations of the immunoglobin regions
get.immunog.locs <- function(build=c("hg18","hg19"),bioC=TRUE,text=FALSE) {
  nchr <- 22
  build <- ucsc.sanitizer(build[1])
  if(build[1]=="hg19") {
    #hg19
    chr <- c(22,14,2,14)
    stz <- c(22385572,105994256,89156874,22090057)
    enz <- c(23265082,107281230,89630187,23021097)
  } else {
    # hg18
    chr <- c(22,14,2,14)
    stz <- c(20715572,105065301,88937989,21159897)
    enz <- c(21595082,106352275,89411302,22090937)
  }
  nmz <- c("ig_c22","ig_c14_a","ig_c2","ig_c14_b")
  reg.dat <- rep("immunoglobin",length(chr))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chr,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- vector("list",nchr); names(outData) <- paste("chr",1:nchr,sep="")
    for (cc in 1:nchr) {
      if(cc %in% chr) {
        outData[[cc]] <- list(start=stz[chr==cc],end=enz[chr==cc])
      }
    }
  }
  return(outData)
}

# iFunctions
## get the locations of the centromeres
get.centromere.locs <- function(dir="",build=c("hg18","hg19"),bioC=TRUE,text=FALSE,autosomes=FALSE)
{
  dir <- validate.dir.for(dir,c("ano"),warn=FALSE); success <- TRUE
  build <- ucsc.sanitizer(build)
  local.file <- cat.path(dir$ano,"cyto")
  tt <- get.cyto(build=build,bioC=FALSE,dir=dir)
  chrn <- paste(1:22)
  if(!autosomes) { 
    chrn <- c(chrn,c("X","Y"))
  }
  nchr <- length(chrn)
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",chrn,sep="")
  for (cc in 1:nchr) {
    just.centros <- tt[paste(tt[,5])=="acen",]
    just.chr <- just.centros[which(paste(just.centros[,1])==names(my.chr.range)[cc]),]
    my.chr.range[[cc]] <- list(start=min(just.chr[,2]), end=max(just.chr[,3]))
  }
  reg.dat <- rep("centromere",nchr)
  nmz <- paste(reg.dat,chrn,sep="_")
  stz <- sapply(my.chr.range,"[[",1)
  enz <- sapply(my.chr.range,"[[",2)
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=TRUE)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=gsub("chr","",chrn),
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


# iFunctions
## get the locations of each cytoband (karotype)
# if dir left blank won't leave a trace
get.cyto <- function(build="hg18",dir=NULL,bioC=TRUE,refresh=FALSE) {
  build <- ucsc.sanitizer(build)
  local.file="cyto"
  if(is.null(dir)) {
    local.file <- cat.path("",local.file,suf=build,ext="tar.gz")
  } else { 
    local.file <- cat.path(dir,local.file,suf=build,ext="tar.gz")
  }
  if(!file.exists(local.file) | refresh) {
    golden.path <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",build,"/database/cytoBand.txt.gz",sep="")
    success <- tryCatch( download.file(url=golden.path,local.file,quiet=T),error=function(e) { F } )
    if(is.logical(success)) {
      if(!success) { warning("couldn't reach ucsc website! try sourcing cytoband data elsewhere"); return(NULL) } }
    tt <- reader(local.file,header=FALSE)
    if(is.null(dir)) { unlink(local.file) }
  } else {
    tt <- reader(local.file)
  }
  colnames(tt) <- c("chr","start","end","band","negpos")
  write.table(tt,file=local.file,col.names=T,row.names=F,sep="\t",quote=F)
  mychr <- gsub("chr","",tt$chr,fixed=T)
  fullbands <- paste(mychr,tt$band,sep="")
  if(bioC ) {
    st <- as.numeric(tt$start)
    en <- as.numeric(tt$end)
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=st,end=en,names=fullbands),space=mychr,
                          negpos=tt$negpos,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    #if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- tt 
    if("band" %in% colnames(outData)) {
      ## make 'chr-band' rownames to be consistent with the RangedData object if bioC=T
      rownames(outData) <- fullbands
      #outData <- outData[,-which(colnames(outData) %in% "band")]
    }
  }
  return(outData)
}


# iFunctions
#' Get HapMap recombination rates for hg18 (build 36)
#' 
#' Recombination rate files can be used to calculate recombination distances
#' for genome locations, in centimorgans. This function downloads these reference
#' files from the hapmap NCBI website. At the time of writing they were only 
#' availble for build 36. If using a more recent build I suggest using the
#' conversion function conv.37.36(), then recwindow(), then conv.36.37() to 
#' get recombination distances for other builds. If getOption("save.annot.in.current")
#' is <=0 then no files will be kept. Otherwise an object containing this mapping data
#' will be saved in the local directory if dir=NULL, or else in the directory specified.
#' Allowing this reference to be saved will greatly increase the speed of this function
#' for subsequent lookups
#' @param dir character, location to store binary file with the recombination maps for
#' chromosomes 1-22. If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param verbose logical, if the binary file is not already downloaded, when verbose
#' is TRUE, there will be some output to the console indicating the progress of the
#' download. If FALSE, all output is suppressed.
#' @param refresh logical, if you already have the binary file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param compress logical, this argument is passed to 'save' and will result in a larger
#' binary file size, but quicker loading times, so 'FALSE' is recommended for faster retrieval.
#' @export
#' @return Returns a list object of length 22, containing the recombination map files
#' as 22 separate data.frame's.
#' @example
#' ## not run as it takes roughly 2 minutes to download and read-in ##
#' ## uncomment the following 3 lines to run:
#' ## rec.map <- get.recombination.map(getwd())
#' ## file.on.disk <- "rrates_genetic_map_chr_1_22_b36.RData"
#' ## if(file.exists(file.on.disk)) { unlink(file.on.disk) } # remove the downloaded file
get.recombination.map <- function(dir=NULL,verbose=TRUE,refresh=FALSE, compress=FALSE) {
  n.chr <- 22
  hap.dir <- "http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/"
  temp.dir <- "recombinationratesGF13fDR1er119"
  local.file <- "rrates_genetic_map_chr_1_22_b36.RData"
  if(!file.exists(temp.dir)) { dir.create(temp.dir) } 
  local.files=paste0(temp.dir,"/genetic_map_chr",1:n.chr,"_b36.txt")
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    local.files <- cat.path(dir,local.files,ext="txt")
    local.file <- cat.path(dir,local.file,ext="RData")
  }
  if(!file.exists(local.file) | refresh) {
    if(verbose) { cat("Downloading recombination data from: ",hap.dir,"\n") }
    hapmap.urls <- cat.path(dir=hap.dir,fn=basename(local.files))
    success <- TRUE
    for (cc in 1:n.chr) {
      #print(hapmap.urls[cc])
      success <- tryCatch( download.file(url=hapmap.urls[cc],local.files[cc],quiet=T),error=function(e) { F } )
      if(verbose) { loop.tracker(cc,n.chr*2) }
    }
    if(is.logical(success)) {
      if(!success) { warning("couldn't download at least one of the files from: ",hap.dir); return(NULL) } }
    map.files.list <- vector("list",n.chr)
    for (cc in 1:n.chr) {
      map.files.list[[cc]] <- read.table(local.files[cc],header=TRUE)
      if(is.data.frame(map.files.list[[cc]])) { 
        unlink(local.files[cc]) 
      } else { warning("downloaded map file was corrupt for chr",cc) }
      if(verbose) { loop.tracker(n.chr+cc,n.chr*2) }
    }
    if(file.exists(temp.dir)) { file.remove(temp.dir) }   # delete the temporary directory
  } else {
    map.files.list <- reader(local.file)
  }
  if(!is.null(dir)) { save(map.files.list,file=local.file,compress=compress) }
  if(length(map.files.list)!=n.chr) { stop("Unfortunately the object derived seems corrupted") }
  names(map.files.list) <- paste0("chr",1:n.chr)
  return(map.files.list)
}


# iFunctions
## get list of all exons, names, starts, ends
get.exon.annot <- function(dir=NULL,build="hg18",bioC=T,list.out=F) {
  build <- ucsc.sanitizer(build)
  ## load exon annotation (store locally if not there already)
  from.scr <- T
  if(!is.null(dir)) {
    ex.fn <- cat.path(dir$ano,"exonAnnot.RData")
    if(file.exists(ex.fn)) {
      ts <- get(paste(load(ex.fn)))
      if((bioC|!list.out) & is.data.frame(ts)) { from.scr <- F }
      if((list.out & !bioC) & is.list(ts)) { from.scr <- F }
    }
  } 
  must.use.package("GenomicFeatures",T)
  if(from.scr) {
    must.use.package("gage",T)
    # get transcripts from build table 'knownGene'
    success <- tryCatch(txdb <- makeTranscriptDbFromUCSC(genome=build,
                                                         tablename="knownGene")  ,error=function(e) { F } )
    if(is.logical(success)) { 
      if(!success) {
        warning("Couldn't reach build website! try again later or, \n",
                "if in europe/uk, there may still be a bug in rtracklayer; \n",
                "Installing the latest version of R and bioconductor\n",
                "and running biocLite('rtracklayer'), should fix this")
        return(NULL) }
    }
    ts = transcriptsBy(txdb, by="gene")
    data(egSymb)
    select <- match(names(ts),egSymb[,1])
    names(ts)[!is.na(select)] <- egSymb[,2][select[!is.na(select)]]
  }
  if(!list.out | bioC) {
    if(from.scr) {
      ts <- as.data.frame(ts)
      chrs <- paste(ts$seqnames); chrs <- gsub("chr","",chrs)
      ts$seqnames <- chrs
      if(all(colnames(ts)==c("element","seqnames","start","end","width","strand","tx_id","tx_name"))) {
        colnames(ts) <- c("gene","chr","start","end","width","strand","txid","txname")
      } else {
        cat(" unexpected colnames found using makeTranscriptDbFrombuild()\n")
        if(bioC) { cat(" therefore returning data.frame instead of RangedData object\n") ;
                        bioC <- F }
      }
      if(exists("ex.fn")) { save(ts,file=ex.fn) }
    }
    if(bioC) {
      ts <- RangedData(ranges=IRanges(start=ts$start,end=ts$end),
                       space=ts$chr,gene=ts$gene, strand=ts$strand,
                       txid=ts$txid, txname=ts$txname,universe=build)
      ts <- toGenomeOrder2(ts,strict=T)
    }
  } else {
    if(from.scr) { if(exists("ex.fn")) { save(ts,file=ex.fn) } }
  }
  return(ts)
}

# iFunctions
# internal, tidy chromosome names using extra chromosomal annotation into rough chromosomes
tidy.extra.chr <- function(chr,select=FALSE) {
  # most relevant to hg18
  SEL_c6 <- grep("c6",chr,ignore.case=T)
  SEL_c5 <- grep("c5",chr,ignore.case=T)
  SEL_NT <- grep("NT",chr,ignore.case=T)
  # most relevant to hg19
  SEL_LRG <- grep("LRG",chr,ignore.case=T)
  SEL_HG <- grep("HG",chr,ignore.case=T)
  SEL_GL <- grep("GL",chr,ignore.case=T)
  SEL_HS <- grep("HSCHR",chr,ignore.case=T)
  if(select) {
    # create TRUE/FALSE as to whether list elements have weird chromosome codes
    all <- unique(c(SEL_c6,SEL_c5,SEL_NT,SEL_LRG,SEL_HG,SEL_GL,SEL_HS))
    return(!1:length(chr) %in% all)
  } else {
    # transform weird chromosomes into more palatable codes
    chr[SEL_c6] <- 6  # prevent issues with c6_COX, c6_QBL  
    chr[SEL_c5] <- 5  # prevent issues with c5_H2  
    chr[SEL_NT] <- "Z_NT"  # merge all NT regions to one label
    chr[SEL_LRG] <- "Z_LRG"  # merge all NT regions to one label
    chr[SEL_HG] <- "Z_HG"  # merge all NT regions to one label
    chr[SEL_GL] <- "Z_GL"  # merge all NT regions to one label
    X <- names(table(chr))
    X <- X[grep("HSCHR",X)]
    if(length(X)>0) {
      HSC <- gsub("_","",substr(gsub("HSCHR","",X),1,2))
      for(cc in 1:length(X)) {
        #cat("replacing ",X[cc]," with ",HSC[cc],"\n",sep="")
        chr[chr==X[cc]] <- HSC[cc]
      }
    }
    return(chr)
  }  
}

# iFunctions
## get list of all genes, names, starts, ends, optionally bands
#' @param ens.id logical, whether to include the ensembl id in the dataframe
get.gene.annot <- function(dir=NULL,build="hg18",bioC=TRUE,duplicate.report=FALSE,
                  one.to.one=FALSE,remap.extra=FALSE,discard.extra=TRUE,ens.id=FALSE,
                           refresh=FALSE,GRanges=TRUE) {
  # faster than exon, but only contains whole gene ranges, not transcripts
  # allows report on duplicates as some might be confused as to why some genes
  # have more than one row in the listing (split across ranges usually)
  # run with dir as NULL to refresh changes in COX
  must.use.package(c("biomaRt","genoset","gage"),T)
  build <- ucsc.sanitizer(build)
  from.scr <- T
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    utxt <- ""; if(one.to.one) { utxt <- "_unq" }
    if(ens.id) { utxt <- paste(utxt,"ens",sep="_") }
    gn.fn <- cat.path(dir$ano,"geneAnnot",pref=build,suf=utxt,ext="RData")
    if(file.exists(gn.fn) & !refresh) {
      dat <- get(paste(load(gn.fn)))
      from.scr <- F
    }
  }
  # colnames for output
  nm.list <- c("gene","chr","start","end","band")
  if(ens.id & bioC) { warning("ens.id=TRUE only has an effect when bioC=FALSE") }
  if(from.scr) {
    if(build=="hg18") {
      ens <- useMart("ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     host="may2009.archive.ensembl.org",
                     path="/biomart/martservice",
                     archive=FALSE)
    } else {
      ens <- useMart("ensembl")
    }
    ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
    #data(egSymb)
    attr.list <- c("hgnc_symbol", "chromosome_name",
      "start_position", "end_position", "band")
    if(ens.id) { attr.list <- c(attr.list,"ensembl_gene_id") }
    #dat <- getBM(attributes = attr.list, filters = "hgnc_symbol",
    #             values = egSymb[,2], mart = ens)
    dat <- getBM(attributes = attr.list,  mart = ens)
    if(exists("gn.fn")) { save(dat,file=gn.fn) }
  } 
  if(ens.id) { nm.list <- c(nm.list,"ens.id") }
  if(remap.extra) {
    dat$chromosome_name <- tidy.extra.chr(dat$chromosome_name)
    #dat$chromosome_name[grep("c6",dat$chromosome_name,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
    #dat$chromosome_name[grep("c5",dat$chromosome_name,ignore.case=T)] <- 5  # prevent issues with c5_H2  
    #dat$chromosome_name[grep("NT",dat$chromosome_name,ignore.case=T)] <- "Z_NT"  # merge all NT regions to one label
  }
  if(discard.extra) {
    ## http://www.lrg-sequence.org/ ##
    # note that if remapping is already done, then these won't be discarded unless remapping failed
    tt <- tidy.extra.chr(dat$chromosome_name,select=TRUE)
    badz <- which(!tt)
    if(length(badz)>0) { dat <- dat[-badz,] } # remove LRG, GS, HG, NT, COX, etc annotation from set
  }
  if(bioC) {
    outData <- RangedData(ranges=IRanges(start=dat$start_position,end=dat$end_position),
                          space=dat$chromosome_name,gene=dat$hgnc_symbol, band=dat$band, universe=build)
    outData <- toGenomeOrder2(outData,strict=T)
    if(duplicate.report) {
      dup.genes <- gene.duplicate.report(outData)
      # but haven't done anything about them or removed them! 
    }
    if(GRanges) { outData <- as(outData, "GRanges") }
  } else {
    outData <- dat; colnames(outData) <- nm.list
  }
  return(outData)
}


# iFunctions
## create telomere locations - artibrary number of kb at start and end of each CHR
get.telomere.locs <- function(dir="",kb=10,build=c("hg18","hg19"),bioC=TRUE,text=FALSE,autosomes=FALSE,mito.zeros=FALSE)
{
  # the actual telomeres are typically about 10kb, but
  # for cnv-QC purposes want to exclude a larger region like 500kb
  # Mt have no telomeres, are circular, but for some purposes might want zero values in there
  build <- ucsc.sanitizer(build)
  chr.lens <- get.chr.lens(dir=dir,build=build[1],autosomes=FALSE,mito=mito.zeros)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y") } # Mt have no telomeres, are circular
  nchr <- length(n) #default
  if(mito.zeros) { n <- c(n,"M") }
  my.chr.range <- vector("list",length(n))
  names(my.chr.range) <- paste("chr",n,sep="")
  for (cc in 1:nchr) {
    one <- force.chr.pos(Pos=c(1,kb*1000),Chr=cc,build=build) # f..c..pos() makes sure is a valid range
    two <- force.chr.pos(Pos=chr.lens[cc]+c(-kb*1000,0),Chr=cc,build=build)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  if(mito.zeros) {
    # add null values for the Mitochondrial chromosome
    cc <- cc+1; one <- c(1,1);  two <- chr.lens[cc]+c(0,0)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  reg.dat <- rep("telomere",length(n)*2)
  chrz <- rep(n,each=2)
  nmz <- paste(reg.dat,chrz,rep(c("a","b"),times=length(n)),sep="_")
  stz <- as.vector(sapply(my.chr.range,"[[",1))
  enz <- as.vector(sapply(my.chr.range,"[[",2))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chrz,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


# iFunctions
#' Get chromosome lengths from build database
#' 
#' Quick and easy way to retrieve human chromosome lengths. Can select from hg18/hg19 (ie, 
#'  build 36/37), or any future builds (hg20, etc) stored in the same location on the build website.
#'  Default is to return lengths for 22 autosomes, but can also retrieve X,Y 
#'  and Mitochondrial DNA lengths by 'autosomes=FALSE' or n=1:25. Even if not connected to 
#'  the internet can retrieve hard coded lengths for hg18/hg19.
#'  @param dir directory to retrieve/download the annotation from/to (defaults to current getwd())
#'  if dir is NULL then will automatically delete the annotation text file from the local directory
#'   after downloading
#'  @param build string, currently 'hg17','hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 17,18,19,35,36,37 as alternative arguments.
#'  @param autosomes logical, if TRUE, only load the lengths for the 22 autosomes, else load X,Y,[MT] as well
#'  @param len.fn optional file name to keep the lengths in
#'  @param whether to include the length of the mitochondrial DNA (will not include unless autosomes is also FALSE)
#'  @param delete.after logical, if TRUE then delete the text file that these lengths were downloaded to. 
#'  If FALSE, then the file will be kept, meaning future lookups will be faster, and available offline.
#'  @examples
#'  get.chr.lens(delete.after=TRUE) # delete.after simply deletes the downloaded txt file after reading
#'  get.chr.lens(build=35,autosomes=TRUE,delete.after=TRUE) # only for autosomes
#'  get.chr.lens(build="hg19",mito=TRUE,delete.after=TRUE) # include mitochondrial DNA length
get.chr.lens <- function(dir="",build=c("hg18","hg19")[1],autosomes=FALSE,len.fn="humanChrLens.txt",
                         mito=FALSE,delete.after=FALSE, verbose=FALSE)
{
  # retrieve chromosome lengths from local annotation file, else download from build
  if(is.null(dir)) { dir <- getwd() ; delete.after <- TRUE }
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  chrlens.f <- cat.path(dir$ano,len.fn) # existing or future lengths file
  build <- ucsc.sanitizer(build)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y","M") }
  hg18.backup <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
                      146274826,140273252,135374737,134452384,132349534,114142980,106368585,
                      100338915,88827254,78774742,76117153,63811651,62435964,46944323,
                      49691432,154913754,57772954,16571)
  hg19.backup <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                      146364022,141213431,135534747,135006516,133851895,115169878,107349540,
                      102531392,90354753,81195210,78077248,59128983,63025520,48129895,
                      51304566,155270560,59373566,16571)
  
  # backups for offline use
  if(build=="hg18") { offline.backup <- hg18.backup  } else { offline.backup <- hg19.backup }
  if(file.exists(chrlens.f))
  {
    # file seems to be in annotation directory already
    chrLens <- readLines(chrlens.f)
    if (length(chrLens)!=length(n))
    {
      #warning("Length of existing chromosome file didn't match expected:",length(n))
      notGot <- T
    } else {
      notGot <- F
      # we have the right length, but do we have the right version?
      if(build=="hg18" & length(which(chrLens %in% hg19.backup))>2) { notGot <- T }
      if(build=="hg19" & length(which(chrLens %in% hg18.backup))>2) { notGot <- T }
    }
  } else { notGot <- T }
  if (notGot | (!build %in% c("hg18","hg19"))) {
    #download from build
    if(verbose) { cat("attempting to download chromosome lengths from genome build ... ") }
    urL <- switch(build,
                  hg17="http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/chromInfo.txt.gz",
                  hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
                  hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",
                  hg20="http://hgdownload.cse.ucsc.edu/goldenPath/hg20/database/chromInfo.txt.gz")
    success <- T
    success <- tryCatch(download.file(urL, chrlens.f,quiet=T),error=function(e) { F } )
    if(!is.logical(success)) { success <- T }
    if(success) {
      if(verbose) {  cat("download successful\n") }
      chrL.f <- readLines(chrlens.f)
      len.lst <- strsplit(chrL.f,"\t")
      nmz <- sapply(len.lst,"[",1)
      lnz <- sapply(len.lst,"[",2)
      want.chr.names <- match(paste("chr",n,sep=""),nmz)
      want.chr.names <- want.chr.names[!is.na(want.chr.names)]
      chrLens <- lnz[want.chr.names]
      names(chrLens) <- want.chr.names
    } else {
      warning("couldn't reach build website, so have used offline versions of chr lengths")
      if(!build %in% c("hg18","hg19")) { warning("no offline version for build version:",build) }
      chrLens <- paste(offline.backup)[1:length(n)]; names(chrLens) <- n
      delete.after <- F
    }
    if(length(dir)==1 & dir[1]=="" & delete.after) {
      unlink(chrlens.f)
    } else {
      writeLines(chrLens,con=chrlens.f) # save file for future use
    }
  }
  if(!mito & length(chrLens)>22) { chrLens <- chrLens[-grep("M",n)] }
  return(as.numeric(chrLens))
}


# iFunctions?
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
chip.coverage <- function(snp.info,targ.int=100000,by.chr=F,min.snps=10,
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

# internal function used by chip.coverage() to calculate coverage gaps
calc.cov <- function (snp.info, targ.int, min.snps) {
  cur.chr <- chr2(snp.info)
  if(nrow(snp.info)<10) { warning("not enough snps to calculate coverage"); return(NA) }
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
        while(dif[cc]<gap) {
          # see how many snps we need to skip ahead until the gap to 10th [min.snp] snp is <= targ.int
          # record this distance as a gap where a CNV can't be detected, then start again from the 'skipnum' snp
          skipnum <- skipnum+1
          dif[cc] <- startz[cc+skipnum] - startz[cc]   ## old idea, don't think it works : dif2 <- startz[cc+min.snps+skipnum+1] - startz[min.snps+cc]
        }
        beg[cc] <- startz[cc]; fin[cc] <- startz[cc+skipnum]; Chr[cc] <- cur.chr[cc]
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

# import the DGV into a rangedData object
get.dgv.ranges <- function(dir=NULL,build="hg18",bioC=TRUE,text=FALSE,shortenID=TRUE, compact=FALSE, alt.url=NULL)
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
  return(outData)
}



# do abline at start and end of all CNV positions
# ie to make a graph of a CNV easier to read
draw.cnv.bounds <- function(cnv,chr.offset=NA,pos=NA,cnv.lty="dashed",cnv.col="orange",plot.scl=10^6) {
  done <- F
  if(is(cnv)[1]=="RangedData") {
    # highlight CNVs in rangedData object
    if(all(is.na(chr)) & all(is.na(pos))) {
      cnv <- unique.in.range(cnv,chr,pos,full.overlap=F, rmv.dup=T)
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

# import a ped file to pedData format used by snpStats
get.pedData <- function(file,correct.codes=TRUE,silent=FALSE) {
  rr <- read.ped.file(file,keepsix = TRUE)
  want <- c("familyid","individual","father","mother","sex","affected")
  #ignore <- c("sex","gender","mf")
  #to.ignore <- tolower(colnames(rr)) %in% ignore
  #if(any(to.ignore)) {
  #  rr <- rr[,-which(to.ignore)]
  #}
  have <- colnames(rr)
  if(!silent) { cat(paste("mapping column",have,"==>",want,"\n"),sep="") }
  if(ncol(rr)>6) { warning("expected 6 columns in pedfile, unexpected behaviour could result") }
  if(ncol(rr)<6) { stop("need at least 6 columns in pedfile to map to headings ",paste(want,collapse="")) }
  colnames(rr)[1:length(want)] <- want
  rr <- rr[order(rr$familyid),]
  rr[["member"]] <- unlist(tapply(rep(1,nrow(rr)),factor(rr$familyid),cumsum))
  rr <- rr[,c(2,1,7,3:6)]
  rr <- shift.rownames(rr,T)
  rr[["father"]] <- rr$member[match(rr$father,rownames(rr))]
  rr[["mother"]] <- rr$member[match(rr$mother,rownames(rr))]
  if(correct.codes) {
    badz <- with(rr,which(father==1 & mother==1))
    rr[badz,"mother"] <- 2
    badz <- with(rr,which(father==2 & mother==2))
    rr[badz,"father"] <- 1
  }
  #print(head(rr))
  return(rr)
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
