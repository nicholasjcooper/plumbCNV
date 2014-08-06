
get.borders.from.stats <- function(stat.dat,stat.frame,lo.stat.nm=c("LB",NA,NA,NA),
                                   hi.stat.nm=c("UB","UB","UB","UB"),lo.stat.nm2=c("-2SD",NA,NA,NA),
                                   hi.stat.nm2=c("+2SD","+2SD","+2SD","+2SD"),writeFiles=T)
{
  # CALC SAMPLES falling between possible cutoff thresholds 
  headers <- colnames(stat.dat)
  excl.list <- list()
  sample.List <- row.names(stat.frame)
  for (dd in 1:length(headers))
  {
    excl.list[[dd]] <- character()
    x <- stat.frame[,headers[dd]]
    if (!is.na(lo.stat.nm[dd])) {
      lo.stat <- stat.dat[lo.stat.nm[dd],headers[dd]]
      lo.stat2 <- stat.dat[lo.stat.nm2[dd],headers[dd]]
      excl.list[[dd]] <- c(excl.list[[dd]],sample.List[x < lo.stat2 & x > lo.stat]) }
    
    if (!is.na(hi.stat.nm[dd])) {
      hi.stat <- stat.dat[hi.stat.nm[dd],headers[dd]] 
      hi.stat2 <- stat.dat[hi.stat.nm2[dd],headers[dd]] 
      excl.list[[dd]] <- c(excl.list[[dd]],sample.List[x > hi.stat2 & x < hi.stat]) }
  }
  minn <- 4  # as regardless of how many stats, for later functions, rez needs just length 3
  fail.table <- matrix(logical(),nrow=nrow(stat.frame),ncol=min(minn,length(headers)))
  for (cc in 1:ncol(fail.table)) {
    fail.table[,cc] <- sample.List %in% excl.list[[cc]] }
  
  rownames(fail.table) <- sample.List; colnames(fail.table) <- headers
  return(fail.table)
}


rescale.dens <- function(dnsty,V,W,X=NULL,Y=NULL,ret=F,lim=T)
{
  # modify the units of a density object so the plot can be put into a
  # plot area defined by parameters V,W,X,Y
  if(is.null(X) ) { X <- range(dnsty$x) }
  if(is.null(Y) ) { Y <- range(dnsty$y) }
  dnsty$x <- ((V[2]-V[1])/(X[2]-X[1]))*(dnsty$x - X[1]) + V[1]
  dnsty$y <- ((W[2]-W[1])/(Y[2]-Y[1]))*(dnsty$y - Y[1]) + W[1]
  if (lim)
  { 
    dnsty$x[dnsty$x < V[1]] <- V[1]
    dnsty$x[dnsty$x > V[2]] <- V[2]
    dnsty$y[dnsty$y < W[1]] <- W[1]
    dnsty$y[dnsty$y > W[2]] <- W[2]
  }
  if (ret)
  { return(list(dnsty,X,Y)) } else { return(dnsty) }
}


### temporary functions to send files to myself at home
save.to.pwf <- function(...,fn) {
  save(...,file=fn)
  cmd <- paste("scp",fn,"njc72@linux.ds.cam.ac.uk:/home/njc72/")
  system(cmd)
}


send.to.pwf <- function(fn,dir="") {
  cmd <- paste("scp",cat.path(dir,fn),"njc72@linux.ds.cam.ac.uk:/home/njc72/")
  system(cmd)
}


Anova_by_Mn_Sd <- function(mns,sds,n)
{
  ## computes anova with F and pvalue without using the individual observations, just
  ##  the group means, SDs and 'n'(s). converts n to harmonic mean
  ##  if group n's are different
  harmonic_mean <- function(array) {  1/(mean((1/array),na.rm=TRUE)) }
  if (length(n)>1) { n=harmonic_mean(n) }
  k <- length(mns)
  if (length(sds)==k)
  {
    Fstat <- (n*var(mns))/(sum(sds^2)/k)
    dF <- c((k-1),((n*k)-k))
    pVal <- pf(Fstat,dF[1],dF[2],lower.tail=FALSE)
  } else {
    print(paste("means and sd vectors need to be of the same length"))
    return(NULL)
  }
  outlist <- list(Fstat,pVal,dF)
  names(outlist) <- c("F","p","df")
  return(outlist)
}





count.combs <- function(tab,grp.range=2,char="+",mode="integer")
{
  # create table of combinations for a dataframe, with 'grp.range' vars at a time
  nn <- ncol(tab)
  combs <- binary.combs(nn,colnames(tab),colchar=char)
  combs <- matrix(as(as.matrix(combs),mode),ncol=ncol(combs),
                  dimnames=list(rownames(combs),colnames(combs)))
  combs <- combs[,colSums(combs) %in% grp.range]
  # combs <- apply(combs,2,as.logical)
  return(combs)
}


binary.combs <- function(n,name=paste(c(1:n)),colchar="+")
{
  # return all possible logical combinations of a set of binary variables
  if(!is.numeric(n[1]) | n[1]<1 | n[1]>20 | round(n[1])!=n[1] | length(n)!=1) { 
    warning("n must be a positive integer <= 20") ; return(NULL) }
  nc <- 2^n
  store <- character(nc)
  for (cc in ((1:nc)-1)) 
  { 
    store[cc+1] <- (paste(rev(as.integer(intToBits(cc))), collapse="")) 
  }
  maxlen <- nchar(store[1]) #32 usually?
  firstNonZero <- min(which(strsplit(tail(store,1),"")[[1]]==1))
  store <- substr(store,firstNonZero,maxlen)
  store.fr <- as.data.frame(strsplit(store,""),stringsAsFactors=F)
  for (dd in 2:nc) {
    colnames(store.fr)[dd] <- paste(name[as.logical(as.integer(store.fr[,dd]))],collapse=colchar) }
  colnames(store.fr)[1] <- "none"
  return(store.fr)
}



circle <- function(x,y,r,plot=T,l=100,...) {
  # but can easily use graphics:symbols
  # make coords for a circle, location x,y, radius r, use length points
  l <- round(l[1]) 
  if(!is.numeric(x) | !is.numeric(y) | !is.numeric(r)) { stop("args must be numeric" ) }
  theta <- seq(0, 2 * pi, length=l)
  xo <- x + r * cos(theta)
  yo <- y + r * sin(theta)
  out <- cbind(xo,yo)
  if(plot) {
    lines(out,...)
  } else {
    return(out)
  }
}


colMedianQS <- function(mmm)
{
  # apply 'quicksort' C implemented median function to each column of a matrix
  len <- as.integer(ncol(mmm)) 
  apply(mmm,2,medianQS, len)
}


medianQS <- function(X,len=NULL,na.rm=T)
{
  # calculate median using 'quicksort' algorithm written in C
  if(na.rm) { X <- X[!is.na(X)]; len <- NULL }  
  if(is.null(len)) { len <- length(X) }
  C.out = .C("quickerselect",len,X,0.01,NAOK=T,DUP=F)
  return((C.out[[3]]))
}



search.tree2 <- function(full.line,matchlist,rng)
{
  # attempt at fast way to search through long list (not used???)
  shorter <- substr(full.line,1,rng[2])
  spli <- strsplit(shorter," ",fixed=T)[[1]][1]
  return(spli %in% matchlist)
}


create.buffer <- function(matchlist) {
  ## cool search but takes twice as long as conventional method
  char.r <- range(nchar(matchlist))
  n.snp <- length(matchlist)
  find.tree <- matrix(character(),ncol=(diff(char.r)+1),nrow=n.snp)
  dd <- char.r[1]
  find.tree[,1] <- substr(matchlist,1,dd)
  for (dd in (char.r[1]+1):char.r[2])
  {
    find.tree[,dd-char.r[1]+1] <- substr(matchlist,dd,dd)
  }
  rownames(find.tree) <- paste(1:nrow(find.tree))
  return(find.tree)
}




search.tree <- function(full.line,ftree)
{
  # attempt at fast way to search through long list (not used???)
  nf <- T; coln <- 1
  len1 <- nchar(ftree[1,1])
  next.search <- substr(full.line,1,len1)
  nC <- ncol(ftree)
  while(nf & (coln<=nC))
  {
    less.rows <- which(ftree[,coln]==next.search)
    ftree <- ftree[less.rows,]
    if(length(less.rows)==1) { 
      nf <- F
    } else { 
      next.search <- substr(full.line,len1+coln,len1+coln) 
      coln <- coln + 1
    }
  }
  if(!nf) {
    str.fnd <- (paste(ftree,collapse=""))
    nc <- nchar(str.fnd)
    if (substr(full.line,1,nc)==str.fnd)
    { return(T) } else { (return(F))  }
  } else {
    return(F)
  }
}



rg2 <- function(pos,ws,chrs)
{
  # not really sure what this little function does?
  targ.pos <- chrs+pos   #-1
  reg.num <- targ.pos %/% ws + 1
  return(reg.num)
}



check.hdr <- function(hdr,chr.list) 
{
  # for parsing headers in sanger GCBASE5 GC content UCSC files
  #if first line contains 'random' skip whole file
  ifrandom <- grep("random",hdr)
  ifchrom <- grep("chrom",hdr)
  next.chr <- 0
  if(length(ifrandom)==0) {
    if(length(ifchrom)>0) {
      hdrF <- gsub("variableStep chrom=chr","",hdr)
      hdrF <- gsub(" span=5","",hdrF)
      hdrF <- substr(hdrF,1,2)
      hdrF <- gsub("_","",hdrF)
      if (hdrF %in% paste(chr.list))
      { 
        next.chr <- as.numeric(hdrF)
      }
    }
  }
  return(next.chr)
}





rmv.dir.old <- function(X) {
  # return file name (without directory if present)
  file.segs <- strsplit(X,"/",fixed=T)[[1]]
  lss <- length(file.segs)
  if (lss>1) { out <- paste(file.segs[lss]) } else { out <- X }
  return(out)
}



chr.lab.to.num <- function(chr.col)
{
  # some annotation files will have some chromosome information entered
  # as text (e.g mitochondrial = MT, X/Y, etc) ; converts these to numbers
  chromogrps <- list("1","2","3","4","5","6","7","8","9","10","11","12"
                     ,"13","14","15","16","17","18","19","20","21","22",
                     c("X","Y","MT","0","23","24","25","26","27","28","XY"))
  chrgrp <- numeric(length(chr.col))
  for (cc in 1:length(chromogrps))
  {
    chrgrp[paste(chr.col) %in% chromogrps[[cc]]] <- cc
  }
  if (0 %in% chrgrp)
  {
    uniq.bad <- (unique(chr.col[chrgrp==0]))
    cat(paste("Error: unknown chromosome label in file:",uniq.bad,"\n"))
    cat("change file or add value to 'chromogrps' in 'chr.lab.to.num'\n")
    stop(print(chr.col[which(!paste(chr.col) %in% unlist(chromogrps))]))
  }
  return(chrgrp)
}


old.chr.list.to.txt <- function(chr.list) {
  ## takes a specific list object from some annotation lookup functions and converts to standard text positions
  nchr <- 22
  if(length(chr.list)!=nchr) { warning("less than 22 chromosomes found in list") }
  if(!is.list(chr.list)) { stop("Error: not a list") }
  text.out <- NULL
  for (cc in 1:length(chr.list))
  {
    if(length(chr.list[[cc]])>0) {
      text.out <- c(text.out,paste(names(chr.list)[cc],":",
                                   format(chr.list[[cc]]$start,scientific=F,trim=T),"-",
                                   format(chr.list[[cc]]$end,scientific=F,trim=T),sep=""))
    }
  }
  return(text.out)
}



