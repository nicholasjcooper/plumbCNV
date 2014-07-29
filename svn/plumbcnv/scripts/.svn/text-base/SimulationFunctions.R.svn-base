
#### SIMULATION FUNCTIONS

spike.in.cnv <- function(chr,pos,sample="",dir=NULL,big.lrr=NULL,big.baf=NULL,snp.info=NULL,add=T,type=1,noise=1) {
  ## add a simulated CNV to an existing dataset (bigMatrix)
  if(!is.null(dir)) { 
    dir <- validate.dir.for(dir,c("big","ano"))
    DT <- read.data.tracker(dir)
    snp.info <- read.snp.info(dir)
    if(!is.null(DT)) {
      if(is.null(big.lrr)) { big.lrr <- "big.lrr" }
      if(is.null(big.baf)) { big.baf <- "big.baf" }
      lst <- get.tracker.for.plot(DT=DT,dir=dir,PREPOSTPC=F,baf.file=big.baf,lrr.file=big.lrr,
                                  LRR=T,BAF=T,samples=sample)
      big.lrr <- lst$lrr.file; big.baf <- lst$baf.file
    } else {
      warning("invalid data tracker or 'dir'")
    }
  } else {
    if(is.null(big.lrr) | is.null(big.baf) | is.null(snp.info)) {
      stop("if dir is empty, then big.lrr, big.baf and snp.info arguments must be provided")
    }
  }
  if(is.null(big.lrr) | is.null(big.baf)) { stop("must have valid big.lrr and big.baf big.matrix objects or 'dir' which contains these")}
  if(is(snp.info)[1]!="RangedData") { stop("snp.info must be type RangedData") }
  lrr <- get.big.matrix(big.lrr,dir$big)
  baf <- get.big.matrix(big.baf,dir$big)
  incompletes <- which( (!sample %in% colnames(lrr)) | (!sample %in% colnames(baf)) )
  if(length(incompletes)>0) {
    sample <- sample[-incompletes]    
    if(length(sample)<1) { stop("no valid samples") } else {
      warning(length(incompletes),"sample(s) removed from list as not in both baf.file and lrr.file")
    }
  }
  pos.warn <- "pos must be a vector of length 2, or a matrix with 2 columns (start, end)"
  if(!is.null(dim(pos))) {
    if(dim(pos)[2]==2) {
      if(nrow(pos)!=length(chr)) {
        chr <- rep(chr[1],times=nrow(pos)); warning("chr, pos had different lengths, only chr[1] will be used")
      }
      st <- pos[,1]; en <- pos[,2]
    } else {
      stop(pos.warn)
    }
  } else {
    if(length(pos)!=2) { stop(pos.warn) } else {
      st <- pos[1]; en <- pos[2]; chr <- chr[1]
    }
  }
  sample.num.lrr <- match(sample,colnames(lrr))
  sample.num.baf <- match(sample,colnames(baf))
  rnlrr <- rownames(lrr)
  rnbaf <- rownames(baf)
  prv.chr <- 0; rnsnp <- character()
  for (cc in 1:length(chr)) {
    rst <- start(snp.info[chr[cc]]); ren <- end(snp.info[chr[cc]])
    force.chr.pos(Pos=c(st[cc],en[cc]),Chr=chr[cc],snp.info=snp.info)
    select <- (rst >= st[cc] & ren <= en[cc])
    if(length(which(select))<1) { warning("in set #",cc,"no snps in range"); next }
    if(chr[cc]!=prv.chr) { rnsnp <- rownames(snp.info[chr[cc]]) }
    snp.list <- rnsnp[select,]
    snp.list <- snp.list[(snp.list %in% rnlrr) & (snp.list %in% rnbaf)]
    n.snps <- length(snp.list)
    baf.y <- replicate(length(sample),sim.snp.baf(n.snps,type=type,noise=noise))
    lrr.y <- replicate(length(sample),sim.snp.lrr(n.snps,type=type,noise=noise))
    snp.poz.lrr <- match(snp.list,rnlrr)
    snp.poz.baf <- match(snp.list,rnbaf)
    if(add) {
      # add on top of current LRR (e.g, keeps organic noise, assuming already normal copy)
      lrr[snp.poz.lrr,sample.num.lrr] <- lrr[snp.poz.lrr,sample.num.lrr]+lrr.y
    } else {
      # completely replace current data with simulated data
      lrr[snp.poz.lrr,sample.num.lrr] <- lrr.y
    }
    baf[snp.poz.baf,sample.num.baf] <- baf.y
  }
  flush(lrr); flush(baf)
}



#spike.in.cnv(chr=1,pos=c(1.0e+07,1.2e+07),sample="5426924068_R04C02",dir=dir,add=F,type="normal",noise=0.5)

#cnv.plot(dir=dir,samples="5426924068_R04C02",LRR=T,BAF=T,Chr=1,Pos=c(1.0e+07,1.2e+07)+c(-10^6,10^6),tag.cnv=T,lrr.file="big.lrr",Cnv=c(1.0e+07,1.2e+07))




Xcn <- function(n,m,s,f,baf=rep(NA,n),cn=3,...) {
  #simulating a DUP with cn=3/4 copies
  # simulate n points with mean m, sd s
  # f - function type, baf= beta allele freq
  # retrieve p's for each genotype based on Baf
  cn <- length(f)-1
  baf.p <- baf.ps(baf,cn=cn,...) 
  X <- numeric(n)
  for (cc in 1:n)
  {
    r <- sample(cn+1,1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


Xh <- function(n,m,s,f,baf=rep(NA,n),p.fail=0,...) {
  #simulating a copy-1-DEL or LOHz/ROHm (no hets)
  # simulate n points with mean m, sd s, p=probability of erroneous het , 
  # f - function type, baf= beta allele freq 
  baf.p <- baf.ps(baf,...) # retrieve p's for each genotype based on Baf
  baf.p[,2] <- p.fail   # replace the middle (het) value with an error probability
  X <- numeric(n)
  for (cc in 1:n)
  {
    r <- sample(3,1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


X2 <- function(n,m,s,f,baf=rep(NA,n),...) {
  # simulate a BA.Freq for each snp, then place BAF-signal accordingly
  X <- numeric(n)
  baf.p <- baf.ps(baf,...) # retrieve p's for each genotype based on Baf
  for (cc in 1:n)
  {
    r <- sample(length(f),1,prob=baf.p[cc,])
    argz <- list(m[r],s[r])
    X[cc] <- min(1,max(0,do.call(f[r],argz)))
  }
  return(X)
}


Xn <- function(n,m,s) {
  # simulate simple normal distribution
  X <- rnorm(n,m,s)
  return(X)
}


baf.ps <- function(baf=NA,cn=2,mn=.13,stdev=.15) {
  # random BAF if NA, else use the value
  n <- length(baf); sel <- which(is.na(baf)) ; m <- length(sel)
  if(m>0) {
    baf[sel] <- pmin(.5,pmax(0.0001,rnorm(m,mn,stdev)))   # m,s derived from ichip
    TF <- as.logical(round(runif(m)))
    baf[sel][TF] <- 1-baf[sel][TF]
  } 
  baf <- as.numeric(baf)
  baf <- pmin(1,pmax(0,baf))
  if(cn==3) {
    baf.p <- cbind((baf^3), 3*((baf^2)*(1-baf)), 3*(baf*((1-baf)^2)), ((1-baf)^3))
  } else {
    if(cn==4) {
      baf.p <- cbind((baf^4), 4*((baf^3)*(1-baf)), 6*((baf^2)*((1-baf)^2)), 4*(baf*((1-baf)^3)),((1-baf)^4))
    } else { 
      #i.e, cn==2#
      baf.p <- cbind(baf*baf, 2*baf*(1-baf), (1-baf)*(1-baf))
    }
  }
  return(baf.p)
}


ceils <- function(m,s) { 1-abs(rnorm(1,m,s)) }
floors <- function(m,s) { abs(rnorm(1,m,s)) }
middles <- function(m,s) { rnorm(1,m,s) }


insert.sim <- function(x=c(1:length(lrr.vec)),lrr.vec,baf.vec,type="normal",noise=1){
  if(length(lrr.vec)!=length(baf.vec)) { stop("lrr and baf need to have the same length") }
  lrr.vec[x] <- sim.snp.lrr(n=length(x),type=type,noise=noise)
  baf.vec[x] <- sim.snp.baf(n=length(x),type=type,noise=noise)
  return(list(LRR=lrr.vec,BAF=baf.vec))
}


format.cnv <- function(type) {
  if(is.character(type)) { 
    type <- tolower(type) 
    if(type=="del1") { type <- "del" }
    if(type=="dup3") { type <- "dup" }
  } else {
    if(is.numeric(type)) {
      type <- max(0,min(type,4))
    } else {
      stop("invalid type")
    }
    type <- c("del0","del","normal","dup","dup4")[type+1]
  }
  return(type)
}


sim.snp.lrr <- function(n=1,type="normal",noise=1,mode="LRR",...) {
  ## allow lower/uppercase, del/dup without number defaults to copy 1/3,
  ## numeric values can also be entered as copy number 0,1,2,3,4
  type <- format.cnv(type=type)  
  gen.fn <- switch(type,normal=normal,loh=loh,del0=del0,del=del1,dup=dup3,dup4=dup4,...)
  out <- gen.fn(n=n,noise=noise,mode=mode)
  return(out)
}


sim.snp.baf <- function(n=1,type="normal",noise=1,...) {
  return(sim.snp.lrr(n=n,type=type,noise=noise,mode="BAF",...))
}


## states for simulation

del0 <- function(n=1,noise=1,mode="LRR"){
  if(mode=="LRR") {
    m1 <- c(-3.52)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(1.32)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    out <- runif(n,0,1)
  }
  return(out)
}


del1 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- c(-0.66)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.28)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.01,.5,0.01)
    s1 <- c(.025,.5,.025)*noise
    p.fail <- .0075*noise
    f1 <- c("floors","middles","ceils")
    out <- Xh(n,m1,s1,f1,baf=baf,p.fail=p.fail,...)
  }
  return(out)
}



normal <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...) {
  if(mode=="LRR") {
    m1 <- 0 + rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- .16*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.015,.5,0.015)
    s1 <- c(.0225,.025,.0225)*noise
    f1 <- c("floors","middles","ceils")
    out <- X2(n,m1,s1,f1,baf=baf,...)
  }
  return(out)
}


loh <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- 0 + rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- .21*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0.01,.5,0.01)
    s1 <- c(.02,.5,.02)*noise
    p.fail <- .0075*noise
    f1 <- c("floors","middles","ceils")
    out <- Xh(n,m1,s1,f1,baf=baf,p.fail=p.fail,...)
  }
  return(out)
}


dup3 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...){
  if(mode=="LRR") {
    m1 <- c(0.395)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.21)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0,.35,.65,0)
    s1 <- c(.0175,.025,.025,.0175)*noise
    #p1 <- c(1,3,3,1)/8
    f1 <- c("floors","middles","middles","ceils")
    out <- Xcn(n,m1,s1,f1,baf=baf,cn=3,...)
  }
  return(out)
}




dup4 <- function(n=1,noise=1,mode="LRR",baf=rep(NA,n),...) {
  if(mode=="LRR") {
    m1 <- c(.678)+ rnorm(1,0,.002*noise) # ? SE of the plate mean
    s1 <- c(.19)*noise
    out <- Xn(n=n,m=m1,s=s1)
  } else {
    m1 <- c(0,.25,.5,.75,0)
    s1 <- c(.0175,.02,.025,.02,.0175)*noise
    #p1 <- c(1,4,6,4,1)/16
    f1 <- c("floors","middles","middles","middles","ceils")
    out <- Xcn(n,m1,s1,f1,baf=baf,cn=4,...)
  }
  return(out)
}


get.baf.lik <- function(dir="",n=100000,noise=2,recalc=F) {
  dir <- validate.dir.for(dir,"ano")
  store.path <- cat.path(dir$ano,"BAFlikelihood.RData")
  if(file.exists(store.path) & !recalc) {
    sim.bafs <- get(paste(load(store.path)))
    if(all(dim(sim.bafs)==c(101,5))) {
      return(sim.bafs)
    }
  }
  sim.bafs <- matrix(nrow=101,ncol=5); 
  rwz <- c(0:100/100)
  colnames(sim.bafs) <- paste(0:4); rownames(sim.bafs) <- paste(rwz)
  for (cc in 1:ncol(sim.bafs)) {
    sim.bafs[,cc] <- table(c(round(sim.snp.baf(n-101,type=cc-1,noise=noise),2),rwz))
  }
  likelihood.matrix <- log(sim.bafs / rowSums(sim.bafs))
  save(likelihood.matrix,file=store.path)
  return(likelihood.matrix)
}
