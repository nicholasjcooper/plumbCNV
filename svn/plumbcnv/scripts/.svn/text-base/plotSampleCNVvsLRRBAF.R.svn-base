library(xtable)
library(bigmemory)
library(biganalytics)

source("/chiswick/data/ncooper/metabochipCNVanalysis2012/Scripts/DefineDirectoriesLocations.R")
source(paste(dir.scr,"FunctionsCNVAnalysis.R",sep=""))

# load main datafile
sample.fl <- "MedianNoOtherFail.txt"   #"Median.txt"  #"borderSampsOnlyMean.txt"
rngOn <- F
prepost <- 2
medSmooth <- T  # use medium smooting on the LRR
Lrr <- T   # include a plot of LRR
Baf <- F  # include a plot of BAF
One <- F  # plot just a specific id entered at cmd line (rather than list file)
One.id <- ""  # 130203_E07_1958_BC1040792  # that id ^
big.descr <- "LRRFiltSortdescrFile"
baf.descr <- "BAFFiltSortdescrFile"
out.fl <- "MedianExtremes.pdf"  #"SampleCNV_BAF_LRR.pdf"
targetChr <- 1:22  # chromosomes that the graph should be centered around
nC <- 22  # number of chromosome
med.rat <- 10  # smoothing rate for median smooth. 10=normal, lower subtle, higher very smooth
gcOverlay <- F
mbs <- 40 # megabase start position (if zooming into partial chromosome)
mbe <- 50 # megabase start position (if zooming into partial chromosome)
zoom <- F # whether to view a partial chromosome, will override other conflicting settings if necessary
customLims <- F 
L1 <- L2 <- NA
c.ylim <- ""
cust.sub <- "" 

#### COMMAND ARGS ####
# process command line options (if present)
arg.list <- parse.args(commandArgs(),coms=c("LOCAL","MDN","LRR","BAF","GC","ZOOM","PREPOST","TARG","SMOOTH","MBS","MBE","LFILE","BFILE","SFILE","OFILE","YLIM","TITLE","SAMP"),  def=c(paste(as.integer(c(rngOn,medSmooth,Lrr,Baf,gcOverlay,zoom,prepost,0,
                       med.rat))),paste(c(mbs,mbe)),
                        big.descr,baf.descr,sample.fl,out.fl,c.ylim,cust.sub,One.id))
if (!is.null(arg.list)) {
  if (arg.list["LOCAL",]=="1") { rngOn <- T } else { rngOn <- F }
  if (arg.list["MDN",]=="1") { medSmooth <- T } else { medSmooth <- F }
  if (arg.list["LRR",]=="1") { Lrr <- T } else { Lrr <- F }
  if (arg.list["BAF",]=="1") { Baf <- T } else { Baf <- F }
  if (arg.list["GC",]=="1") { gcOverlay <- T } else { gcOverlay <- F }
  if (arg.list["ZOOM",]=="1") { zoom <- T } else { zoom <- F }
  if (arg.list["PREPOST",1] %in% paste(c(0:22))) { prepost <- as.integer(paste(arg.list["PREPOST",][1])) }
  if (arg.list["TARG",1] %in% paste(c(1:22))) { targetChr <- as.integer(paste(arg.list["TARG",][1])) }
  if (arg.list["SMOOTH",1] %in% paste(c(2:1000))) { med.rat <- as.integer(paste(arg.list["SMOOTH",][1])) }
  if (!is.na(as.numeric(arg.list["MBS",1]))) { mbs <- as.numeric(paste(arg.list["MBS",][1])) }
  if (!is.na(as.numeric(arg.list["MBE",1]))) { mbe <- as.numeric(paste(arg.list["MBE",][1])) }
  if (nchar(arg.list["LFILE",])>1) { big.descr <- (paste(arg.list["LFILE",][1])) }
  if (nchar(arg.list["BFILE",])>1) { baf.descr <- (paste(arg.list["BFILE",][1])) }
  if (nchar(arg.list["SFILE",])>1) { sample.fl <- (paste(arg.list["SFILE",][1])) }
  if (nchar(arg.list["OFILE",])>1) { out.fl <- (paste(arg.list["OFILE",][1])) }
  if (length(grep(",",arg.list["YLIM",]))>0) { c.ylim <- (paste(arg.list["YLIM",][1])) }
  if (nchar(arg.list["TITLE",])>1) { cust.sub <- (paste(arg.list["TITLE",][1])) }
  if (nchar(arg.list["SAMP",])>1) { One.id <- (paste(arg.list["SAMP",][1])) }
}

######################

if(One.id!="")
{
  print("Using single sample specified in command line (rather than a list file)")
  samplelist <- One.id	
  One <- T
} else {
  samplelist <- readLines(paste(dir.ano,sample.fl,sep=""))
}

## over-ride settings incompatible with viewing a partial chromosome if ZOOM is on:
if(zoom) 
{ 
	rngOn <- T ; targetChr <- targetChr[1] 
	med.rat <- max(10,med.rat)
	prepost <- 0 ;	zxlim <- range(c(mbs,mbe))
	zxlim[zxlim<0] <- 0; zxlim[zxlim>250] <- 250  # chromosome length restrictions 
	print(paste("using zoomed plot of Chr",targetChr,"from",mbs,"mB to",mbe,"mB"))
} else { 
	zxlim <- NULL 
}


if(Lrr) { 
	bigMat2 <- getBigMat(big.descr,dir.big) ;
	print("LRR")
	print(bigMat2[1:10,"130205_H03_1958_BC1040923"])
	not.in.lrr <- !(samplelist %in% colnames(bigMat2))
	if(any(not.in.lrr))
	{
		samplelist <- samplelist[!not.in.lrr]
		print(paste("Excluding",length(which(not.in.lrr)),"samples not in the LRR matrix"))
	}
}
if(Baf) { 
	bigMat3 <- getBigMat(baf.descr,dir.big) 
	print("BAF")
	print(bigMat3[1:10,"130205_H03_1958_BC1040923"])
	print("should be different..!")
	not.in.baf <- !(samplelist %in% colnames(bigMat3))
	if(any(not.in.baf))
	{
		samplelist <- samplelist[!not.in.baf]
		print(paste("Excluding",length(which(not.in.baf)),"samples not in the BAF matrix"))
	}
}

if(c.ylim!="") {
	c.ylim <- gsub("'","",c.ylim,fixed=T)
	c.ylim <- gsub("\"","",c.ylim,fixed=T)
	l1l2 <- strsplit(c.ylim,",",fixed=T)[[1]]
	L1 <- as.numeric(l1l2[1]); 	L2 <- as.numeric(l1l2[2])
	if(!is.na(L1) & !is.na(L2)) { customLims <- T } 
}

if(cust.sub!="") {
	cust.sub <- gsub("_"," ",cust.sub,fixed=T)
}

auto.lab <- Lrr & !Baf
numplotsin1 <- sum(as.numeric(c(Lrr,Baf)))

if(!rngOn){
	targetChr <- 1:22
}

chrN <- integer(length(samplelist)); chrLab <- character(length(samplelist))
names(chrN) <- names(chrLab) <- samplelist
chrN[samplelist] <- targetChr[1]
chrLab[samplelist] <- paste("Chromosome",targetChr[1])
if(rngOn){
	chrLab[samplelist] <- ""
	 # titl <- c(paste("Sample",samplelist[cc]),paste(chrLab[cc]))
}
c.post.n <- c.pre.n <- prepost 

## plot all samples in samplelist
cnvPlotFileName <- paste(dir.qc.bor,out.fl,sep="")

CHR.INFO <- get.chr.filt.info(rownames(bigMat2),dir.ano,ret=T)

# set default to select all chromosomes for plotting
chr.select <- 1:nC

# if smoothing option in use then adjust y axis limits accordingly
if (medSmooth) { limz <- c(-.5,.5) } else { limz <- c(-2,2) }
blimz <- c(0,1)
if (gcOverlay & medSmooth) { limz <- c(-.5, .7) }
if (customLims) { limz <- c(L1,L2) }
xl <- paste(c("Genome","Chromosome")[1+as.numeric(rngOn)],"Position (Megabases)")

ofn <- paste(cnvPlotFileName,sep="")
pdf(ofn)
## loop through each bad sample plotting the desired plots
print(paste("plotting",length(samplelist),"samples..."))
for (cc in 1:length(samplelist))
{
  titl <- c(paste("Sample",deident(samplelist[cc],dir.ano)))
  if(cust.sub!="") { titl <- c(titl,cust.sub) }
  par(mfrow=c(numplotsin1,1))
  # start plotting 
  if(Lrr) {
  	  if (Baf) { xl1 <- "" ; xt <- "n" } else { xl1 <- xl; xt <- "s" }
  	  par(mar=list(c(5, 4, 10, 2),c(0.2, 4, 5, 2))[[numplotsin1]])
	  # full LRR data for chromosomes (full colour)
	  out.list <- col.plot.lrr(samplelist[cc], bigMat2, b.dir=dir.big, centre.chr=chrN[[samplelist[cc]]], 
	    CHR.INFO, plotAdj=rngOn, NULL, set.chr.ref=rngOn, c.pre.n, c.post.n, medSmooth, 
	    ratio=med.rat, ylim=limz, xlim=zxlim, main = titl,  xlab = xl1, ylab="Log-R Ratio", xaxt=xt)
	  x.coords <- out.list$x.coords; chr.select <- out.list$chr.select
	  # add chromosome labels with bad ones coloured red
	  null.result <- add.chr.top.axis(chr.select, badChrN= NULL, x.coords, nC=nC, sub=F)
	 # legend("topleft",legend=c("Raw LRR data [coloured by chromosome]"),cex=.7,
	  #         pch=c(21),col=c("grey"),bty="n")
	 if(auto.lab) { mtext ("Chromosome number",side=3,line=2.5,cex=.9) }
  }
  if (gcOverlay) {
  	# load average data for 10^6 windows prepared using 'extractGCwindowsFromSangerTxt.R'
	gc.avs6 <- reader("GCAv6.RData",dir.gc,qt=T)[[1]]/100
	gc.avs6[gc.avs6<.1] <- NA  # erase zero values as likely just due to gaps in annotation
	ofs <- out.list$offset
	xes <- c(1:length(gc.avs6))-ofs
	lines(xes,gc.avs6,col="black")
	legend("top",legend="GC Percentage",lwd=1,bty="n")
  }
  if(Baf) {
  	  par(mar=list(c(5, 4, 10, 2),c(5, 4, 0.2, 2))[[numplotsin1]])
  	  if (Lrr) { titl <- "" }
	  # full LRR data for chromosomes (full colour)
	  out.list <- col.plot.lrr(samplelist[cc], bigMat3, b.dir=dir.big, 
	    centre.chr=chrN[[samplelist[cc]]], CHR.INFO, plotAdj=rngOn, NULL, set.chr.ref=rngOn, 
        c.pre.n, c.post.n, m.smooth=F, ratio=med.rat, ylim=blimz, xlim=zxlim,
         main = titl, xlab=xl , ylab="B Allele Frequency")
	  x.coords <- out.list$x.coords; chr.select <- out.list$chr.select
	  # add chromosome labels with bad ones coloured red
	  if (!Lrr) { null.result <- add.chr.top.axis(chr.select, badChrN= NULL, x.coords, nC=nC, sub=F) }
	  #legend("topleft",legend=c("Raw BAF data [coloured by chromosome]"),cex=.7,
	   #           pch=c(21),col=c("grey"),bty="n") 
  }
  cat(".")  # progress dots
}
dev.off()
print(paste("produced file:",ofn))
