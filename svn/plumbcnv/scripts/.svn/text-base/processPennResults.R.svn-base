source("/chiswick/data/ncooper/metabochipCNVanalysis2012/Scripts/DefineDirectoriesLocations.R")
source(paste(dir.scr,"FunctionsCNVAnalysis.R",sep=""))

cnv.fn <- "PENNRESULTS/ALL.m2.cnv"
makeFamFile <- F  # make a plink format family file to match a .cnv file
statsOnCNVfile <- F  # produce a set of statistics on a CNV file
removeSamps <- F  # use to remove samples from a pink .cnv file (ie, with too many CNVs)
myPlinkFile <- "Main/ALL.m2.rmv.im.cent.plink"  # file to remove from (without final suffix)
ref.fn <- "Main/all.cnv"  # based on counts in
suf.for.rmv <- "TMcnvs"   # add this suffix after modifying the plink file
rare <- ""  #"DEL" or "DUP" to use specific poisson derived outlier cutoffs per sample
	
#### COMMAND ARGS ####
# process command line options (if present)
arg.list <- parse.args(commandArgs(),coms=c("FAM","STAT","RMV","FILE","TAKEFR","BASEON","SUF","RARE"),
                       def=c(paste(as.integer(c(makeFamFile,statsOnCNVfile,removeSamps))),
                        cnv.fn,myPlinkFile,ref.fn,suf.for.rmv,rare))
if (!is.null(arg.list)) {
  if (arg.list["FAM",]=="1") { makeFamFile <- T } else { makeFamFile <- F }
  if (arg.list["STAT",]=="1") { statsOnCNVfile <- T } else { statsOnCNVfile <- F }
  if (arg.list["RMV",]=="1") { removeSamps <- T } else { removeSamps <- F }
  if (nchar(arg.list["FILE",])>1) { cnv.fn <- (paste(arg.list["FILE",][1])) }
  if (nchar(arg.list["TAKEFR",])>1) { myPlinkFile <- (paste(arg.list["TAKEFR",][1])) }
  if (nchar(arg.list["BASEON",])>1) { ref.fn <- (paste(arg.list["BASEON",][1])) }
  if (nchar(arg.list["SUF",])>1) { suf.for.rmv <- (paste(arg.list["SUF",][1])) }
  if (nchar(arg.list["RARE",])>1) { rare <- (paste(arg.list["RARE",][1])) }
}
######################

cnv.fl.fn <- paste(dir.cnv,cnv.fn,sep="")
fn <- paste(dir.cnv,ref.fn,sep="")

# establish cutoffs for rare deletions/duplications using rates in all samples *2
del.rate <- 0.403; dup.rate <- 0.178; ndels <- 2023; ndups <- 895
dup.cutoff <- which(round(ppois(c(1:10),dup.rate,lower.tail=F)*ndups,4)<.05)[1]-1
del.cutoff <- which(round(ppois(c(1:10),del.rate,lower.tail=F)*ndels,4)<.05)[1]-1
print(paste("cutoff for duplicates:",dup.cutoff,"; cutoff for deletions:",del.cutoff))

do.lens.summary <- function(lenz,dat=NULL)
{
	print(paste("Length >0 and <=1kb ",length(which(lenz>0 & lenz<=1000))))
	print(paste("Length >1kb and <=10kb ",length(which(lenz>1000 & lenz<=10000))))
	print(paste("Length >10kb and <=50kb ",length(which(lenz>10000 & lenz<=50000))))
	print(paste("Length >50kb and <=100kb ",length(which(lenz>50000 & lenz<=100000))))
	print(paste("Length >100kb and <=500kb ",length(which(lenz>100000 & lenz<=500000))))
	print(paste("Length >500kb and <=5mb ",length(which(lenz>500000 & lenz<=5000000))))
	print(paste("Length >5mb ",length(which(lenz>5000000))))
	if(!is.null(dat)) {
		dat.want <- dat[,c("IID","CHR","BP1","BP2")]
		print("IDs below 1kb")
		print(cbind(dat.want,lenz)[(which(lenz>0 & lenz<=1000)),])
		print("IDs >500kb and <=5mb")
		print(cbind(dat.want,lenz)[(which(lenz>500000 & lenz<=5000000)),])
		print("IDs above 3mb")
		print(cbind(dat.want,lenz)[(which(lenz>3000000)),])
	}
}

if(makeFamFile) {
	##### create plink family file
	
	# tricky format, white space delimited, generally different numbers of spaces
	cnv.dat <- read.table(cnv.fl.fn,header=F)

	SUBINFO <- read.delim(subinfo.fn,header=T)
	
	#Family ID
	#Sample ID
	#Paternal ID
	#Maternal ID
	#Sex (1=male; 2=female; other=unknown)
	#Affection (0=unknown; 1=unaffected; 2=affected)
	
	PED.FAM.TEMP <- SUBINFO[,c(3,1,3,3,4,3)]
	PED.FAM.TEMP[,c(1,6)] <- 1
	PED.FAM.TEMP[,c(3,4)] <- 0
	ofn <- paste(dir.cnv,"plink.fam",sep="")
	write.table(PED.FAM.TEMP,file=,quote=F,col.names=F,row.names=F)
	print(paste("wrote fam file:",ofn))
}

###########

if(removeSamps)
{
	#myPlinkFile <- "Main/ALL.m2.rmv.im.cent.plink"  # file to remove from
	#fn <- paste(dir.cnv,"Main/all.cnv",sep="")  # based on counts in
	dat <- read.table(fn,header=T)  #; colnames(dat) <- 
	tt <- as.numeric(table(dat$IID))[order(as.numeric(table(dat$IID)))]
	ss <- names(table(dat$IID))[order(as.numeric(table(dat$IID)))]
	print("Mean, SD, threshold:")
	print(mean(tt)); print(sd(tt)) ; print(mean(tt)+3*sd(tt))  # mean/sd/3SD for dels/dups per subj
	thr <- round(sd(tt)*3 + mean(tt))
	thr2 <- round(sd(tt)*2 + mean(tt))
	print("plates with the most outliers (2SD)")
	print(table(substr(ss[tt>thr2],1,6))) # plates with most outliers #PLATE 7  = 130196 BAD!!
	print("mean number of CNVs in all plates")
	cnv.tab <- tapply(tt,factor(substr(ss,1,6)),mean)
	print(cnv.tab) # number of dels/dups per plate
	print("HISTO of CNVs in all plates")
	print(textogram(cnv.tab),breaks=10)
	print("CNV length stats")
	lenz <- (dat$BP2 - dat$BP1)
	print(summary(lenz)) #which(lenz==10048355)
	do.lens.summary(lenz,dat)
	ord.list <- cbind(ss,tt)[order(tt),] # ordered lists of dels/dups by SubID
	print("50 highest CNV samples")
	print(tail(ord.list,50))
	print("HISTO of CNVs per sample")
	print(textogram(as.numeric(ord.list[,2])))
	
	if (toupper(rare)=="DEL") { thr <- del.cutoff ; print("using rare del cutoff") }
	if (toupper(rare)=="DUP") { thr <- dup.cutoff ; print("using rare dup cutoff") }
	toomanyers <- ss[(which(tt>thr))]
	print(paste("removing",length(toomanyers),"samples with too many CNVs from",myPlinkFile))
	print(toomanyers)
	remove.samp.from.plink(toomanyers,paste(dir.cnv,myPlinkFile,sep=""),why=suf.for.rmv)
	print("done.")
}


if(statsOnCNVfile) {
	dat <- read.table(fn,header=T)  #; colnames(dat) <- 
	tt <- as.numeric(table(dat$IID))[order(as.numeric(table(dat$IID)))]
	# headers in file
	#FID	IID	CHR	BP1	BP2	TYPE	SCORE	SITES
	lenz <- (dat$BP2 - dat$BP1)
	print(paste("Summary of CNVs in file:",fn))
	print("Summary of CNV lengths (ALL)")
	print(summary(lenz),digits=8) #which(lenz==10048355)
	do.lens.summary(lenz,dat)
	print("log histogram") 
	textogram(log(1+lenz))
	print("Summary of CNV lengths (by chromosome)")
	print("Mean:");	print(tapply(lenz,factor(dat$CHR),mean),digits=2)
	print("Summary of CNV lengths (by cnv type)")
	print("Mean:");	print(tapply(lenz,factor(dat$TYPE),mean),digits=2)
	print("Range:"); print(tapply(lenz,factor(dat$TYPE),range),digits=2)	
	denz <- lenz/dat$SITES
	#print(cbind(lenz,dat$SITES)[lenz<20,])
	print("Summary of CNV densitys (ALL, length per SNP)")
	print(summary(denz))
	print("log histogram") 
	textogram(log(1+denz))
	print("Summary of CNV densitys (by chromosome)")
	print("Mean:");	print(tapply(denz,factor(dat$CHR),mean),digits=2)
	print("Summary of CNV densitys (by cnv type)")
	print("Mean:");	print(tapply(denz,factor(dat$TYPE),mean),digits=2)	
	print("Range:"); print(tapply(denz,factor(dat$TYPE),range),digits=2)	
	print("Summary of number of SNPS per CNV (ALL)")
	print(summary(dat$SITES)) 
	print("log histogram") 
	textogram(log(1+dat$SITES))
	print("Summary of number of SNPS per CNV (by chromosome)")
	print("Mean:");	print(tapply(dat$SITES,factor(dat$CHR),mean),digits=2)
	print("Summary of number of SNPS per CNV (by cnv type)")
	print("Mean:");	print(tapply(dat$SITES,factor(dat$TYPE),mean),digits=2)
	print("Range:"); print(tapply(dat$SITES,factor(dat$TYPE),range),digits=2)	
	
	#stats on cnvs
	print("Summary of number of CNVs per Sample")
	print(summary(tt),digits=2)
}	



