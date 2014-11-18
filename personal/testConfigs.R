if(!exists("do.qual")) { do.qual <- TRUE }

source("~/github/NCmisc/NCmisc.R")

source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

ped.file <- "~/Documents/necessaryfilesICHIPFam/t1dgc-pedfile-2011-08-05.tab"

(load(paste0("RESULTS/TDT_results",suffix,".RData"))) # these don't seem to have the right data...

dir <- make.dir(getwd())

DT <- read.data.tracker(dir)
DT$CnvResults$RESULT[[1]] <- paste0("RESULTS/",suffix,"24.RData")
if(do.qual) {
  qs.results <- process.quality.scores(DT,suffix=suffix,dir,n.pcs=NA,restore=do.qual)
}

cnvResult <- reader(paste0("RESULTS/fullresult",suffix,".RData"))

final.tables <- compile.qs.results.to.cc.toptable(qs.results,dir=dir,suffix=suffix,cnvResult)

outlist <- annotate.denovos(qs.results=qs.results,sum.del=sum.del,tdt3=tdt3,ped.file=ped.file)
#list.to.env(outlist)

DELK <- outlist[["kids"]]
DELP <- outlist[["parents"]]


filt <- rep(TRUE, nrow(DELK))
cat("filt = ALL\n")
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1]))

cat("filt = rm decline is NA\n")
filt <- with(DELK,!is.na(decline) & abs(decline)<.99)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm decline <.5\n")
filt <- with(DELK,!is.na(decline) & abs(decline)<.5)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm decline <.2\n")
filt <- with(DELK,!is.na(decline) & abs(decline)<.2)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm decline <.1\n")
filt <- with(DELK,!is.na(decline) & abs(decline)<.1)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm QS>.95\n")
filt <- with(DELK,score>.95)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm QS>.75\n")
filt <- with(DELK,score>.75)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))

cat("filt = rm QS>.5\n")
filt <- with(DELK,score>.75)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))



cat("filt = rm QS>.0001\n")
filt <- with(DELK,score>.0001)
print(FETable(table(DELK$denovo[filt],DELK$phenotype[filt])[,2:1],Ns=c(649,165)))
