

metabo <- F
plumber <- T
dgv.valid <- F
samp.excl <- T
eval <- F

source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
#mode <- 4
#source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")

#system(paste("cp","/chiswick/data/ncooper/ImmunochipReplication/Scripts/getDataGS.sh",
#             "/home/ncooper/Documents/necessaryfilesMCHIP/getDataGS.sh"))
#system(paste("cp","/chiswick/data/ncooper/ImmunochipReplication/Scripts/getDataGS.sh",
#             "/home/ncooper/Documents/necessaryfilesICHIP/getDataGS.sh"))



#plumbCNV(dir.raw="/ipswich/data/Immunochip/FinalReports/",
#         dir.base="/chiswick/data/ncooper/ImmunochipReplication/",


auxdirI <- "/home/ncooper/Documents/necessaryfilesICHIP"
auxdirM <- "/home/ncooper/Documents/necessaryfilesMCHIP"

dir_rawM <- "/chiswick/data/store/metabochip/FinalReports/"
dir_rawI <- "/ipswich/data/Immunochip/FinalReports/"

dir_baseM <- "/chiswick/data/ncooper/metabochipRunTest/"
dir_baseI <- "/chiswick/data/ncooper/immunochipRunTwo/"

s.supM <- "/chiswick/data/store/metabochip/PLINK/Metabo_20100426_58C.bim" # support file name (eg, bim)
s.supI <- "/ipswich/data/Immunochip/FinalReports/sanger-controls.txt"

gsfM <- F
gsfI <- T

##############
## SETTINGS ##
##############

chip.settings <- list(dir.raw=dir_rawI,dir.base=dir_baseI,snp.support=s.supI,
                      gsf=gsfI,aux.files.dir=auxdirI)
#chip.settings <- list(dir.raw=dir_rawM,dir.base=dir_baseM,snp.support=s.supM,
#                     gsf=gsfM,aux.files.dir=auxdirM)


base.settings <- list(dt.name="datatracker",
                      delete.as.we.go=F,
                      grps=c(1,2,3),snp.fields=NULL,geno.file=NULL,
                      big.lrr=paste("LRR",1:3,"descrFile",sep=""),big.baf="BAFdescrFile",
                      run.mode=c("scratch","normal","big")[3],
                      snp.run.mode=c("normal","skip","plink")[1],
                      plink.imp=F,HD.mode=F,
                      n.cores=22,low.ram=F,
                      hide.penn.plink=T,penn.path="/usr/local/bin/penncnv64/",
                      ucsc="hg18",erase.previous=F,verbose=F)


penn.settings <- list(hmm="hh550.hmm",relative=T,run.manual=F,print.cmds=F,q.cores=100,
                      grid.id="all.q",cluster.fn="q.cmd",use.penn.gc=F)



if(comp) {
  snp.settings <- list(callrate.samp.thr=.94,callrate.snp.thr=.97,hwe.thr=0.0001,
                       snp.grp.miss=T,grp.hwe.z.thr=4,grp.cr.thr=.001,het.lo=.17,het.hi=.25)
} else {
  snp.settings <- list(callrate.samp.thr=.94,callrate.snp.thr=.95,hwe.thr=0.0001,
                       snp.grp.miss=F,grp.hwe.z.thr=400,grp.cr.thr=10^-20,het.lo=.05,het.hi=.50)
}




if(samp.set=="none") {
samp.settings <- list(nSD=99,mean.thr=c(NA,NA),dlrs.thr=c(NA,NA),gc.thr=c(NA,NA),
                      badPlateThresh=0.99,skip.chr.ab=T,lob=2,hib=2.5,pctile.bound=0.01,cohort.pc.correct=F,
                      batch="plate",other.batch=list(),
                      lrr.report=T,chr.ab.report=T,plate.report=T)
} else {
  if(samp.set=="light") {

    samp.settings <- list(nSD=3.5,mean.thr=c("-3.5SD","+3.5SD"),dlrs.thr=c(NA,"+3.5SD"),gc.thr=c("-3.5SD","+3.5SD"),
                      badPlateThresh=0.50,skip.chr.ab=F,lob=2.5,hib=3,pctile.bound=0.01,cohort.pc.correct=F,
                      batch="plate",other.batch=list(),
                      lrr.report=T,chr.ab.report=T,plate.report=T)
  } else {
    #heavy
    samp.settings <- list(nSD=3,mean.thr=c("LB","UB"),dlrs.thr=c(NA,"UB"),gc.thr=c("LB","UB"),
                      badPlateThresh=0.30,skip.chr.ab=F,lob=2,hib=2.5,pctile.bound=0.01,cohort.pc.correct=F,
                      batch="plate",other.batch=list(),
                      lrr.report=T,chr.ab.report=T,plate.report=T)
  }
}




if(pca.set==0) {
  pca.settings <- list(num.pcs=0,pc.to.keep=.20,assoc=F,n.store=50,correct.sex=F,
                       comparison=F,comp.gc=F,comps="plate")
} else {
  if(pca.set==6) {
    
    pca.settings <- list(num.pcs=6,pc.to.keep=.15,assoc=F,n.store=20,correct.sex=T,
                         comparison=T,comp.gc=T,comps="plate")
  } else {
    #24
    pca.settings <- list(num.pcs=24,pc.to.keep=.30,assoc=F,n.store=50,correct.sex=T,
                         comparison=T,comp.gc=T,comps="plate")
  }
}



if(!do.cnv) {
  #none
  cnv.settings <- list(out.format="Ranges",results="everything",print.summary.overlaps=T,
                     cnv.qc=F,rare.qc=F,plate.qc=F,pval=0.01,del.rate=2,dup.rate=2,thr.sd=5,plate.thr=5,
                     rmv.low.plates=F,min.sites=4,rare.olp=0.5,rare.pc=0.01,rmv.bad.reg=F)
} else {
  #normal
  cnv.settings <- list(out.format="Ranges",results="everything",print.summary.overlaps=T,
                     cnv.qc=T,rare.qc=T,plate.qc=T,pval=0.01,del.rate=0.4,dup.rate=0.2,thr.sd=3,plate.thr=3,
                     rmv.low.plates=F,min.sites=6,rare.olp=0.5,rare.pc=0.01,rmv.bad.reg=T)
}

settings <- c(chip.settings,base.settings,snp.settings,samp.settings,pca.settings,penn.settings,cnv.settings)

###################

dir <- make.dir("/chiswick/data/ncooper/immunochipRunTwo/")

if(plumber) {
  if(samp.excl) {
    ## add other sample exclusion ##
    extra.excl <- "home/ncooper/Documents/necessaryfilesICHIP/excl.samples.support.txt"
    system(paste("cp",extra.excl,dir$excl))
  }

  cnvResult <- plumbcnv(settings,start.at=my.st,pause.after=my.end,restore.mode=restore,
                      result.pref=suffix)
  save(cnvResult,file=cat.path(dir$res,"fullresult",suf=suffix,ext="RData"))
}

if(eval) {
  #source("validationWorking.R")
  do.plots <- F # don't redo fp/fn plots each time
  do.metabo <- F  # don't redo the metabochip analysis each time
  do.qs <- T
  plot.range <- 1:20
  min.snps <- 5
  print.common <- F
  ind.examples <- F
  n.pcs <- pca.set; n.pcs2 <- NA
  suffix <- suffix
  dir2 <- make.dir("/chiswick/data/ncooper/metabochipRunTest/")
  if(metabo) { source("metabochip.validation.R") }
  from.scratch <- F  # first time only set to T, otherwise F
  suffix <- suffix
  validate <- T
  dir <- make.dir("/chiswick/data/ncooper/immunochipRunTwo/")
  DT <- read.data.tracker(dir)
  thr <- .95
  n.cores <- 1
  if(dgv.valid) { source("create.dgv.validation.set.for.ichip.R") }
}
