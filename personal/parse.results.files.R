source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/validation.functions.R")
library(genoset)

patch.on <- T
patch.file <- "/chiswick/data/ncooper/ImmunochipReplication/Scripts/fixDUPsPatch.RData"

# ONE-OFF to correct accidental deletion of DUPs, length =4/5
# write count correction file, to reference in all future runs from raw
# update the actual RData files each time to add the missing dups, these will
#  be counted and integrated by various things further downstream.
# will need to rerun once the full set of RUNs completes
if(F) {
  run.count.correction <- ids.to.readd <- run.samp.correction <- list()
  for (pp in c(0,1)) {
    idf <- all.runs.dup[all.runs.dup$RUN %% 2!=0 & all.runs.dup$phenotype==pp,][["id"]]
    runf <- all.runs.dup[all.runs.dup$RUN %% 2!=0 & all.runs.dup$phenotype==pp,][["RUN"]]
    summary.by.run <- tapply(idf,factor(runf),function(X) { table(table(X)) })
    run.count.correction[[pp+1]] <- sapply(summary.by.run, function(X) { sum(X[5]*5,X[4]*4) })
    run.samp.correction[[pp+1]] <- sapply(summary.by.run, function(X) { sum(X[5],X[4]) })
    #print(run.count.correction[[pp+1]])
    #print(run.samp.correction[[pp+1]])
    ids.to.readd[[pp+1]] <- tapply(idf,factor(runf),
                                   function(X) { y <- table(X); names(y)[y==4 | y==5] })
    ## now grab CNVs from these IDs/RUNs, copy them +1 to RUN
    ## also add counts to the respective phenotype tables
  }
  names(run.count.correction) <- names(ids.to.readd) <- c("0","1")
  save(run.count.correction,run.samp.correction,ids.to.readd,
       file=patch.file)
}





## read in each phenosummary and final summary and convert to a more usable format


# build paths to data files (these will have QS already)
dr <- character(2)
dr[2] <- "/chiswick/data/ncooper/immunochipRunTest/RESULTS/Runz/"         
dr[1] <- "/chiswick/data/ncooper/immunochipRunTwo/RESULTS/Runz/"
dir.ord <- rep(c(1,2,2,1,2,2,1,2,2),each=6); sf <- rep(rep(c(24,6,0),each=2),8)
all.locs <- paste(dr[dir.ord],"RUN",pad.left(1:54,"0"),"/RUN",pad.left(1:54,"0"),sf,".RData",sep="")

rdir <- "/chiswick/data/ncooper/immunochipRunTest/RESULTS/Runz/"
rdir2 <- "/chiswick/data/ncooper/immunochipRunTwo/RESULTS/Runz/"
all.r <- c(1:54)  ; rns <- pad.left(all.r,char="0")



## this is the initial loop to load in key data from each datafile
# save in the set of RUNs
DEL.thr <- .90
DUP.thr <- .90

RUNS <- vector("list",54); names(RUNS) <- rns

for(cc in 1:length(RUNS)) {
  cat(cc,"..")
  RUNS[[cc]] <- get.nxt.run(names(RUNS)[cc],rdir)
  if(length(unlist(RUNS[[cc]]))==1) {
    cat("-")
    RUNS[[cc]] <- get.nxt.run(names(RUNS)[cc],rdir2)
    if(length(unlist(RUNS[[cc]]))>1) { cat("++") } else { cat("--")}
  } else {
    cat("+")
  }
}
cat("\n")


## THIS IS TO PATCH THE FACT I CHANGED THE CODE AT SOME STAGE
# THAT CALCULATED THE SUMMARIES, and this helps make the data
# consistent for the parsing functions. This is all the latest ones
# plus a special patch for 37 which is particularly rooted, must have changed mid-run?
listo <- c(18,31:34,38:42,49:54) #'PROPER' column in excel ss
for (dd in listo) {
  if(!"samples failed due to other criteria" %in% names(RUNS[[dd]][[2]]$SAMP)) {
    RUNS[[dd]][[2]]$SAMP[["samples failed due to other criteria"]] <- 0
  }
  RUNS[[dd]][[2]]$SAMP[["samples failed Sample-QC"]] <- RUNS[[dd]][[2]]$SAMP[["samples failed due to other criteria"]]
  RUNS[[dd]][[2]]$SAMP[["samples passed pre-CNV QC"]] <- RUNS[[dd]][[2]]$SAMP[["samples passed QC"]]
  RUNS[[dd]][[2]]$SAMP <- RUNS[[dd]][[2]]$SAMP[-which(names(RUNS[[dd]][[2]]$SAMP)=="samples failed due to other criteria")] 
  RUNS[[dd]][[2]]$SAMP <- RUNS[[dd]][[2]]$SAMP[-which(names(RUNS[[dd]][[2]]$SAMP)=="samples passed QC")] 
}
RUNS[[37]][[2]]$SAMP[["samples failed Sample-QC"]] <- RUNS[[37]][[2]]$SAMP[["samples failed due to other criteria"]]
RUNS[[37]][[2]]$SAMP <- RUNS[[37]][[2]]$SAMP[-which(names(RUNS[[37]][[2]]$SAMP)=="samples failed due to other criteria")] 

#run.count.correction,ids.to.readd
load(patch.file)

## THIS LOOP IS CRUCIAL
## it loads each CNV genoset object in turn, appends to global objects, 
# used extensively later, then calculates odds ratios and denominators
# fixes sample count inconsistency between adding the CNV-QC samples to counts or not
for (cc in 1:54) {
  cat(cc,"..")
  if(file.exists(all.locs[cc])) {
    XX <- get(paste(load(all.locs[cc])))
    if(patch.on) { if(cc %% 2!=0) {
      #odd: save the extra DUPs
      rescue0 <- which(XX[[5]]$id %in% ids.to.readd[[1]][[paste(cc)]])
      rescue1 <- which(XX[[5]]$id %in% ids.to.readd[[2]][[paste(cc)]])
      allofem <- sort(c(rescue0,rescue1))
      if(length(allofem)>0) {
        rescue.seg <- XX[[5]][allofem,]
        rescue.c <- cc
      } else { rescue <- -1 }
    } else {
      #even: add in the extra DUPs
      if(rescue.c==cc-1) {
        XX[[5]] <- rbind(XX[[5]],rescue.seg)
      }
    } }
    flist <- extract.qs.qc(XX,DEL.thr=DEL.thr,DUP.thr=DUP.thr)
    XX[[4]][["RUN"]] <- cc; XX[[5]][["RUN"]] <- cc
    if(cc!=1) {
      all.runs.del <- rbind(all.runs.del,XX[[4]])
      all.runs.dup <- rbind(all.runs.dup,XX[[5]])
    } else { 
      all.runs.del <- XX[[4]]
      all.runs.dup <- XX[[5]]
    }
    rm(XX)
  } else {
    flist <- list()
  }
  if(length(unlist(RUNS[[cc]]))>1) {
    #RUNS <- update.qs.fails(RUNS,cc)
    if(cc %% 2!=0) {    RUNS <- update.samplefail.counts(RUNS,cc) }
    RUNS <- update.divisors(RUNS,cc,flist)
    cat("\n"); print(RUNS[[cc]][["divisors"]])
  }
}
       
## need to add a step to filter on quality score.
## get the number failing QS in case vs controls
## add to current summary to see whether this changes the situation

## create ways of visualising, graph for different thresholds


### NEED THIS TO GET 'Mnz' average quality scores per run, 
## otherwise most of this data helps identify which CNVs are common across
## runs, which is mainly important for the grid plots

## append quality scores to rDELs and rDUPs data, plus 'run' number or QC cats
## rbind them all together
## match which CNVs appear in each run.
## look at % of singletons in each RUN, or % below a threshold
## could rank each CNV within a RUN as to how many times it appears in other RUNs
## then enrich this by adding weight to RUNs which have more CNVs appearing more often
## could add quality scores to the weight
if(T) {
  #all.runs.del, all.runs.dup
  
  all.runs.d <- all.runs.del  #del
  #tail(sort(table(all.runs.d$id)))
  reso <- get.correspondence(all.runs.d)
  res <- reso[[1]]
  mnz <- reso[[2]]
  
if(F) {
### THESE ARE THE BIG GRID HEATMAPS
colpal2 <- colorRampPalette(c("blue","white","red"),interpolate="linear")
colz2 <- (do.call(colpal2,list(1001)))
pcs <- c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23) # the ones with PC correction

evens <- rownames(res) %in% paste(-0+((1:24)*2))
odds <- rownames(res) %in% paste(-1+((1:24)*2))
#M <- matrix(runif(10000),nrow=100) - .5
M <- res
M <- M[evens,evens]
M <- M[pcs,pcs]
pdf("testo7.pdf",width=10,height=10)
linscale <- round((round(M,3)*500)+500)
#sqscale <- round((round(sign(M)∗sqrt(abs(M)),3)∗500)+500)
colo <- colz2[linscale] #colo<colz2[sqscale]

plot(x=rep(as.numeric(colnames(M)),nrow(M)),y=rep(as.numeric(colnames(M)),each=nrow(M)),
     pch=22,bg=colo,col=NA,ylim=c(54,0),
     bty="n",ylab="rows of correlation matrix",
     xlab="columns of correlation matrix",cex=6)

#legend("top",legend=c("r1","r=+1"),bty="n",ncol=3,fill=c("blue","red"))
dev.off()
print(getwd())

}
# for DELs and DUPs
# this shows that PC 24 and 6 are quite similar, and more related
# that CNV qc is essential
# that some SNP qc is essential
# some sample QC seems good but the least effect


## IMPORT TO HAVE IN WORKSPACE FOR BOTH R-DOCs
## this imports the configuration of each run from a file:
runconfig <- read.csv("/chiswick/data/ncooper/ImmunochipReplication/Scripts/runconfigs.csv")
sel <- which(runconfig[,1] %in% paste("RUN",pad.left(names(mnz),"0"),sep=""))
runconfig <- runconfig[sel,]


## START BEST PLOT CALCs ##

# change to 1 or 2 to plot graphs for DELs/DUPs and run from here down to # ENDPLOT #
dow <- c("DEL","DUP")[1]

### setup the vector of odds ratios for each run del,dup
# vary '#' as to whether we look at the ratio considering quality score or not
vec <- numeric()
for (cc in as.numeric(names(mnz))) {
  if(dow=="DEL") {
   vec[cc] <- RUNS[[cc]][["divisors"]][7,10] # DELs with QS
  } else {
   vec[cc] <- RUNS[[cc]][["divisors"]][8,11] # DUPs with QS
  }
  #vec[cc] <- RUNS[[cc]][["divisors"]][5,8] # DELs without QS
  #vec[cc] <- RUNS[[cc]][["divisors"]][6,9] # DUPs without QS
}

if(F) {
########BEST GRAPHS########
# want to colour/shape for PC, CNV-QC, SNP-QC, Sample-QC
###PLOTS OF RUN-mean-QS versus ORs ###
#colnames(runconfig)
if(dow=="DEL") { xl <- c(0.74,1.26) } else { xl <- c(0.85,1.15) }
lwdz <- c(.7,1,2) #c(1,1,1);
pdf(paste(dow,"full_qs_v_or.pdf",sep="_"))
plot(narm(vec),mnz,pch=c(1,4,3)[runconfig[,5]+2], 
     col=c("red","blue","green")[runconfig[,2]+2],
     cex=c(1,2,3)[runconfig[,3]+2],bty="l",lwd=lwdz[runconfig[,4]+2],
     xlab="case:control CNV-ratio",ylab="whole run mean quality score") # sample-qc
#plot(narm(vec),mnz,pch=c("o","+","x")[runconfig[,5]+2], 
#     col=c("blue","orange","red")[runconfig[,3]+2]) # snp-qc
#plot(narm(vec),mnz,pch=c("o","+","x")[runconfig[,5]+2], 
#     col=c("blue","orange","red")[runconfig[,4]+2]) # pca
dev.off()
abline(h=.75,lty="dashed",col="darkgrey")
##ZOOM##
pdf(paste(dow,"zoom_qs_v_or.pdf",sep="_"))
plot(narm(vec),mnz,pch=c(1,4,3)[runconfig[,5]+2], 
     col=c("red","blue","green")[runconfig[,2]+2],ylim=c(.8,1),xlim=xl,
     cex=c(1,2,3)[runconfig[,3]+2],bty="l",lwd=c(.7,1,2)[runconfig[,4]+2],
     xlab="case:control CNV-ratio",ylab="whole run mean quality score") # sample-qc
dev.off()
### ENDPLOT ###
}
    
# plot for DELs, DUPs with and without QS
# extract sens, spec and plot similarly
# plot nice versions of the grids

# boxplots of the number of snps per run - hopefully look better once 4,5's patched in
mydata <- cbind(all.runs.del[["RUN"]],all.runs.del[["numSnps"]])
colnames(mydata) <- c("RUN","numSnps")
boxplot(numSnps~RUN,data=mydata,ylim=c(0,40))
mydata <- cbind(all.runs.dup[["RUN"]],all.runs.dup[["numSnps"]])
colnames(mydata) <- c("RUN","numSnps")
boxplot(numSnps~RUN,data=mydata,ylim=c(0,40))

    


## these aren't really necessary now we have the extraction functions
# what db?   CGHSeq = 1; bead = 9; all = 17


nsnp <- 6  # 6 8 10

## start to setup the graphs we want here:
    
do.all.plots.qs.or <- F  # can't do this until freq.dgv has been run
    
if(do.all.plots.qs.or) {
  
  ## GET DATA FOR A PARTICULAR CONDITION#
  ## PLOT FOR mnz [average quality] and 
  ##   'vec' [odds ratio] - which can be calc incl/not incl QS
#   dat <- get.dgv.result(cnv="DEL",an="SENS",db="BEAD",nsnp=6)
#   
#   plot(dat,mnz,pch=c(1,4,3)[runconfig[,5]+2], main="DEL [6 SNP]",
#        col=c("red","blue","green")[runconfig[,2]+2], 
#        cex=c(1,2,3)[runconfig[,3]+2],bty="l",lwd=lwdz[runconfig[,4]+2],
#        xlab="dgv sensitivity",ylab="whole run mean quality score") # sample-qc
#   
#   
#   plot(dat,narm(vec),pch=c(1,4,3)[runconfig[,5]+2], main="DEL [6 SNP]", 
#        col=c("red","blue","green")[runconfig[,2]+2],
#        cex=c(1,2,3)[runconfig[,3]+2],bty="l",lwd=lwdz[runconfig[,4]+2],
#        xlab="dgv sensitivity",ylab="case:control CNV-ratio") # sample-qc
  lwdz <- c(0.75,1.00,2.00)
  
  pdf("~/Documents/goodmaterials/QS.OR90.pdf") #,height=7,width=5)
  par(mfrow=c(2,1))
  for(aa in 1:length(cnv)) {
    for(bb in 1:length(db)) {
      for(cc in 1:length(sdbs)) {
        for(dd in 1:length(an)) {
          # plot QS and OR side by side
          plot.dgv.result(qs=T,cnv=cnv[aa],an=an[dd],db=db[bb],nsnp=sdbs[cc])
          plot.dgv.result(qs=F,cnv=cnv[aa],an=an[dd],db=db[bb],nsnp=sdbs[cc])
        }
      }
    }
  }
  dev.off()
}
 ##plot.dgv.result(T,cnv="DEL",an="SPEC",db="ALL",nsnp=10)


###LEGEND
plot(c(0,1,0,1),c(0.1,0.1,1,1),
     xlab="",ylab="",main="LEGEND",col="white",bty="l",xaxt="n",yaxt="n",bty="n")
legend("top",legend=c("CNV-QC off","CNV-QC on",
                      "Sample-QC off","Sample-QC relaxed","Sample-QC stringent",
                      "SNP-QC off","SNP-QC call rates + HWE","SNP-QC extensive",
                      paste("PC-correction",c("off","6 components","24 components"))),
       pch=c(1,3,15,15,15,0,0,0,rep(NA,3)),
       col=c(rep("black",2),"red","blue","green",rep("black",6)),
       pt.cex=c(1.5,1.5,1.5,1.5,1.5,.75,1.25,1.75,NA,NA,NA),
       lty=c(rep(NA,8),rep("solid",3)),
       lwd=c(rep(NA,8),0.75,1.00,2.00),bty="n")
    
    ### these are great - now at the point of using freq.dgv to check whether the low sensitivity 
    # is due to not filtering for frequency.
    # surely lots of bugs to iron out here.