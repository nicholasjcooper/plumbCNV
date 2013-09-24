source("FunctionsCNVAnalysis.R")
load.all.libs()
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/validation.functions.R")

thr <- "00"
sites <- 6
datab <- "CGHSEQ"
DELs <- F
ord <- 6   # 8 = SENS, 12= SPEC, 6 =QS, 7 = OR


load(paste("/chiswick/data/ncooper/ImmunochipReplication/Scripts/SaveThr",thr,".RData",sep=""))

vec2 <- numeric()
for (cc in as.numeric(names(mnz))) {
    vec2[cc] <- RUNS[[cc]][["divisors"]][8,11] # DUPs with QS
  #vec[cc] <- RUNS[[cc]][["divisors"]][5,8] # DELs without QS
  #vec[cc] <- RUNS[[cc]][["divisors"]][6,9] # DUPs without QS
}

if(DELs) {
  dat1 <- get.dgv.result(an="SENS",cnv="DEL",db=datab,nsnp=sites,stat="COUNT")
  dat2 <- get.dgv.result(an="SENS",cnv="DEL",db=datab,nsnp=sites,stat="PC")
  dat3 <- get.dgv.result(an="SPEC",cnv="DEL",db=datab,nsnp=sites,stat="COUNT")
  dat4 <- get.dgv.result(an="SPEC",cnv="DEL",db=datab,nsnp=sites,stat="PC")
} else {  
  dat1 <- get.dgv.result(an="SENS",cnv="DUP",db=datab,nsnp=sites,stat="COUNT")
  dat2 <- get.dgv.result(an="SENS",cnv="DUP",db=datab,nsnp=sites,stat="PC")
  dat3 <- get.dgv.result(an="SPEC",cnv="DUP",db=datab,nsnp=sites,stat="COUNT")
  dat4 <- get.dgv.result(an="SPEC",cnv="DUP",db=datab,nsnp=sites,stat="PC")
}
all.runs.d <- all.runs.dup
mnz2 <- tapply(all.runs.d$score,factor(all.runs.d$RUN),function(x) { round(mean(x,na.rm=T),3) })

#write.csv(cbind(runconfig,mnz,vec,dat1,dat2,mnz2,vec2,dat3,dat4),file="mysummary.csv")
if(DELs) { mnmn <- mnz; vecz <- vec } else { mnmn <- mnz2; vecz <- vec2}
mytab <- cbind(runconfig,mnmn,vecz,dat1,dat2,dat3,dat4)

good.ones <- mytab[,6]>.9 & mytab[,7]<1.2
#good.ones <- mytab[,5]>.9 & mytab[,9]>.9 & mytab[,6]<1.2
minitab <- mytab[good.ones,]

#cat("results based on best across DUPs and DELs:\n")
#print(minitab[order(rowSums(minitab[,c(8,9,12,13)])),])

minitab <- minitab[,-1]
colnames(minitab) <- c("Sample-QC","SNP-QC","PC-comps","CNV-QC","QS","OR","DGV-hits",
                       "Sens","DGV-valid","Spec")

minitab <- as.data.frame(minitab)
minitab[,1] <- c("none","relaxed","stringent")[minitab[,1]+2]
minitab[,2] <-c("none","basic","full")[minitab[,2]+2]
minitab[,3] <-c(0,6,24)[minitab[,3]+2]
minitab[,4] <-c("none",NA,"full")[minitab[,4]+2]
xtable(minitab[order((minitab[,5])),],digits=c(0,0,0,0,0,3,2,0,3,0,3))

#runconfig[good.ones]
if(DELs) {
  cat("results based on best across DELs:\n")    
  print(minitab[order((minitab[,ord])),])
} else {
  cat("results based on best across DUPs:\n")
  print(minitab[order((minitab[,ord])),])
}

if(DELs) { 
  #res <- get.correspondence(all.runs.del)[[1]]
} else {
  #res <- get.correspondence(all.runs.dup)[[1]]
  res <- res2
}
print(dim(res))
print(res[which(good.ones),which(good.ones)])
     
test.set <- c(10,22,4,8,20,2)
print(res[test.set,test.set])

print((sum(res[test.set,test.set])-6)/30)

rss <- ((rowSums(res[test.set,test.set]))-1)/5
css <- ((colSums(res[test.set,test.set]))-1)/5
print(rss)
print(css)
print(rss+css)



if(F) {
  ## REGRESSION RESULTS FOR FIRST YEAR REPORT ##
  
  oth <- cbind(runconfig,mnz,mnz2,vec,vec2)
  colnames(oth)[6:9] <- c("DEL.qs","DUP.qs","DEL.OR","DUP.OR")
  o1 <- lm(DEL.qs ~ Sample + SNP + PCA + CNV,data=oth)
  o2 <- lm(DUP.qs ~ Sample + SNP + PCA + CNV,data=oth)
  o3 <- lm(DEL.OR ~ 1+ Sample + SNP + PCA + CNV,data=oth)
  o4 <- lm(DUP.OR ~ 1+ Sample + SNP + PCA + CNV,data=oth)
  xtable(o1)
  xtable(o2)
  xtable(o3)
  xtable(o4)
  
  
  fm1 <- lmer(pc.del ~ db + src + Sample + SNP + PCA + CNV + (1|run) + (0+db|run), 
              data=runconfig2)
  
  fm1 <- lmer(pc.del ~ db + src + Sample + SNP + PCA + CNV + (1|run) + (0+db|run), 
              data=res.list2)
  fm2 <- lmer(pc.dup ~ db + src + Sample + SNP + PCA + CNV + (1|run) + (0+db|run), 
              data=res.list2)
  fm3 <- lmer(pc.del2 ~ db + src + Sample + SNP + PCA + CNV + (1|run) + (0+db|run), 
              data=res.list2)
  fm4 <- lmer(pc.dup2 ~ db + src + Sample + SNP + PCA + CNV + (1|run) + (0+db|run), 
              data=res.list2)
  
  for (cc in 1:4) {
    f1 <- attr(summary(get(paste("fm",cc,sep=""))),"coefs")
    f1[,2] <- round(f1[,2],4); f1[,3] <- round(f1[,3],2)
    f1[,1] <- round(f1[,1],3)
    print(xtable(f1,digits=c(1,3,4,2)))
  }
  
  
  lm(,data=res.list)
  
  res.list2 <- cbind(res.list,runconfig[res.list$run,])
}
