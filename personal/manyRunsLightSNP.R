#~/R-2.15.3/bin/R --no-save --slave < manyRunsLightSNP.R > domany.1s.log

ones12 <- F
tweny20s <- F
firdy7s <- T
six6 <- F

if(ones12) {
  suffix <- "RUN01"
  
  comp <- F
  samp.set <- "heavy"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN02"
  
  comp <- F
  samp.set <- "heavy"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN03"
  
  comp <- F
  samp.set <- "heavy"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN04"
  comp <- F
  samp.set <- "heavy"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN05"
  
  comp <- F
  samp.set <- "heavy"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
}

if(six6) {  
  suffix <- "RUN06"
  
  comp <- F
  samp.set <- "heavy"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")

}



if(tweny20s) {
  suffix <- "RUN19"
  
  comp <- F
  samp.set <- "light"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- F
  do.cnv <- F
  restore <- F
 
 
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")

  suffix <- "RUN20"

  comp <- F
  samp.set <- "light"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- F
  do.cnv <- T
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")

  suffix <- "RUN21"
 
  comp <- F
  samp.set <- "light"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")

  suffix <- "RUN22"
  comp <- F
  samp.set <- "light"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")

  suffix <- "RUN23"
 
  comp <- F
  samp.set <- "light"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
 
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN24"
  
  comp <- F
  samp.set <- "light"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
}

if(firdy7s) {

  suffix <- "RUN37"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
 # source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN38"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN39"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN40"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN41"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
  suffix <- "RUN42"
  
  comp <- F
  samp.set <- "none"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRunsLightSNP.R")
  
}
