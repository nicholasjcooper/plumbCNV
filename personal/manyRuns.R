elev12 <- F
tweny25s <- F
fordy3s <- T

if(elev12) {
  
  suffix <- "RUN07"
  
  comp <- T
  samp.set <- "heavy"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN08"
  
  comp <- T
  samp.set <- "heavy"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN09"
  
  comp <- T
  samp.set <- "heavy"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN10"
  comp <- T
  samp.set <- "heavy"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN11"

  comp <- T
  samp.set <- "heavy"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")

  suffix <- "RUN12"

  comp <- T
  samp.set <- "heavy"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
}



if(tweny25s) {
  suffix <- "RUN25"
  
  comp <- T
  samp.set <- "light"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
 
 
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")

  suffix <- "RUN26"

  comp <- T
  samp.set <- "light"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")

  suffix <- "RUN27"
 
  comp <- T
  samp.set <- "light"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")

  suffix <- "RUN28"
  comp <- T
  samp.set <- "light"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F


  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")

  suffix <- "RUN29"
 
  comp <- T
  samp.set <- "light"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
 
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN30"
  
  comp <- T
  samp.set <- "light"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
}

if(fordy3s) {

#then
#RUN43
#RUN44
#RUN45
#RUN46
#RUN47
#RUN48

  suffix <- "RUN43"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 24
  my.st <- 3
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN44"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 24
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN45"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 6
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN46"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 6
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN47"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 0
  my.st <- 4
  my.end <- 6
  eval <- T
  do.cnv <- F
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
  suffix <- "RUN48"
  
  comp <- T
  samp.set <- "none"
  pca.set <- 0
  my.st <- 6
  my.end <- 6
  eval <- T
  do.cnv <- T
  restore <- F
  
  
  source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/multiRuns.R")
  
}
