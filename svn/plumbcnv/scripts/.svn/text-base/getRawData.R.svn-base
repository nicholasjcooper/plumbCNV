source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R")
source("/chiswick/data/ncooper/ImmunochipReplication/Scripts/DefineDirectoriesLocations.R")


#snp.info.fn <- paste(dir$sup,"sanger-controls.txt",sep="")
snp.info.fn <- paste(dir$sup,"Metabo_20100426_58C.bim",sep="") # support file name (eg, bim)
genome.stud.file <- F # T whether snp pos,chr,id file is a genome studio file (ie., not vcf/bim/map)

init.data.read <- function(dir,doLRR=T,doBAF=F,plink.imp=F,scr.name="getDataGS.sh",scr.dir="",
                           snp.info.fn,genome.stud.file=F,combine.files=F) 
{
  dir <- validate.dir.for(dir,c("scr","baf.col","lrr.dat","base","col","sup","raw"))
  
  if(combine.files) {
    combin <- "CD"
  } else {
    combin <- ""
  }
  if(is.null(scr.dir) | scr.dir==""){
    scr.dir <- dir$scr
  }
  if (genome.stud.file) { 
    gs <- "G" 
  } else { 
    gs <- "" 
    nc <- nchar(snp.info.fn)
    suf <- substr(snp.info.fn,nc-3,nc)
    if(toupper(suf) %in% c(".BIM",".VCF",".MAP")) {
      cat("Snp file is",suf,"type, which should be supported\n")
    } else {
      cat("Not expecting snp file with type",suf,", snp information may not be ")
      cat("imported successfully. Check results of this script carefully.\n")
    }
  }
  
  if(doLRR) {
    if(!plink.imp) {
      #import LRR for QC in R:     
      cat("Running bash data import script\n")
      my.cmd <- paste(" -lSM",gs,combin," -T 'LRR' -F '",dir$raw,"' -O '",
                      dir$base,"' -m '",snp.info.fn,"'",sep="")
      scr.file <- dir.fn.cmb(scr.dir,scr.name,must.exist=T)
      cmd <- paste(scr.file,my.cmd,collapse="",sep="")
      cat(cmd,"\n")
      system(command=cmd)
    } else {
      #import LRR and do plink QC:  
      cat("Running bash data import script and calling plink for snp QC\n")
      my.cmd <- paste(" -SLMRPlf",gs,combin," -T 'LRR' -F '",dir$raw,"' -O '",
                      dir$base,"' -m '",snp.info.fn,"'",sep="")
      scr.file <- dir.fn.cmb(scr.dir,scr.name,must.exist=T)
      cmd <- paste(scr.file,my.cmd,collapse="",sep="")
      system(command=cmd)
    }
  }
  
  lrr.filez <- c(list.files(dir$col),list.files(dir$lrr.dat))
  
  if(doBAF) {
    #import BAF: 
    cat("Running bash data import script for BAF data\n")
    if(combin!="") { combin <- paste("-",combin,sep="") }
    my.cmd <- paste(combin," -N 3 -T 'BAF'-F '",dir$raw,"' -O '",dir$base,"'",sep="")
    scr.file <- dir.fn.cmb(scr.dir,scr.name,must.exist=T)
    cmd <- paste(scr.file,my.cmd,collapse="",sep="")
    system(command=cmd)
  }
  
  baf.filez <- list.files(dir$baf.col)
  
  cat("complete!\n")
  out <- list(lrr.filez,baf.filez)
  names(out) <- c("lrr.list.files","baf.list.files")
  return(out)
}


REZ <- init.data.read(dir,doLRR=T,doBAF=F,plink.imp=F,scr.name="getDataGS.sh",
                      scr.dir="/chiswick/data/ncooper/ImmunochipReplication/Scripts/",
                      snp.info.fn,genome.stud.file=F,combine.files=F) 
  
print(REZ)

## write smart file,dir combination function!

system("R --slave < excludeSampsQC > F0.txt")