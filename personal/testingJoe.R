geno.files <- c("/chiswick/data/ncooper/JoeTest/LRRDATA/RAWDATA/plumbcnv_GENICA_input_test.txt",
"/chiswick/data/ncooper/JoeTest/LRRDATA/RAWDATA/plumbcnv_RBCS_input_test.txt")


dir <- make.dir("/chiswick/data/ncooper/JoeTest/",dir.raw="/chiswick/data/ncooper/JoeTest/LRRDATA/RAWDATA/")


snpMatLst <- import.snp.matrix.list(snpIDs,dir=dir,data=NULL,samples=NULL,n.cores=1,snp.fields=field.list)


snpMatLst <- import.snp.matrix.list(snpIDs,dir=dir,data=NULL,samples=unlist(sample.fn),n.cores=1,snp.fields=field.list)

SM <- read.snps.long(cat.path(dir$col,geno.files[1]),sample.id="92372",snp.id=readLines(cat.path(dir$ano,"snpNames.txt")),fields=fl,skip=2,in.order=T,every=42,codes="nucleotide",sep="\t",comment="#",simplify=c(F,F))