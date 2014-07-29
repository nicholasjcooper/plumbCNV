#SCRIPT to load in all directories used across scripts#

mode <- 2

if(mode==1) {
  dir <- def.dirs.fn(dir.raw="/ipswich/data/Immunochip/FinalReports/",
                   dir.sup="/ipswich/data/Immunochip/FinalReports/support/",
                   dir.base="/chiswick/data/ncooper/ImmunochipReplication/")
}

if(mode==2) {
  dir <- def.dirs.fn(dir.raw="/ipswich/data/Immunochip/FinalReports/",
                   dir.sup="/ipswich/data/Immunochip/FinalReports/support/",
                   dir.base="/dunwich/scratch/ncooper/ImmunochipReplication/")
  nuthin <- init.dirs.fn(dir,overwrite=F,silent=T,
              info.dir="/home/ncooper/Documents/necessaryfilesICHIP")
}

if(mode==3) {
  dir <- def.dirs.fn(dir.raw="/chiswick/data/store/metabochip/FinalReports/",
                   dir.sup="/chiswick/data/store/metabochip/PLINK/",
                   dir.base="/chiswick/data/ncooper/metabochipCNVanalysis2012/")

}

# create directories required (if not already present)
nuthin <- init.dirs.fn(dir,overwrite=F,silent=T)  #,ignore=c("sup","raw"))


rm(mode)