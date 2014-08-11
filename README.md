========
plumbCNV
========

## plumbCNV pipeline for comprehensive quality control and CNV detection using SNP arrays (R package) ##

https://github.com/nicholasjcooper/plumbCNV

This package has not yet been completely encased in the usual R format and directory structure, and so has not been submitted to CRAN.

To use this now, you must install by copying the files to a directory on your machine. You must also copy the file iFunctions.R from the repository: nicholaslcooper/iChip, and it's probably easiest to copy this into the same directory as the plumbCNV scripts.


Running instructions
====================

At the moment you need to source() the individual function files: FunctionsCNVAnalysis.R, geneticFunctions.R, validation.functions.R, iFunctions.R and tdtFunctions.R so that all required functions are loaded into the session (if you modify the path in the header of FunctionsCNVAnalysis.R, the other scripts should be sourced automatically when this file is sourced). Package installation may or may not happen automatically - try 'load.all.libs()' and make sure everything works, and ensure that packages are up to date, particularly for the bioconductor packages, especially for the human genome BSgenome.Hsapiens.UCSC.hg18/hg19 reference which has changed since july 2014.

To prepare for running, work your way through the instructions below. Feel free to contact me if you can't find your answer in this document [nick.cooper@cimr.cam.ac.uk]. Once it's all in place it should be quite easy to use, you just need to use one function 'plumbCNV()' which does everything really; it takes the input files as parameters, as well as various thresholds.

'/personal/multiRunsLightSnp.R' and '/personal/multirunsFam.R' are examples of how I run the pipeline once everything is installed and the input files are in place. It might be easiest to just modify these files to match your own filenames, etc.

These run on my system with 22 cores ('n.cores' parameter) for most tasks, and automatically submit to a 100-core 'grid' for others ('q.cores' parameter), but you don't have to use a grid if you don't have one, in which case you would set q.cores=0.


Public software required to be installed on your system:
--------------------------------------------------------

PennCNV - http://www.openbioinformatics.org/penncnv/penncnv_download.html

Plink - http://pngu.mgh.harvard.edu/~purcell/plink/

BLAS - you may have BLAS already as part of R or linux, if not: http://www.openblas.net/


R packages required:
--------------------
NB: need at least R v3.0, but this code should install everything you need

#reader and NCmisc
install.packages(reader,dependencies=TRUE)

install.packages(bigpca,dependencies=TRUE)

require(reader)

#Other packages
further.packages <- c("bigmemory","biganalytics","multicore","lattice","compiler","NCmisc")

bioC.packages <- c("IRanges","BiocGenerics","Biobase","GenomicRanges","genoset","bigalgebra")

#these NCmisc functions should install all these packages with this command
must.use.package(further.packages) 

must.use.package(bioC.packages,TRUE)



INPUT FILES REQUIRED
====================

MAIN DATA FILES
---------------

The main data for CNV analysis is the Log-R-Ratio (LRR) and Beta Allele Frequency (BAF) intensity data. These are calculated from X,Y or R,Theta values (the same raw intensities used for genotype calling), combined with an allele frequency reference. Genome studio files for a reference dataset will often have the LRR and BAF data already included, but if not, it may require a bit of work to calculate these yourself. This is not something supported by plumbCNV. Note that plink .bim/.bam files do not contain the intensity data for CNV calling, so if your dataset is in this format you may need to contact a data manager or administrator to obtain the raw intensity data that was used to genotype your samples.

The LRR and BAF intensity data can be read into plumbCNV from several formats.

* Genome studio file - this is recommended as this file contains snp ids, sample ids and allele codes, providing everything the program needs in a single source file. Create a file.spec.txt that indicates the column number where each FILE = raw.data.file.name; GRP = cohort.of.samples.genotyped.at.same.centre; TYPE = zip/txt file; SAMP=sample-id column number; SNP = snp-id column number, A1: allele-1 column number; A2: allele-2 column number; LRR: LRR column number; BAF: BAF column number.
* Long format text file - Similar to a genome studio file, but perhaps with different columns/ordering, it may even be a single column as long as you have files showing the consistent order of SAMPLE-ID and SNP-ID, when this text file is arranged by SAMPLE-ID and SNP-ID, with SNP-ID cycling through every ID before the next SAMPLE-ID is listed. Create a file spec just like for the genome studio file, however if this file does not contain allele data you will need to skip the SNP-QC step, or provide your own snq-qc exclusions files, by placing them (list of sample/snp ids) into the appropriate subdirectories of /ANNOTATION/
* Matrix format text file - Separate matrix files for LRR and BAF; if you are using this format you will not be able to use file.spec.txt to input your data information, you will need to read more carefully the input parameters for the plumbCNV() function, to make sure you are manually inputting all the required file names. See above for the long format text with regards to your options regarding SNP-QC.
* big.matrix objects containing LRR and BAF intensities - These can be inputted manually allowing skipping of the initial data import step (would start the pipeline from step 1, instead of step 0). 

You can split the raw input files into as many subsets as you like and still analyse them as if they are one file, for instance you may have cohorts:  Cases N=12000, Controls N=4000; in this case you can split the cases into 3 separate files, ensuring that the last sample is not split between the end of one and start of the next file, then in file.spec.txt you would have 3 row entries for the case dataset and 1 for controls, and you would indicate that the cases are all from the same group using the 'GRP' column and allocating the same code. If cases or controls come from multiple sources, code the GRP as different numbers as this will handle batch effects separately for each source, then at the end stage each will be combined and the pheno.lookup.txt file in /ANNOTATION will be the way that cases and controls are denoted.

PERFORMANCE AND TIMINGS
-----------------------

This pipeline is build to handle very large microarray datasets. Initial testing has been done on 16,000 samples for a 200K SNP array. There is no reason foreseeable why 200,000 samples for a million SNPs would not work using the same code, as everything is scalable and NOT linearly dependent on RAM (it is linearly dependent in hard-disk space). 
4GB of RAM is a minimum requirement (use low.ram=TRUE) and things will run more quickly the more memory you have. 20GB or more is recommended for the best performance. Terabytes of Hard Disk space may be required for large datasets. The code is strongly parallel, and running with multicores is ideal. Most of the multicore operation assumes that the cores are immediately accessible (for instance via the R-package 'multicore'). Use of GRIDs/Clusters is supported for limited operations within the pipeline, mainly to run PennCNV. To utilize this functionality you may need to follow the instructions for CLUSTERS/GRIDs just below. Most of my running of the pipeline at the DIL lab has used specs like:

Linux 64bit server
* Largest dataset: 16,000 samples x 200K SNPs
* 80GB RAM [typically ~20gb will be used by the program]
* 32 local cores (multicore) [I set n.cores=22]
* 150 cluster cores [I set q.cores=100]

Using this setup, the data import usually takes a few hours, and then steps 3-6 a few hours, in total, typically 5-8 hours for a full run through. The first time you run is slower because various types of annotation need to downloaded and calculated. If settings remain the same, you can also run with option 'restore=TRUE' will always attempt to use existing datafiles and calculations rather than regenerate new ones, often saving much time. Efficiencies are mostly designed for large datasets, so this pipeline may seem unusually slow if only run for several hundred samples or less.


CUSTOMISING THE CLUSTER/GRID SUBMISSION FUNCTION
------------------------------------------------
Regarding use with a cluster or GRID, one of the arguments for plumbCNV() is 'cluster.fn'. The default value is 'q.cmd', which is the function that runs on my queue. But, you can change this to anything you like. So if you write and
'source()' a function something like the one below, then you can just set cluster.fn="my.slurm".

my.slurm <- function(file.name,output.dir="",id="") {

    return(paste("sbatch -C nickC --wrap -o ",output.dir,"
",file.name,sep=""))

}

You just need to make sure that the function has the same arguments as the original function 'q.cmd' and that the output makes a call to the cluster that results in the output file being written to 'output.dir'.

Alternatively you can set q.cores=0 to avoid using the cluster, and PennCNV will still be run in parallel using however many cores you have available. Another alternative is setting the 2 arguments: run.manual=TRUE,print.cmds=TRUE; which will print the bash/putty penn cnv commands to the console and allow you to run the CNV calling manually using whatever method you like. Although this option is largely untested.


SYSTEM LIMITATIONS AND REQUIREMENTS
-----------------------------------

* this pipeline will not work in MS Windows (except perhaps via putty to a linux server), it requires multiple linux commands installed. It WILL work on MAC OS X and Linux.

* make sure PennCNV is installed for 64bit if your machine is 64-bit (default download is currently 32-bit)

* if any single raw data files contain more than roughly 1,000,000,000 samples x snps, you may need to use the option to run SNP-QC in plink, as the SnpMatrix object is limited by the maximum permitted size of R-objects. Alternatively you can split the raw input files into as many subsets as you like and still analyse them as if they are one file, see instructions above




SUPPORT AND ANNOTATION FILES
----------------------------

# NB: Please view this section as plain text (download the README.md text file), using a fixed-width font, as otherwise it will not give a correct impression of the correct input file formats!

Note also whether each has a header line (or not)

SNP support data (namely: snp-id,chr,pos) e.g, bim,vcf,map,map3, etc. does not matter what it is called but you eventually need to enter the file location into the plumbcnv() function as a parameter.

==> snpdata.map <==
1	rs642690	205835793
1	rs6427082	165562462
1	rs6427160	167158805
1	rs6427196	167747847
.        ...              ...


"plate.lookup.txt" - a list of what microarray plate each ID in your study was acquired on ['well' is optional]

==> plate.lookup.txt <==
id	        plate	well
WTCCCT542884	125454 	A12
WTCCCT542896	125454 	B12
WTCCCT542908	125454 	C12
WTCCCT542920	125454 	D12
...             ..     ..


"pheno.lookup.txt" - a list of which phenotype code applies to each ID in your study:

==> pheno.lookup.txt <==
phenotype
WTCCCT542884 0
WTCCCT542896 0
WTCCCT542908 1
WTCCCT542920 1
...


"sex.lookup.txt" - a list of which sex applies to each ID in your study (optional but recommended)

==> sex.lookup.txt <==
sampleid	sex
WTCCCT542884	male
WTCCCT542896	female
WTCCCT542908	male
WTCCCT542920	male



FILES TO SUPPORT IMPORT OF BAF AND LRR DATA
-------------------------------------------

If you will be importing your data from raw genome studio files, then you'll need "file.spec.txt" to show, file formats, which columns contain the variables of interest, and which files belong to which group (doesn't have to be the same as phenotype groupings, this is more for separate data sources, e.g, a different centre, or data acquired in different years, etc.

==> file.spec.txt  <==
FILE	GRP	TYPE	SAMP	SNP	A1	A2	LRR	BAF
sanger-controls.txt	1	txt	1	2	6	7	13	12
cases-part1.txt.gz	2	gzip	1	2	6	7	13	12
cases-part2.txt.gz	2	gzip	1	2	6	7	13	12
uva-controls.txt	3	txt	1	2	5	6	12	11



OR if importing from text files of just BAF/LRR, you'll need : "snpNames.txt" and "*.ids" instead, with contents corresponding to the order of the data in the files, and file names for the ID files being the same as the datafiles, but with '.ids' appended. If the IDs are already inside the raw files as row/column names, you don't need these files.

==> snpNames.txt <==
chr1:109457160
chr1:109457233
chr1:109457614
...

==> sanger-controls.txt.ids <==
WTCCCT542884
WTCCCT542896

==> cases-part1.txt.gz.ids <==
WTCCCT542908
WTCCCT542920

==> cases-part2.txt.gz.ids <==
WTCCCT542901
WTCCCT542927

==> uva-controls.txt.ids <==
59360811
50385211
52660511
55018311





Some important input parameters for the plumbCNV() function
-----------------------------------------------------------

dir.base="/data/ncooper/ImmunochipRunTest/" # base plumbCNV directory where all the subdirectories of results will be written

dir.raw="/data/Immunochip/FinalReports/" # location of the raw data, genome studio, long file or matrix format

aux.dir="/home/ncooper/Documents/otherFilesICHIP" # location of optional/additional support files not already in ANNOTATION, this can lie outside the main folder structure, which can be particularly convenient when running the pipeline on a dataset for the first time

snp.support="/data/store/metabochip/PLINK/MetaboChip.map" # snp support file name (eg, map, bim, etc), when using gsf=FALSE

gsf=TRUE  # TRUE/FALSE, whether you have prepared a specific snp annotation file (use FALSE), or whether the 
program should try to generate this file using a genome studio file contain SNP id, Chr, and Pos information (use TRUE)

run.mode="scratch" # can be "scratch","normal" or "big", depending on whether you are reading from: i) raw files; 2) processed, long format files or LRR/BAF matrices; or (iii) LRR/BAF big.matrix objects

snp.run.mode="normal" # can be "normal", "skip", or "plink"; depending on whether you are happy to run the QC in snpStats (maximum 1 billion genotypes), or whether you want to use plink SNP-QC, or none. If using PLINK, also make sure that 'plink.imp=TRUE'.

ped.file="ped.fam" # a plink ped file (family file) giving the family structure if trios are being used, must also use 'trio=TRUE' to use proper trio-calling, or also add 'joint=TRUE' to use 'joint' calling right from the beginning (this option is extremely slow and may take weeks for very large datsets).

start.at=0  # where to start the pipeline, there are 7 steps, 0:6, and you can start from any step you have reached previously

pause.after=6 # where to pause the pipeline (default is step 6, to run through to the end, but you can also opt to pause and review progress, appropriateness of thresholds and results at an earlier step)

penn.path="/usr/local/exports/bin/penncnv/"  ## the location of your PennCNV installation

hwe.thr=10^-8   # the hardy weinberg threshold for SNP-QC, note that for common CNVs this should be a very low value

build="hg18" # build 36/hg-18 is the default, so if you are using build 37/hg19, make sure you modify this

rare.pc=0.01  # the upper threshold for CNVs to detect. This is 1%. 5% may also be a commonly desired threshold, or use something like .95 if looking for common CNVs/CNPs.

num.pcs=24  # the number of principle components to correct for. The more you use the cleaner the data will be, our testing showed that 24 was best for our dataset. Yours may be different. This will have a huge effect on quality scores. If you are seeking common CNVs/CNPs you may need to set this threshold lower to avoid cleaning away the variants you are seeking.

pc.to.keep=.15  # the percentage of data (randomy selected) to use for the PCA analysis - this should always be substantially less than 100% to avoid removing all variation. I would suggest lower percentages particular if looking for common CNVs/CNPs.

exclude.bad.reg=TRUE # whether to exclude known regions with high false positives (immunoglobin, telomeres, centromeres) these are known to produce spurious results that are hard to QC, so default action is to exclude them. You may want them in which case, set this to FALSE.

min.sites=6 # the minimum number of SNPs required to call a CNV. This might be as low as 3 (limit of PennCNV), but in the literature 10 and 20 have been more common thresholds, as increasing the number of sites increases the confidence in the calls. My testing showed that using 6 sites gave a good balance of sensitivity and specificity and false positives can be cleaned up by using a quality score filter at the final stages of analysis.

rare.qc=TRUE # one of the steps recommended is to discard samples and places with too many rare CNVs, as this is a good indicator of poor quality data. However, there is some chance of discarding real data, so you may wish to turn this step off, if so set to FALSE.

restore=FALSE # this option will reuse any available pre-calculated data which can speed up subsequent runs of the same pipeline step(s). Note that you should set this to FALSE if you have made a substantial change to the settings since the last run, as the old settings may be inadvertantly used instead.

n.cores=22 # the number of parallel cores to use (these must be accessible by the r-package 'multicore')

q.cores=100 # number of cluster cores to use (allows running of the penn-cnv HMM command in massive parallel), this also requires that the function name passed to 'cluster.fn' passes a call to your cluster with the right syntax, e.g, see the default cluster function 'q.cmd' for example.

hmm="/usr/local/bin/penncnv/lib/hh550.hmm" # PennCNV comes with some built-in hidden markov model files, like hh550.hmm and hhall.hmm; there doesn't seem to be much difference in performance, the main parameter differences are for the LOH state. I also recommend testing the sensitivity of using different parameters, particularly as the PC-correction reduces the variance somewhat, which i think allows increase of sensitivity without increasing false positives - I have had good results moving the parameter means for DELs and DUPs closer to zero than the default HMM file. You can test this iteratively by setting start.at=5, pause.after=6, restore=TRUE, and varying this parameter file.


An example call to the main function for a rare CNV pipeline
------------------------------------------------------------
cnv.result <-
plumbCNV(dir.raw="/data/Immunochip/FinalReports/", dir.base="/data/ncooper/ImmunochipRunTest/", snp.support="/data/store/metabochip/PLINK/MetaboChip.map",aux.dir="/home/ncooper/Documents/otherFilesICHIP", gsf=FALSE,start.at=0,pause.after=6,
penn.path="/usr/local/bin/penncnv64/",hwe.thr=10^-8,hmm="/usr/local/bin/penncnv/lib/hh550.hmm",
build="hg19",rare.pc=0.05,num.pcs=24,pc.to.keep=.20,exclude.bad.reg=T,min.sites=6,rare.qc=TRUE,restore=FALSE)


An example call to the main function for a common CNP pipeline
--------------------------------------------------------------

cnv.result <-
plumbCNV(dir.raw="/store/ccge_vol1/icogs/cnvdata/",allele.codes=c("A","B"),
dir.base=base.dir, snp.support="/store/ccge_vol2/research/icogs/cnv/plumb_CNV/icogs_map_file_b37.map",
gsf=FALSE, aux.files.dir="/store/ccge_vol2/research/icogs/cnv/plumb_CNV/aux1/", start.at=0, pause.after=6,
penn.path="/usr/local/exports/bin/penncnv/", hwe.thr=10^-100, hmm="/usr/local/bin/penncnv/lib/hh550.hmm",
build="hg19",rare.pc=0.95, num.pcs=4, pc.to.keep=.10, exclude.bad.reg=T, min.sites=8, rare.qc=FALSE, restore=FALSE)


Note that another option is to use the alternative command 'plumbcnv()' lower case, which is the same as plumbCNV(), except that it allows the use of an input object called 'settings' which is simply a list containing your settings, this can be tidier. So you can nicely organise your settings in lists:

main.settings <- list(dir.raw="/store/ccge_vol1/icogs/cnvdata/", allele.codes=c("A","B"),dir.base= ... , etc)

penn.settings <- list(penn.path="/usr/local/exports/bin/penncnv/", hmm="/usr/local/bin/penncnv/lib/hh550.hmm",q.cores=0)

pca.settings <- list(num.pcs=12, pc.to.keep=.15)

cnv.qc.settings <- list(rare.pc=0.05, exclude.bad.reg=FALSE, min.sites=10, rare.qc=TRUE)

ALL.SETTINGS <- c(main.settings,penn.settings,pca.settings,cnv.qc.settings)


#And then use the simplified input:

cnv.result <- plumbcnv(settings=ALL.SETTINGS)


Trouble-shooting: Problems encountered by others
------------------------------------------------

* make sure you enter the correct value for build 36/37  (hg18/hg19) into the main function

* make sure file.spec.txt is in tab separated format with no extra spaces

* if your alleles are not coded in the genome studio file as A,C,G,T, make sure you enter a value for allele.codes=c() in the main function

* plumbCNV is designed for reasonably large datasets, if you are using less than 100 samples or less than 10,000 SNPs you may encounter unexpected behaviour, as there has been little testing of these sorts of datasets.


Description and Examples of output
----------------------------------

As far the output that you should expect, there is quite a lot. So while running, it is quite verbose letting you know what is going on, there are 7 main steps:

* 0) raw LRR, BAF data conversion to long file or matrix format

* 1) import of matrix or long file LRR, BAF, data

* 2) import of genotype data and analysis of SNP-QC

* 3) sample QC analyses

* 4) principle components correction

* 5) running of PennCNV

* 6) CNV-qc, summaries, overlaps, results


When you run the plumbCNV function you can select the start and end step, so you can do the analysis a step at a time to monitor progress, sometimes repeating a step adjusting thresholds, etc.

For instance:

plumbCNV(start.at=1, pause.after=2,...)

will just run steps 1 and 2 above, or:

plumbCNV(start.at=0, pause.after=6,...)

... will run from start to finish.


At the same time as doing the processing, various tables and figures are generated automatically, if you look through the output you will see lots of lines like '~wrote file: ImmunochipFamilies/SAMPLEQC/MEAN_DLRS_GC/LRRMean+DLRS+GCWave.pdf ', showing that a table or figure was produced. These output files are always stored in the same standardized folder structure (the directories are automatically created by the plumbCNV() function in your base directory). So there are tables of means, SDs, thresholds, QC figures, tables of call rates, graphs of outlying LRR plots, etc.

The final step, number 6 produces the bulk of results. This includes R object files containing the CNVs called (listed separately for deletions and duplications), CNV regions, summaries of overlaps with genes, exons and the database of genomic variants (if those options are selected), and it will also perform regional and global association tests with your phenotype, and by length, although obviously all the data is saved so you can do your own analyses. There are quality scores in the files showing how reliable each CNV will be. There is absolutely no worthwhile data produced by the pipeline that is not saved and accessible to you, the analyst.

The best way to run would be to look at my files like /personal/multiRunsLightSNP.R from github, showing the script i use to setup a run of plumbCNV, with all the settings I use accessible by editing the file. This script contains some if-statements that control multiple settings at once (which was how I automated the 54 runs in my upcoming (hopefully) HMG publication).

All of the raw data is saved in bigmemory format so can be accessed quickly, and there are plotting functions allowing you to plot raw LRR and BAF data for CNVs, showing overlapping genes, or allows plotting of CNVs against genome location colouring for cases and controls, e.g, plot.pheno.cnvs(), cnv.plot(), 



##############################################################################

