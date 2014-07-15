#!/bin/bash
  
#######################################################################
##               Bash Script to make a Plink fam file                ##
#######################################################################
  outdir=/chiswick/data/ncooper/metabochipRunTest/
  # get number of subject ids
  lll=$(wc -l "$outdir"/ANNOTATION/subIdsALL.txt | cut -d ' ' -f 1)
  # generate column of 1s
  printf '1\n%.0s' {$(seq 1 $lll)} > 1s.txt
  printf '0\n%.0s' {$(seq 1 $lll)} > 0s.txt
  ## have sex
  sexfile="$outdir"/ANNOTATION/sex.lookup.txt
  if [ -f $sexfile ]
  then
   awk -F, 'NR==FNR{a[$2]=$1;next}{print a[$1] $1;}' "$outdir"/ANNOTATION/subIdsALL.txt "$sexfile" | cut -f 2 | tail -"$lll" > sex.vec.txt
   ppp=$(wc -l sex.vec.txt | cut -d ' ' -f 1)
  fi
  if [ $ppp -eq $lll ] ; then echo sex file had correct length ; else cp 1s.txt sex.vec.txt ; fi
  ## have pheno
  phfile="$outdir"/ANNOTATION/pheno.lookup.txt
  if [ -f $phfile ]
  then
   awk -F, 'NR==FNR{a[$2]=$1;next}{print a[$1] $1;}' "$outdir"/ANNOTATION/subIdsALL.txt "$phfile" | cut -d ' ' -f 2 | tail -"$lll" > ph.vec.txt
   ppp=$(wc -l ph.vec.txt | cut -d ' ' -f 1)
  fi
  if [ $ppp -eq $lll ] ; then echo pheno file had correct length ; else cp 1s.txt ph.vec.txt ; fi
  ## generate family file
  paste subIdsALL.txt subIdsALL.txt 0s.txt 0s.txt sex.vec.txt ph.vec.txt  > snpdataCustom.fam
  #rm sex.vec.txt ph.vec.txt 0s.txt 1s.txt