#!/bin/bash

echo "Track system usage for 2880 x 30s = 24 hours"
top -b -d 1 -n 1 -u ncooper > temp.top  
head -1 temp.top
grep 'ncooper' temp.top | head -1 > progress.top

for i in $(seq 1 2880);
do
  top -b -d 30 -n 1 -u ncooper > temp1.top  
  grep ncooper temp1.top | head -1 > temp.top
  cat progress.top temp.top > progress2.top
  rm progress.top
  mv progress2.top progress.top
done

rm temp.top