#!/bin/bash
  
cd /data/

module load samtools

for file in $(ls *sortedByCoord.out.bam)
do
samtools index $file
done
exit
