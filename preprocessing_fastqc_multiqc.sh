#!/bin/bash
  
module load fastqc/0.11.9
cd /data/

for file in $(ls *.fastq.gz)
do
#  echo "My file name is $file"
#  echo "Running fastqc....."
  fastqc $file
done

module load multiqc

multiqc *
