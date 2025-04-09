#!/bin/bash

#Fisrt make the intervals list, the header from this file must be extracted from the bam files
#and the body of the file is extracted from the gtf file .
#So the first thing you need to do is get the header from your SAM/BAM file:
#samtools view -H [your.bam] > intervalList.txt
#If your GTF file is standard and we assume that it contains only ribosomal intervals, 
#then we need the first, fourth, fifth, seventh, and ninth fields from the file. We can append them onto our text file which contains the header:
#cut -s -f 1,4,5,7,9 [your.gtf] >> intervalListBody.txt

cd /data/Sorted_BAMS/

for file in $(ls sorted_aligned*)
do
  echo "My file name is $file"
  echo "Running CollectRnaSeqMetrics....."
  picard CollectRnaSeqMetrics \
  I=$file \
  O=RNA_Metrics_$file.txt \
  REF_FLAT=/data/Genomes/Mouse/refFlat_chr_rem.txt \
  STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=data/Genomes/Mouse/intervalList.txt \
  CHART_OUTPUT=RNA_Metrics_$file.pdf
done
echo "------------------------------Picard RNA Metrcis Done------------------------------"
