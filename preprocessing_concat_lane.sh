#!/bin/bash

cd /data/

#Concat R1
for f in $(ls *L003_R1_001.fastq.gz)
do
r1=$f;
r2=${f/L003_R1_001.fastq.gz/}L004_R1_001.fastq.gz;
#22F2TFLT4_19996445_S115_L003_R1_001.fastq.gz 
n1="${r1}"; n1=$(echo $r1| cut -d '_' -f 3);
echo "The sample name part1 is $n1"
n2="${r1}"; n2=$(echo $r1| cut -d '_' -f 2);
echo "The sample name part2 is $n2"
n3=$n1'_'$n2
echo $n3
echo "The R1 file is $r1, R2 file is  $r2";
cat $r1 $r2 > lane_merged_$n3.R1.fastq.gz
echo "The final filename is lane_merged_$n3.R1.fastq"
done

#Concat R2
for f in $(ls *L003_R2_001.fastq.gz)
do
r1=$f;
r2=${f/L003_R2_001.fastq.gz/}L004_R2_001.fastq.gz;
#echo "The R1 file is $r1, R2 file is  $r2";
n1="${r1}"; n1=$(echo $r1| cut -d '_' -f 3);
echo "The sample name part1 is $n1"
n2="${r1}"; n2=$(echo $r1| cut -d '_' -f 2);
echo "The sample name part2 is $n2"
n3=$n1'_'$n2
echo $n3
cat $r1 $r2 > lane_merged_$n3.R2.fastq.gz
echo "The final filename is lane_merged_$n3.R1.fastq"
done

