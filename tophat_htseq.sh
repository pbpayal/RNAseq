cd /data3/Payal/

#----------------------------------------TOPHAT---------------------------------------------
for file in $(ls *001.fastq.gz)
do
  echo "My file name is $file"
  echo "Running ----HISAT2-----"
  hisat2 --threads 8 --time -x /data/Genomes/Mouse/GRCm38/Hisat2_Index/Mus_musculus.GRCm38.dna.primary_assembly -U $file -S aligned_$file.sam --summary-file summary_$file.txt
  echo "The output is aligned_$file.sam"
done
echo "------------------------------Finished Running TOPHAT------------------------------"

#----------------------------------------Samtools Sort and Convert BAM---------------------------------------------
for file in $(ls aligned_*)
do
  echo "Running ----SAMTOOLS SORT AND CONVERSION TO BAM----"
  echo "My file name is $file"
  samtools sort -o sorted_$file -O bam $file
  echo "The output is sorted_$file"
  
done
echo "------------------------------Finished sorting and converting files with SAMTOOLS------------------------------"

#----------------------------------------Htseq Counts---------------------------------------------
for file in $(ls sorted_*) 
do
  echo "Running ----HTSEQ COUNTS----"
  echo "My file name is $file"
  htseq-count --format bam --order pos -t exon $file /data/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.91.gtf > counts_$file.txt
  echo "The output is counts_$file.txt"
done
echo "------------------------------Finished HTSEQ Counts------------------------------"
echo "------------------------------Pipeline Over------------------------------"
