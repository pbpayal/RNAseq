#!/bin/bash

cd /home/pbanerjee/Payal/snakemake_test/griffithlab_brain_vs_uhr/HBR_UHR_ERCC_1_sample

GTF="/data/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.91.gtf"
GENOME_DIR="/data/Genomes/Mouse/STAR_Mouse_Index/"

for f in *.read1.fastq.gz
do
r1=$f;
r2=${f/.read1.fastq.gz/}.read2.fastq.gz;
echo "The R1 file is $r1 and the R2 file is  $r2";
STAR --runThreadN 8 --runMode alignReads --readFilesCommand zcat --limitGenomeGenerateRAM=119000000000 \
--genomeSAsparseD 3 --genomeSAindexNbases 12 --genomeChrBinNbits=16 --outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts --outFileNamePrefix star_out_$r1 \
--genomeDir $GENOME_DIR \
--sjdbGTFfile $GTF \
--readFilesIn $r1 $r2
echo "The output filename is star_out_$r1"
done

for file in $(ls *bam)
do
echo "Running ----HTSEQ COUNTS----"
echo "My file name is $file"
htseq-count --format bam --additional-attr=gene_name --order pos -t exon $file $GTF  > counts_$file.txt
echo "The output is counts_$file.txt"
done
echo "------------------------------Finished HTSEQ Counts------------------------------"
