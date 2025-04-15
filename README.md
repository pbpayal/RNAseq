**General Pipeline for RNAseq**

* Fastqc > STAR_2.7.8a > Samtools(1.21)Index > htseq_2.0.7 > Deseq2

* Fastqc > STAR_2.7.8a > Samtools(1.21)Index > RSEM > Deseq2

* Fastqc > STAR_2.7.8a > Samtools(1.21)Index > feature_counts > Deseq2

* Fastqc > HISAT2 > Samtools(1.21)Index > feature_counts > Deseq2

* Fastqc > HISAT2 > Samtools(1.21)Index > htseq > Deseq2

**Quality Trimming depending on library preparation protocol:**
1. Trimmomatic
2. Fastp
3. Cutadapt
4. BBDuK

**3' and 5' bias detection tools** - Picard CollectRnaSeqMetrics, RNAseqQC



