
#### Function ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("pheatmap")
library("cowplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library("biomaRt")
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library("EnhancedVolcano")
library(readxl)
install.packages("writexl")
library(writexl)




######################################
# Automated DEG for each condition
######################################
jen_deseq2_auto_function <- function(count_matrix,sample_annotation,contrasts, outputPrefix, outputdirectory){
  rownames(count_matrix) <- count_matrix$gene_id
  count_matrix <- subset(count_matrix, select = -c(gene_id))
  #count_matrix <- round(count_matrix)
  rownames(sample_annotation) <- sample_annotation$sample
  sample_annotation <- subset(sample_annotation, select = -c(sample))
  all(rownames(sample_annotation) == colnames(count_matrix))
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = sample_annotation,
                                design = ~ condition)
  dds_deg <- DESeq(dds)
  print(resultsNames(dds_deg))
  res <- results(dds_deg, contrast=contrasts)
  #write.csv(res, file = paste0(outputdirectory,outputPrefix, "-summaryDEG.txt"))
  print(head(res))
  print(summary(res))
  #write.csv(summary(res), file = paste0(outputdirectory,outputPrefix, "-summaryDEG.txt"))
  res_sub = subset(res, padj<0.1)
  res_df <- as.data.frame(res)
  res_sub <- res_sub[order(res_sub$padj),]
  res_sub
  res_data <- as.data.frame(res_sub)
  res_data$gene_id <- rownames(res_data)
  head(res_data)
  ensembl=useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id","external_gene_name", "description"), values=res_data$gene_id, mart= ensembl)
  # rename colname  geneid to ensembl_gene_id to be able to merge on the same column name
  colnames(res_data)[7] <- c("ensembl_gene_id")
  # merge the annotation with your data
  resdata_anno <- merge(res_data, gene_list, by="ensembl_gene_id")
  print(head(resdata_anno))
  write.csv(resdata_anno, file = paste0(outputdirectory,outputPrefix, "-pval_0.1_results.csv"))
  #################################
  # MA Plot
  #################################
  png(paste0(outputdirectory,outputPrefix, "-MAplot.png"), height = 1200, width = 1200, res = 110)
  plotMA(dds_deg, ylim=c(-2,2),main = outputPrefix )
  dev.off()
  #################################
  # PCA Plot
  #################################
  #png(paste0(outputdirectory,outputPrefix, "-PCAplot.png"),height = 800, width = 600, res = 100 )
  dds_vst <- vst(dds_deg, blind = T)
  nudge <- position_nudge(y = 1)
  pcaplot_vst <- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), position = nudge, vjust=1) + geom_point(size=0.1)
  print(pcaplot_vst)
  ggsave(paste0(outputdirectory,outputPrefix, "-PCAplot.png"), pcaplot_vst, height = 10, width = 10, dpi = 320) 
  #dev.off()
  #################################
  #Heatmap of top 100 DEGs
  #################################
  ## Heatmap
  # Make a master GeneID to GeneName matrix
  resdata_anno_ordered = resdata_anno[order(resdata_anno$padj),]
  head(resdata_anno_ordered)
  head(resdata_anno_ordered$ensembl_gene_id)
  head(resdata_anno_ordered$mgi_symbol)
  gene_name_mapping_mat <- cbind(resdata_anno_ordered$ensembl_gene_id,resdata_anno_ordered$mgi_symbol)
  gene_map_matrix <- as.matrix(gene_name_mapping_mat)
  head(gene_map_matrix)
  nrow(gene_map_matrix)
  gene_map_df <- as.data.frame(gene_name_mapping_mat)
  colnames(gene_map_df) <- c("ensembl_gene_id","gene_name")
  mat = assay(dds_vst)[ head(order(res$padj),100), ]
  # nrow(mat)
  mat1 <- as.data.frame(mat)
  # Add a column named ensembl_gene_id to the matrix for mapping later
  mat1$ensembl_gene_id <- rownames(mat1)
  # Map the vsd ENSEMBL IDS to master Gene Name matrix
  p2 <- merge(gene_map_df, mat1, by="ensembl_gene_id", sort = FALSE)
  # Remove the ensembl_id column
  p2_subset <- subset(p2, select = -c(ensembl_gene_id))
  # Add the gene names as rownames
  rownames(p2_subset) <- p2_subset$gene_name
  # Remove the redundant gene name column now
  p2_subset2 <- subset(p2_subset, select = -c(gene_name))
  # convert into matrix as pheatmap requires a matrix as input
  p2_mat <- as.matrix(p2_subset2)
  # Optional to add additional condition subgroups in heatmap for easy visualization
  df = as.data.frame(colData(dds_vst)[,c("condition")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  # and plot the actual heatmap
  pheatmap <- pheatmap(p2_mat, annotation_col=df, main = "Top 100 Gene Expression", fontsize_col = 7, fontsize_row = 5.5, width = 8, height = 12, cluster_rows = FALSE)
  pheatmap
  ggsave(paste0(outputdirectory,outputPrefix, "-Heatmap.png"), pheatmap, height = 10, width = 10, dpi = 320)
  # dev.off()
  #################################
  # Volcano Plot
  #################################
  #png(paste0(outputdirectory,outputPrefix, "-Volcanoplot.png"), height = 900, width = 1200, res = 110)
  volcano_plot <- EnhancedVolcano(resdata_anno,
                                  lab = resdata_anno$mgi_symbol,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  pCutoff = 10e-5,
                                  FCcutoff = 1,
                                  pointSize = 2.0,
                                  labSize = 4.0,
                                  title = outputPrefix,
                                  subtitle = 'Differential expression')
  print(volcano_plot)
  ggsave(paste0(outputdirectory,outputPrefix, "-Volcanoplot.png"), volcano_plot, height = 10, width = 15, dpi = 320) 
  #dev.off()
}


# D_vs_E_Heart
outputdirectory <- "/Documents/2024/D_vs_E_Heart/"
outputPrefix <- "E_vs_D_Heart"
count_matrix <- read.csv(file="/Documents/2024/D_vs_E_Heart/D_vs_E_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/D_vs_E_Heart/D_vs_E_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Heart", "Grp_D_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# D_vs_E_Liver
outputdirectory <- "/Documents/2024/D_vs_E_Liver/"
outputPrefix <- "E_vs_D_Liver"
count_matrix <- read.csv(file="/Documents/2024/D_vs_E_Liver/D_vs_E_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/D_vs_E_Liver/D_vs_E_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Liver", "Grp_D_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_I_Heart
outputdirectory <- "/Documents/2024/C_vs_I_Heart/"
outputPrefix <- "C_vs_I_Heart"
count_matrix <- read.csv(file="/Documents/2024/C_vs_I_Heart/C_vs_I_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_I_Heart/C_vs_I_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_heart", "Grp_I_heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_I_Liver
outputdirectory <- "/Documents/2024/C_vs_I_Liver/"
outputPrefix <- "C_vs_I_Liver"
count_matrix <- read.csv(file="/Documents/2024/C_vs_I_Liver/C_vs_I_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_I_Liver/C_vs_I_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_liver", "Grp_I_liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# A_vs_B_Heart
outputdirectory <- "/Documents/2024/A_vs_B_Heart/"
outputPrefix <- "B_vs_A_Heart"
count_matrix <- read.csv(file="/Documents/2024/A_vs_B_Heart/A_vs_B_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/A_vs_B_Heart/A_vs_B_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_B_Heart", "Grp_A_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# A_vs_B_Liver
outputdirectory <- "/Documents/2024/A_vs_B_Liver/"
outputPrefix <- "B_vs_A_Liver"
count_matrix <- read.csv(file="/Documents/2024/A_vs_B_Liver/A_vs_B_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/A_vs_B_Liver/A_vs_B_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_B_Liver", "Grp_A_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# A_vs_C_Heart
outputdirectory <- "/Documents/2024/A_vs_C_Heart/"
outputPrefix <- "C_vs_A_Heart"
count_matrix <- read.csv(file="/Documents/2024/A_vs_C_Heart/A_vs_C_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/A_vs_C_Heart/A_vs_C_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_A_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# A_vs_C_Liver
outputdirectory <- "/Documents/2024/A_vs_C_Liver/"
outputPrefix <- "C_vs_A_Liver"
count_matrix <- read.csv(file="/Documents/2024/A_vs_C_Liver/A_vs_C_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/A_vs_C_Liver/A_vs_C_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_A_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# B_vs_C_Heart
outputdirectory <- "/Documents/2024/B_vs_C_Heart/"
outputPrefix <- "C_vs_B_Heart"
count_matrix <- read.csv(file="/Documents/2024/B_vs_C_Heart/B_vs_C_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/B_vs_C_Heart/B_vs_C_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_B_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# B_vs_C_Liver
outputdirectory <- "/Documents/2024/B_vs_C_Liver/"
outputPrefix <- "C_vs_B_Liver"
count_matrix <- read.csv(file="/Documents/2024/B_vs_C_Liver/B_vs_C_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/B_vs_C_Liver/B_vs_C_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_B_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# B_vs_D_Heart
outputdirectory <- "/Documents/2024/B_vs_D_Heart/"
outputPrefix <- "B_vs_D_Heart"
count_matrix <- read.csv(file="/Documents/2024/B_vs_D_Heart/B_vs_D_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/B_vs_D_Heart/B_vs_D_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_B_Heart", "Grp_D_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# B_vs_D_Liver
outputdirectory <- "/Documents/2024/B_vs_D_Liver/"
outputPrefix <- "B_vs_D_Liver"
count_matrix <- read.csv(file="/Documents/2024/B_vs_D_Liver/B_vs_D_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/B_vs_D_Liver/B_vs_D_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_B_Liver", "Grp_D_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_E_Heart
outputdirectory <- "/Documents/2024/C_vs_E_Heart/"
outputPrefix <- "C_vs_E_Heart"
count_matrix <- read.csv(file="/Documents/2024/C_vs_E_Heart/C_vs_E_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_E_Heart/C_vs_E_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_E_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_E_Liver
outputdirectory <- "/Documents/2024/C_vs_E_Liver/"
outputPrefix <- "C_vs_E_Liver"
count_matrix <- read.csv(file="/Documents/2024/C_vs_E_Liver/C_vs_E_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_E_Liver/C_vs_E_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_E_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_F_Heart
outputdirectory <- "/Documents/2024/C_vs_F_Heart/"
outputPrefix <- "C_vs_F_Heart"
count_matrix <- read.csv(file="/Documents/2024/C_vs_F_Heart/C_vs_F_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_F_Heart/C_vs_F_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_F_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_F_Liver
outputdirectory <- "/Documents/2024/C_vs_F_Liver/"
outputPrefix <- "C_vs_F_Liver"
count_matrix <- read.csv(file="/Documents/2024/C_vs_F_Liver/C_vs_F_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_F_Liver/C_vs_F_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_F_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_G_Heart
outputdirectory <- "/Documents/2024/C_vs_G_Heart/"
outputPrefix <- "C_vs_G_Heart"
count_matrix <- read.csv(file="/Documents/2024/C_vs_G_Heart/C_vs_G_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_G_Heart/C_vs_G_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_G_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_G_Liver
outputdirectory <- "/Documents/2024/C_vs_G_Liver/"
outputPrefix <- "C_vs_G_Liver"
count_matrix <- read.csv(file="/Documents/2024/C_vs_G_Liver/C_vs_G_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_G_Liver/C_vs_G_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_G_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_H_Heart
outputdirectory <- "/Documents/2024/C_vs_H_Heart/"
outputPrefix <- "C_vs_H_Heart"
count_matrix <- read.csv(file="/Documents/2024/C_vs_H_Heart/C_vs_H_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_H_Heart/C_vs_H_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Heart", "Grp_H_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# C_vs_H_Liver
outputdirectory <- "/Documents/2024/C_vs_H_Liver/"
outputPrefix <- "C_vs_H_Liver"
count_matrix <- read.csv(file="/Documents/2024/C_vs_H_Liver/C_vs_H_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/C_vs_H_Liver/C_vs_H_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_C_Liver", "Grp_H_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# E_vs_H_Heart
outputdirectory <- "/Documents/2024/E_vs_H_Heart/"
outputPrefix <- "E_vs_H_Heart"
count_matrix <- read.csv(file="/Documents/2024/E_vs_H_Heart/E_vs_H_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/E_vs_H_Heart/E_vs_H_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Heart", "Grp_H_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# E_vs_H_Liver
outputdirectory <- "/Documents/2024/E_vs_H_Liver/"
outputPrefix <- "E_vs_H_Liver"
count_matrix <- read.csv(file="/Documents/2024/E_vs_H_Liver/E_vs_H_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/E_vs_H_Liver/E_vs_H_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Liver", "Grp_H_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


# E_vs_I_Heart
outputdirectory <- "/Documents/2024/E_vs_I_Heart/"
outputPrefix <- "E_vs_I_Heart"
count_matrix <- read.csv(file="/Documents/2024/E_vs_I_Heart/E_vs_I_Heart.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/E_vs_I_Heart/E_vs_I_Heart_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Heart", "Grp_I_Heart")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)

# E_vs_I_Liver
outputdirectory <- "/Documents/2024/E_vs_I_Liver/"
outputPrefix <- "E_vs_I_Liver"
count_matrix <- read.csv(file="/Documents/2024/E_vs_I_Liver/E_vs_I_Liver.csv",header =T, sep = "," )
sample_annotation <- read.csv(file = "/Documents/2024/E_vs_I_Liver/E_vs_I_Liver_anno.csv", header =T, sep = ",")
contrasts=c("condition","Grp_E_Liver", "Grp_I_Liver")
# Run function
jen_deseq2_auto_function(count_matrix,sample_annotation,contrasts,outputPrefix, outputdirectory)


######################################
# Liver all samples PCA plot 
######################################

count_matrix1 <- read.csv(file = '/Documents/2024/part1_Liver_counts_matrix.csv', header =T, sep = ",")
row.names(count_matrix1) <- count_matrix1$gene_id
count_matrix1 <- subset(count_matrix1, select = -c(gene_id))
# # In case of RSEM, need to round values
# count_matrix1 <- round(count_matrix1)
sample_annotation1 <- read.csv("/Documents/2024/part1_Liver_anno.csv", header =T, sep = ",")
rownames(sample_annotation1) <- sample_annotation1$sample
sample_annotation1 <- subset(sample_annotation1, select = -c(sample))
all(rownames(sample_annotation1) == colnames(count_matrix1))

dds1 <- DESeqDataSetFromMatrix(countData = count_matrix1,
                               colData = sample_annotation1, 
                               design = ~ condition)
dds_deg1 <- DESeq(dds1)
resultsNames(dds_deg1)
contrasts = c("condition","treated","untreated")
res <- results(dds_deg1, contrast=contrasts)
print(head(res))
print(summary(res))
######################################
# PCA plot 
######################################

dds_vst <- vst(dds1, blind = T)
pcaplot_vst<- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1.5) + geom_point(size=0.1)
pcaplot_vst
outputdirectory <- "/Documents/2024/"
ggsave(paste0(outputdirectory, "PCA_plot_vst_Liver_treated_vs_untreated_1.png"), pcaplot_vst, height = 10, width = 10, dpi = 320) 

######################################
# Heart all samples PCA plot 
######################################

count_matrix1 <- read.csv(file = '/Documents/2024/part1_Heart_counts_matrix.csv', header =T, sep = ",")
row.names(count_matrix1) <- count_matrix1$gene_id
count_matrix1 <- subset(count_matrix1, select = -c(gene_id))
# # In case of RSEM, need to round values
# count_matrix1 <- round(count_matrix1)
sample_annotation1 <- read.csv("/Documents/2024/part1_Heart_anno.csv", header =T, sep = ",")
rownames(sample_annotation1) <- sample_annotation1$sample
sample_annotation1 <- subset(sample_annotation1, select = -c(sample))
all(rownames(sample_annotation1) == colnames(count_matrix1))

dds1 <- DESeqDataSetFromMatrix(countData = count_matrix1,
                               colData = sample_annotation1, 
                               design = ~ condition)
dds_deg1 <- DESeq(dds1)
resultsNames(dds_deg1)
contrasts = c("condition","treated","untreated")
res <- results(dds_deg1, contrast=contrasts)
print(head(res))
print(summary(res))
######################################
# PCA plot 
######################################

dds_vst <- vst(dds1, blind = T)
pcaplot_vst<- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1.5) + geom_point(size=0.1)
pcaplot_vst
outputdirectory <- "/Documents/2024/"
#ggsave(paste0(outputdirectory, "PCA_plot_vst_Liver_treated_vs_untreated.png"), pcaplot_vst, height = 10, width = 10, dpi = 320) 
ggsave(paste0(outputdirectory, "PCA_plot_vst_Heart_treated_vs_untreated1.png"), pcaplot_vst, height = 10, width = 10, dpi = 320) 
