# Installation of the latest released version
install.packages('GOplot')
library(GOplot)
# # Installation of the latest development version
# install_github('wencke/wencke.github.io')
library(cowplot)


# Read gpofiler data
functional_enrichment_data_DHvsH = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/gProfiler_PM__intersections DH vs MH.csv", sep = ",")
functional_enrichment_data_DHvsH <- functional_enrichment_data_DHvsH[c(1,2,3,4,10)]
functional_enrichment_data_reorder_DHvsH <- functional_enrichment_data_DHvsH[, c("source","term_id","term_name","intersections", "adjusted_p_value")]
names(functional_enrichment_data_reorder_DHvsH) <- c("Category","ID", "Term", "Genes", "adj_pval")
# Read Gene DEG data
deg_data_DHvsH = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/DHvsH_deg2_go_plots.csv",sep = ",", header = TRUE)
#deg_data_DHvsH <- deg_data_DHvsH[c(2,4,5,6)]
names(deg_data_DHvsH) <- c("ID","logFC", "p_value", "adj_pvalue")
circ_DHvsH <- circle_dat(functional_enrichment_data_reorder_DHvsH, deg_data_DHvsH)
GOBubble(circ_DHvsH, labels = 1, table.legend = F, ID = F, display = "single")
# reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# GOBubble(reduced_circ, labels = 2.8)
# save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/DHvsH_GO_bubbleplot.png",GOBubble(circ_DHvsH, labels = 1, table.legend = F, ID = F, display = "single"))
# GOCircle(circ)
# save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/DHvsH_GO_circ_plot.png",GOCircle(circ))
# # Generate a circular visualization for 10 terms
GOCircle(circ_DHvsH, nsub =10)

# Read gpofiler data
functional_enrichment_data_CvsL = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/gProfiler_C vs L.csv", sep = ",")
functional_enrichment_data_CvsL <- functional_enrichment_data_CvsL[c(1,2,3,4,10)]
functional_enrichment_data_reorder_CvsL <- functional_enrichment_data_CvsL[, c("source","term_id","term_name","intersections", "adjusted_p_value")]
names(functional_enrichment_data_reorder_CvsL) <- c("Category","ID", "Term", "Genes", "adj_pval")
# Read Gene DEG data
deg_data_CvsL = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/CvsL_foldchange_DEG.csv",sep = ",", header = TRUE)
deg_data_CvsL <- deg_data_CvsL[c(2,4,5,6)]
names(deg_data_CvsL) <- c("ID","logFC", "p_value", "adj_pvalue")
circ_CvsL <- circle_dat(functional_enrichment_data_reorder_CvsL, deg_data_CvsL)
GOBubble(circ_CvsL, labels = 1, table.legend = T, ID = T, display = "single")
# reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# GOBubble(reduced_circ, labels = 2.8)
# save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/CvsL_GO_bubbleplot.png",GOBubble(circ_CvsL, labels = 1, table.legend = T, ID = T, display = "single"))


# EC$genes
# EC$process
# head(EC$david)
# head(EC$genelist)
# circ_test <- circle_dat(EC$david, EC$genelist)
# GOBubble(circ_test, labels = 3, table.legend = T, ID = T, display = "single")

go_plot <- function(infile1, infile2, outfilepath,outname){
  # Read gpofiler data
  functional_enrichment_data = infile1
  functional_enrichment_data <- functional_enrichment_data[c(1,2,3,4,10)]
  functional_enrichment_data_reorder <- functional_enrichment_data[, c("source","term_id","term_name","intersections", "adjusted_p_value")]
  names(functional_enrichment_data_reorder) <- c("Category","ID", "Term", "Genes", "adj_pval")
  # Read Gene DEG data
  deg_data = infile2
  deg_data <- deg_data[c(2,4,5,6)]
  names(deg_data) <- c("ID","logFC", "p_value", "adj_pvalue")
  circ <- circle_dat(functional_enrichment_data_reorder, deg_data)
  print(circ$zscore)
  bubble_plot <- GOBubble(circ, labels = 0.5)
  # bubble_plot
  #save_plot(paste(outfilepath,outname,".png",  sep = ""),bubble_plot)
}
infile1 = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/gProfiler_C vs L.csv", sep = ",")
infile2 = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/CvsL_foldchange_DEG.csv",sep = ",", header = TRUE)
outname = "CvsL_Bubble_plot"
outfilepath = "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/GO_Plots/"
go_plot(infile1, infile2, outfilepath, outname)
