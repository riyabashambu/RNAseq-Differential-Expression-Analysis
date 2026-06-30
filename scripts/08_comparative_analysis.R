################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 08_comparative_analysis.R
#
# Purpose:
# Compare DEGs identified by DESeq2, edgeR and limma-voom
# Generate Venn Diagram
# Generate UpSet Plot
################################################################################

rm(list = ls())

options(scipen = 999)

################################################################################
# Load Packages
################################################################################

library(VennDiagram)
library(UpSetR)

################################################################################
# Create Output Directory
################################################################################

dir.create("plots", showWarnings = FALSE)

################################################################################
# Load Results
################################################################################

deseq_results <- read.csv(
  "results/deseq2/DESeq2_Annotated_Results.csv"
)

edge_results <- read.csv(
  "results/edgeR/edgeR_All_Results.csv",
  row.names = 1
)

limma_results <- read.csv(
  "results/limma_voom/limma_voom_All_Results.csv",
  row.names = 1
)

################################################################################
# Extract Significant Genes
################################################################################

deseq_genes <- subset(
  deseq_results,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)$GeneID

edge_genes <- rownames(
  subset(
    edge_results,
    FDR < 0.05 &
      abs(logFC) > 1
  )
)

limma_genes <- rownames(
  subset(
    limma_results,
    adj.P.Val < 0.05 &
      abs(logFC) > 1
  )
)

################################################################################
# DEG Summary
################################################################################

deg_summary <- data.frame(
  
  Method = c(
    "DESeq2",
    "edgeR",
    "limma-voom"
  ),
  
  Significant_DEGs = c(
    length(deseq_genes),
    length(edge_genes),
    length(limma_genes)
  )
  
)

write.csv(
  deg_summary,
  "results/DEG_Summary.csv",
  row.names = FALSE
)

################################################################################
# Venn Diagram
################################################################################

venn.diagram(
  
  x = list(
    
    DESeq2 = deseq_genes,
    edgeR = edge_genes,
    limma = limma_genes
    
  ),
  
  filename = "plots/DEG_Venn.png",
  
  fill = c(
    "#F8766D",
    "#619CFF",
    "#00BA38"
  ),
  
  alpha = 0.5,
  
  cex = 1.4,
  
  cat.cex = 1.3,
  
  cat.pos = c(-20,20,180)
  
)

################################################################################
# UpSet Plot
################################################################################

gene_lists <- list(
  
  DESeq2 = deseq_genes,
  
  edgeR = edge_genes,
  
  limma = limma_genes
  
)

png(
  "plots/UpSet_Plot.png",
  width = 2200,
  height = 1600,
  res = 300
)

upset(
  
  fromList(gene_lists),
  
  order.by = "freq",
  
  main.bar.color = "steelblue",
  
  sets.bar.color = "gray40",
  
  text.scale = 1.5
  
)

dev.off()

################################################################################
# Shared Genes
################################################################################

shared_genes <- Reduce(
  
  intersect,
  
  list(
    deseq_genes,
    edge_genes,
    limma_genes
  )
  
)

write.csv(
  
  data.frame(
    GeneID = shared_genes
  ),
  
  "results/Shared_DEGs_All_Methods.csv",
  
  row.names = FALSE
  
)

################################################################################
# Finished
################################################################################

cat("\n=========================================\n")
cat("Comparative Analysis Completed\n")
cat("=========================================\n")

cat("DESeq2 DEGs :", length(deseq_genes), "\n")
cat("edgeR DEGs  :", length(edge_genes), "\n")
cat("limma DEGs  :", length(limma_genes), "\n")
cat("Shared DEGs :", length(shared_genes), "\n")

cat("\nResults saved to:\n")
cat("results/\n")
cat("plots/\n")