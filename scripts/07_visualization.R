################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 07_visualization.R
#
# Author  : Riya
#
# Purpose :
# Generate publication-quality figures for RNA-seq differential
# expression analysis.
################################################################################

#------------------------------------------------------------------------------
# Clear Environment
#------------------------------------------------------------------------------

rm(list = ls())

options(scipen = 999)

################################################################################
# Load Required Packages
################################################################################

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(enrichplot)
library(ggplot2)
library(grid)

################################################################################
# Load Analysis Results
################################################################################

load("results/RData/deseq2_results.RData")
load("results/RData/go_results.RData")

################################################################################
# Create Output Directory
################################################################################

dir.create("plots", recursive = TRUE, showWarnings = FALSE)

################################################################################
# Variance Stabilizing Transformation
################################################################################

vsd <- vst(dds, blind = FALSE)

################################################################################
# PCA Plot
################################################################################

png(
  "plots/PCA_Plot.png",
  width = 1800,
  height = 1600,
  res = 300
)

plotPCA(
  vsd,
  intgroup = "condition"
)

dev.off()

################################################################################
# Sample Distance Heatmap
################################################################################

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

png(
  "plots/Sample_Distance_Heatmap.png",
  width = 2000,
  height = 1800,
  res = 300
)

hm <- pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  main = "Sample Distance Heatmap"
)

grid.newpage()
grid.draw(hm$gtable)

dev.off()

################################################################################
# MA Plot
################################################################################

png(
  "plots/MA_Plot.png",
  width = 1800,
  height = 1600,
  res = 300
)

plotMA(
  res,
  ylim = c(-5,5)
)

dev.off()

################################################################################
# Volcano Plot
################################################################################

top_genes <- head(
  rownames(res)[order(res$padj)],
  10
)

png(
  "plots/Volcano_Plot.png",
  width = 2400,
  height = 2000,
  res = 300
)

EnhancedVolcano(
  
  res,
  
  lab = rownames(res),
  
  selectLab = top_genes,
  
  x = "log2FoldChange",
  
  y = "padj",
  
  pCutoff = 0.05,
  
  FCcutoff = 1,
  
  title = "Differential Gene Expression",
  
  subtitle = "Treated vs Untreated",
  
  caption = "DESeq2",
  
  labSize = 4,
  
  pointSize = 2,
  
  drawConnectors = TRUE,
  
  widthConnectors = 0.5,
  
  colConnectors = "grey30",
  
  arrowheads = FALSE
  
)

dev.off()

################################################################################
# GO Dot Plot
################################################################################

png(
  "plots/GO_Dotplot_Exploratory.png",
  width = 2400,
  height = 2000,
  res = 300
)

dotplot(
  go_result_relaxed,
  showCategory = 15
) +
  ggtitle(
    "GO Biological Process Enrichment\n(Exploratory Analysis)"
  )

dev.off()

################################################################################
# GO Bar Plot
################################################################################

png(
  "plots/GO_Barplot_Exploratory.png",
  width = 2400,
  height = 2000,
  res = 300
)

barplot(
  go_result_relaxed,
  showCategory = 15,
  title = "GO Biological Process Enrichment\n(Exploratory Analysis)"
)

dev.off()

################################################################################
# Summary
################################################################################

cat("\n")
cat("=================================================\n")
cat("Visualization Complete\n")
cat("=================================================\n")
cat("Generated Figures:\n")
cat("  PCA Plot\n")
cat("  Sample Distance Heatmap\n")
cat("  MA Plot\n")
cat("  Volcano Plot\n")
cat("  GO Dot Plot\n")
cat("  GO Bar Plot\n")
cat("=================================================\n")