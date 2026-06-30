################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 05_annotation.R
#
# Author  : Riya
#
# Purpose :
# Annotate DESeq2 differential expression results with FlyBase IDs,
# gene symbols, and export complete and significant gene tables.
################################################################################

#------------------------------------------------------------------------------
# Load Required Packages
#------------------------------------------------------------------------------

library(DESeq2)
library(org.Dm.eg.db)
library(AnnotationDbi)

#------------------------------------------------------------------------------
# Load Objects from Previous Analysis
#------------------------------------------------------------------------------

load("results/RData/deseq2_results.RData")

# This file should contain:
# dds
# resOrdered

#------------------------------------------------------------------------------
# Convert DESeq2 Results to Data Frame
#------------------------------------------------------------------------------

res_df <- as.data.frame(resOrdered)

res_df$GeneID <- rownames(res_df)

#------------------------------------------------------------------------------
# Annotate FlyBase IDs with Gene Symbols
#------------------------------------------------------------------------------

res_df$GeneSymbol <- mapIds(
  org.Dm.eg.db,
  keys = res_df$GeneID,
  column = "SYMBOL",
  keytype = "FLYBASE",
  multiVals = "first"
)

#------------------------------------------------------------------------------
# Reorder Columns
#------------------------------------------------------------------------------

res_df <- res_df[, c(
  "GeneID",
  "GeneSymbol",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

#------------------------------------------------------------------------------
# Export Complete Annotated Results
#------------------------------------------------------------------------------

write.csv(
  res_df,
  "results/deseq2/DESeq2_Annotated_Results.csv",
  row.names = FALSE
)

#------------------------------------------------------------------------------
# Extract Significant Genes
#------------------------------------------------------------------------------

sig_genes <- subset(
  res_df,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

#------------------------------------------------------------------------------
# Export Significant Genes
#------------------------------------------------------------------------------

write.csv(
  sig_genes,
  "results/deseq2/Significant_Genes.csv",
  row.names = FALSE
)

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat("========================================\n")
cat("Annotation Complete\n")
cat("========================================\n")
cat("Total genes :", nrow(res_df), "\n")
cat("Significant genes :", nrow(sig_genes), "\n")
cat("Results saved to:\n")
cat("results/deseq2/\n")
cat("========================================\n")