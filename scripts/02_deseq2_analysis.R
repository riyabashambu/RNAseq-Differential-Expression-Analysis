################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 02_deseq2_analysis.R
#
# Author  : Riya
#
# Purpose :
# Perform differential expression analysis using DESeq2.
# Apply log2 fold-change shrinkage, identify significant genes,
# and export normalized counts and DEG tables.
################################################################################

################################################################################
# Clear Environment
################################################################################

rm(list = ls())

options(scipen = 999)

################################################################################
# Load Packages
################################################################################

library(DESeq2)
library(apeglm)

################################################################################
# Load Preprocessed Dataset
################################################################################

dds <- readRDS("data/processed/dds_filtered.rds")

################################################################################
# Run DESeq2 Analysis
################################################################################

dds <- DESeq(dds)

################################################################################
# Extract Results
################################################################################

res <- results(dds)

################################################################################
# Log2 Fold Change Shrinkage
################################################################################

res_shrunk <- lfcShrink(
  dds,
  coef = 2,
  type = "apeglm"
)

################################################################################
# Order Results
################################################################################

res_ordered <- res_shrunk[order(res_shrunk$padj), ]

################################################################################
# Convert to Data Frame
################################################################################

res_df <- as.data.frame(res_ordered)

res_df$GeneID <- rownames(res_df)

################################################################################
# Significant DEGs
################################################################################

sig_res <- subset(
  res_df,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

################################################################################
# Normalized Counts
################################################################################

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

################################################################################
# Create Output Directory
################################################################################

dir.create(
  "results/deseq2",
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
# Export Results
################################################################################

write.csv(
  res_df,
  "results/deseq2/DESeq2_All_Results.csv",
  row.names = FALSE
)

write.csv(
  sig_res,
  "results/deseq2/DESeq2_Significant_DEGs.csv",
  row.names = FALSE
)

write.csv(
  normalized_counts,
  "results/deseq2/Normalized_Counts.csv"
)

################################################################################
# Summary Report
################################################################################

sink("results/deseq2/DESeq2_Summary.txt")

cat("DESeq2 Differential Expression Analysis\n")
cat("--------------------------------------\n\n")

cat("Total genes analysed :",
    nrow(res_df), "\n")

cat("Significant DEGs :",
    nrow(sig_res), "\n\n")

summary(res_shrunk)

sink()

################################################################################
# Save Objects
################################################################################

saveRDS(
  dds,
  "results/deseq2/dds_deseq2.rds"
)

saveRDS(
  res_shrunk,
  "results/deseq2/deseq2_results.rds"
)

################################################################################
# Completion Message
################################################################################

cat("\n========================================\n")
cat("DESeq2 analysis completed successfully.\n")
cat("Results saved in results/deseq2/\n")
cat("========================================\n")