################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 03_edger_analysis.R
#
# Author  : Riya
#
# Purpose :
# Perform differential expression analysis using edgeR.
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
library(edgeR)

################################################################################
# Load Preprocessed Dataset
################################################################################

dds <- readRDS("data/processed/dds_filtered.rds")

################################################################################
# Create edgeR Object
################################################################################

group <- colData(dds)$condition

dge <- DGEList(
  counts = counts(dds),
  group = group
)

################################################################################
# Filter Lowly Expressed Genes
################################################################################

keep <- filterByExpr(dge)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

################################################################################
# Library Size Normalization
################################################################################

dge <- calcNormFactors(dge)

################################################################################
# Design Matrix
################################################################################

design <- model.matrix(~group)

################################################################################
# Estimate Dispersion
################################################################################

dge <- estimateDisp(
  dge,
  design
)

################################################################################
# Differential Expression Analysis
################################################################################

et <- exactTest(
  dge,
  pair = c("treated", "untreated")
)

################################################################################
# Extract Results
################################################################################

res <- topTags(
  et,
  n = Inf
)

res_df <- res$table

res_df$GeneID <- rownames(res_df)

################################################################################
# Significant Genes
################################################################################

sig_res <- subset(
  res_df,
  FDR < 0.05 &
    abs(logFC) > 1
)

################################################################################
# Normalized Counts
################################################################################

norm_counts <- cpm(
  dge,
  normalized.lib.sizes = TRUE
)

################################################################################
# Create Output Directory
################################################################################

dir.create(
  "results/edger",
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
# Export Results
################################################################################

write.csv(
  res_df,
  "results/edger/edgeR_All_Results.csv",
  row.names = FALSE
)

write.csv(
  sig_res,
  "results/edger/edgeR_Significant_DEGs.csv",
  row.names = FALSE
)

write.csv(
  norm_counts,
  "results/edger/edgeR_Normalized_Counts.csv"
)

################################################################################
# Summary Report
################################################################################

sink(
  "results/edger/edgeR_Summary.txt"
)

cat("edgeR Differential Expression Analysis\n")
cat("--------------------------------------\n\n")

cat("Total genes analysed :",
    nrow(res_df),
    "\n")

cat("Significant DEGs :",
    nrow(sig_res),
    "\n\n")

summary(decideTestsDGE(et))

sink()

################################################################################
# Save Objects
################################################################################

saveRDS(
  dge,
  "results/edger/edgeR_DGEList.rds"
)

saveRDS(
  et,
  "results/edger/edgeR_Results.rds"
)

################################################################################
# Completion Message
################################################################################

cat("\n========================================\n")
cat("edgeR analysis completed successfully.\n")
cat("Results saved in results/edger/\n")
cat("========================================\n")

dir.create("results/RData",
           recursive = TRUE,
           showWarnings = FALSE)

save(
  dds,
  res,
  resOrdered,
  file = "results/RData/deseq2_results.RData"
)
