################################################################################
# limma-voom Differential Expression Analysis
################################################################################

# Create a DGEList object from the filtered count matrix
dge_limma <- DGEList(
  counts = counts(dds),
  group = colData(dds)$condition
)

# Calculate TMM normalization factors
dge_limma <- calcNormFactors(dge_limma)

# Construct the design matrix
design_limma <- model.matrix(
  ~ condition,
  data = colData(dds)
)

# Transform count data using the voom method
voom_data <- voom(
  dge_limma,
  design_limma,
  plot = FALSE
)

# Fit the linear model
fit <- lmFit(
  voom_data,
  design_limma
)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract differential expression results
limma_results <- topTable(
  fit,
  coef = 2,
  number = Inf,
  sort.by = "P"
)

# Identify significantly differentially expressed genes
limma_sig <- subset(
  limma_results,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

################################################################################
# Save limma-voom Results
################################################################################

write.csv(
  limma_results,
  "results/limma/limma_All_Results.csv",
  row.names = TRUE
)

write.csv(
  limma_sig,
  "results/limma/limma_Significant_Genes.csv",
  row.names = TRUE
)

################################################################################
# limma-voom Analysis Summary
################################################################################

cat("\n")
cat("=====================================================\n")
cat("          limma-voom Differential Expression\n")
cat("=====================================================\n")
cat("Total genes analysed      :", nrow(limma_results), "\n")
cat("Significant DEGs (FDR<0.05, |log2FC|>1):", nrow(limma_sig), "\n")
cat("Results saved to          : results/limma/\n")
cat("=====================================================\n\n")