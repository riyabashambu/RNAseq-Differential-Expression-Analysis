################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 06_go_enrichment.R
#
# Author  : Riya
#
# Purpose :
# Perform Gene Ontology (GO) and KEGG pathway enrichment analysis
# using clusterProfiler.
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
library(clusterProfiler)
library(org.Dm.eg.db)
library(AnnotationDbi)

################################################################################
# Load DESeq2 Results
################################################################################

load("results/RData/deseq2_results.RData")

################################################################################
# Strict Differentially Expressed Genes
################################################################################

res_df <- as.data.frame(resOrdered)

res_df$GeneID <- rownames(res_df)

sig_genes <- subset(
  res_df,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

################################################################################
# Convert FlyBase IDs to Entrez IDs
################################################################################

gene_entrez <- bitr(
  sig_genes$GeneID,
  fromType = "FLYBASE",
  toType = "ENTREZID",
  OrgDb = org.Dm.eg.db
)

################################################################################
# Create Background Gene Set
################################################################################

background_conversion <- bitr(
  rownames(dds),
  fromType = "FLYBASE",
  toType = "ENTREZID",
  OrgDb = org.Dm.eg.db
)

background_genes <- unique(background_conversion$ENTREZID)

################################################################################
# Gene Ontology Enrichment (Strict Threshold)
################################################################################

go_result <- enrichGO(
  gene = unique(gene_entrez$ENTREZID),
  universe = background_genes,
  OrgDb = org.Dm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(go_result),
  "results/go/GO_Strict_Results.csv",
  row.names = FALSE
)

################################################################################
# Exploratory GO Enrichment
################################################################################

res_relaxed <- results(dds, alpha = 0.10)

sig_genes_relaxed <- subset(
  as.data.frame(res_relaxed),
  padj < 0.10 &
    abs(log2FoldChange) > 0.58
)

flybase_relaxed <- rownames(sig_genes_relaxed)

gene_entrez_relaxed <- bitr(
  flybase_relaxed,
  fromType = "FLYBASE",
  toType = "ENTREZID",
  OrgDb = org.Dm.eg.db
)

go_result_relaxed <- enrichGO(
  gene = unique(gene_entrez_relaxed$ENTREZID),
  universe = background_genes,
  OrgDb = org.Dm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.10,
  qvalueCutoff = 0.50,
  readable = TRUE
)

write.csv(
  as.data.frame(go_result_relaxed),
  "results/go/GO_Exploratory_Results.csv",
  row.names = FALSE
)

################################################################################
# KEGG Pathway Enrichment
################################################################################

kegg_result <- enrichKEGG(
  gene = unique(gene_entrez_relaxed$ENTREZID),
  organism = "dme",
  universe = background_genes,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.10
)

write.csv(
  as.data.frame(kegg_result),
  "results/go/KEGG_Results.csv",
  row.names = FALSE
)

################################################################################
# Summary
################################################################################

cat("\n")
cat("=====================================================\n")
cat(" Functional Enrichment Analysis Completed\n")
cat("=====================================================\n")
cat("Strict GO terms        :", nrow(as.data.frame(go_result)), "\n")
cat("Exploratory GO terms   :", nrow(as.data.frame(go_result_relaxed)), "\n")
cat("KEGG pathways detected :", nrow(as.data.frame(kegg_result)), "\n")
cat("Results saved to       : results/go/\n")
cat("=====================================================\n")

################################################################################
# Save GO Objects
################################################################################

save(
  go_result,
  go_result_relaxed,
  kegg_result,
  file = "results/RData/go_results.RData"
)