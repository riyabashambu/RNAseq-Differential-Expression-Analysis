################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 01_data_preprocessing.R
#
# Author  : Riya
#
# Purpose :
# Import the Pasilla RNA-seq dataset, validate sample metadata,
# perform basic quality assessment, filter low-count genes,
# and prepare the count matrix for downstream analyses.
#
# Dataset :
# Pasilla RNA-seq Dataset
# Organism: Drosophila melanogaster
#
# Input :
# Pasilla Bioconductor experiment package
#
# Output :
# - Filtered count matrix
# - Sample metadata
# - DESeq2 dataset object
# - Preprocessed objects saved for downstream analysis
################################################################################

################################################################################
# Clear Environment
################################################################################

rm(list = ls())

options(scipen = 999)

################################################################################
# Load Required Packages
################################################################################

library(DESeq2)
library(pasilla)

################################################################################
# Import Pasilla Dataset
################################################################################

pasilla_dir <- system.file("extdata", package = "pasilla")

counts <- read.delim(
  file.path(pasilla_dir, "pasilla_gene_counts.tsv"),
  row.names = 1
)

coldata <- read.csv(
  file.path(pasilla_dir, "pasilla_sample_annotation.csv"),
  row.names = 1
)

################################################################################
# Metadata Formatting
################################################################################

rownames(coldata) <- gsub("fb", "", rownames(coldata))

coldata <- coldata[colnames(counts), ]

################################################################################
# Data Validation
################################################################################

cat("\n==============================\n")
cat("DATA VALIDATION\n")
cat("==============================\n\n")

cat("Count Matrix Dimensions:\n")
print(dim(counts))

cat("\nMetadata Dimensions:\n")
print(dim(coldata))

cat("\nSample Names Match:\n")
print(all(colnames(counts) == rownames(coldata)))

cat("\nMissing Values in Count Matrix:\n")
print(sum(is.na(counts)))

cat("\nMissing Values in Metadata:\n")
print(sum(is.na(coldata)))

################################################################################
# Construct DESeq2 Dataset
################################################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

################################################################################
# Filter Lowly Expressed Genes
################################################################################

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

################################################################################
# Preprocessing Summary
################################################################################

cat("\n==============================\n")
cat("PREPROCESSING SUMMARY\n")
cat("==============================\n\n")

cat("Genes before filtering :",
    nrow(counts), "\n")

cat("Genes after filtering  :",
    nrow(dds), "\n")

cat("Samples               :",
    ncol(dds), "\n")

################################################################################
# Create Output Directory
################################################################################

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
# Save Processed Objects
################################################################################

saveRDS(
  dds,
  file = "data/processed/dds_filtered.rds"
)

saveRDS(
  counts(dds),
  file = "data/processed/filtered_counts.rds"
)

write.csv(
  as.data.frame(colData(dds)),
  "data/processed/sample_metadata.csv",
  row.names = TRUE
)

################################################################################
# Completion Message
################################################################################

cat("\n=========================================\n")
cat("Data preprocessing completed successfully.\n")
cat("Processed files saved to data/processed/\n")
cat("=========================================\n")