#!/usr/bin/env Rscript

# DESCRIPTION
# Preprocess scRNA-seq data of MF03 experiment
# - Combine data from 2 captures
# - Keep genes that are expressed in at least 5 cells
# - Keep cells with more than 1700 and less than 8000 genes expressed, more than 5000 UMI and less than 8% mitochondrial counts
# - Normalize data and log1p transform
# - Demultiplex samples and label doublets by LMOs with positive.quantile=0.99
# NB: Cells were super loaded, so doublet percentage of 30% expected.

# INPUT
# - Filtered count matrices of Cell Ranger results
# - Cite-seq-count LMO count matrix
# - List of samples and LMOs

# OUTPUT
# - Seurat object with filtered normalized RNA counts and demultiplexed cells (sample annotation, doublets, negatives)

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data
cat("Loading data\n")
# Load count matrix from Cell Ranger results
be_p1_data <-
  Read10X(
    paste0(result_dir, "MF03/scRNA_seq/cellranger_results/MF03_P1/outs/filtered_feature_bc_matrix/")
  )
be_p2_data <-
  Read10X(
    paste0(result_dir, "MF03/scRNA_seq/cellranger_results/MF03_P2/outs/filtered_feature_bc_matrix/")
  )

# Read in LMO counts from Cite-seq-count results.
lmo_p1 <-
  Read10X(
    paste0(result_dir, "MF03/scRNA_seq/lmo_results/MF03_P1/umi_count/"),
    gene.column = 1
  )
lmo_p2 <-
  Read10X(
    paste0(result_dir, "MF03/scRNA_seq/lmo_results/MF03_P2/umi_count/"),
    gene.column = 1
  )

# Read in table with information which LMO is for which sample
lmo_table <-
  read_tsv(paste0(result_dir, "MF03/scRNA_seq/lmo_list_MF03.tsv"))

## QC and normalization
cat("QC and normalization\n")
### QC

# Initialize Seurat object
be_p1 <- CreateSeuratObject(counts = be_p1_data)
be_p1$pool <- "P1"

be_p2 <- CreateSeuratObject(counts = be_p2_data)
be_p2$pool <- "P2"

# Merge captures, assuming no batch effect
be <- merge(be_p1, y = be_p2, add.cell.ids = c("P1", "P2"))
rm(be_p1_data, be_p2_data)

# Remove genes expressed in less than 5 cells across both pools.
# Doing this here instead of when creating Seurat object to keep genes that are
# potentially only detected in one pool/capture.
be <- be[rowSums(be@assays$RNA@counts) >= 5, ]

# Calculate percentage of reads that map to mitochondrial genome.
# This is an indicator of broken and dead cells.
be[["percent.mt"]] <- PercentageFeatureSet(be, pattern = "^MT-")

# Filter cells based on QC metrics

be <-
  subset(be,
         subset = nFeature_RNA > 1700 &
           nFeature_RNA < 8000 & nCount_RNA > 5000 & percent.mt < 8)

### Normalization

# Normalize cells by sequencing depth per cell and log transform with pseudo count.

be <- NormalizeData(be, normalization.method = "LogNormalize", scale.factor = 10000)

## Demultiplex
cat("Demultiplexing\n")
### Read LMO data

lmo_table <- lmo_table %>%
  mutate(samplename = paste(Name, Sample, sep = "_"))

# Remove row with unknown LMOs, convert rownames to sample names.
# Add pool identifier to column names so that cell names are the same as in Seurat object.
lmo_p1 <- lmo_p1[c(1:8), ]
rownames(lmo_p1) <- lmo_table$samplename
colnames(lmo_p1) <- paste0("P1_", colnames(lmo_p1), "-1")

lmo_p2 <- lmo_p2[c(1:8), ]
rownames(lmo_p2) <- lmo_table$samplename
colnames(lmo_p2) <- paste0("P2_", colnames(lmo_p2), "-1")

# Join the LMO counts of the captures

lmo <- cbind(lmo_p1, lmo_p2)

### Add LMO data to Seurat object

# Select cell barcodes detected by both RNA and LMO.
joint_cells <- intersect(colnames(be), colnames(lmo))

# Since we run Cite-seq-count with a whitelist of all cells detected by cell ranger,
# all cells with RNA counts will have some LMO counts.
# There will be additional cells with LMO counts though.
# Subset both datasets to shared cells.
be <- be[, joint_cells]
lmo <- lmo[, joint_cells]
lmo <- as.matrix(lmo)

# Add LMO counts as assay to Seurat object.
be[["LMO"]] <- CreateAssayObject(counts = lmo)

# Normalize LMO counts
be <- NormalizeData(be, assay = "LMO", normalization.method = "CLR")

### Demultiplex

be <- HTODemux(be, assay = "LMO", positive.quantile = 0.99)

Idents(be) <- "hash.ID"

be@meta.data <- be@meta.data %>%
  mutate(sample = gsub("BC[0-9]+-", "", hash.ID))

cat("Saving results\n")

qs::qsave(
  be,
  paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed.qs")
)
