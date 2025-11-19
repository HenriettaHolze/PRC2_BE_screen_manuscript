#!/usr/bin/env Rscript

# DESCRIPTION
# - Remove doublets identified by LMO demultiplexing and cells without LMO annotation
# - Identify highly variable genes and scale data
# - Cell cycle annotation and regression
# - Regress out cell cycle signature
# - Calcualte PCA, use first 20 PCs for UMAP embedding and Louvain clustering, and 20 nearest neighbours

# INPUT
# - Seurat object with filtered normalized RNA counts and demultiplexed cells (sample annotation, doublets, negatives)

# OUTPUT
# - Seurat object with cc regressed, UMAP reduction and Louvain clustering

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data

cat("Loading data\n")
# Load filtered, normalized and demultiplexed scRNA-seq data.
be <- qs::qread(
  paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed.qs")
)

# Remove doublets (close to 30% of cells) and cells without LMO annotation (likely low quality cells)
be_singlets <- subset(be, LMO_classification.global != "Doublet")
be_singlets <- subset(be_singlets, LMO_classification.global != "Negative")

# Identify highly variable features to use for dimensionality reduction.
be_singlets <- FindVariableFeatures(be_singlets,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Scale gene expression
be_singlets <- ScaleData(be_singlets)

## Cell cycle regression
cat("Cell cycle scoring\n")
# Cell cycle scoring, using Seurat's cell cycle genes
be_singlets <-
  CellCycleScoring(
    be_singlets,
    s.features = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )

# Regress out cell cycle signature
be_singlets <- ScaleData(be_singlets, vars.to.regress = c("S.Score", "G2M.Score"))

## Dimensionality reduction
cat("Dimensionality reduction\n")
# Calculate PCA on scaled data.
be_singlets <- RunPCA(be_singlets, features = VariableFeatures(object = be_singlets))

# Create shared nearest neighbour graph using first 20 principle components and 20 nearest neighbours
be_singlets <- FindNeighbors(be_singlets, dims = 1:20, k.param = 20)

# Perform Louvain clustering
be_singlets <-
  FindClusters(be_singlets, resolution = c(0.2, 0.3, 0.4, 0.5, 0.6))

# Create UMAP embedding using same first 20 PCs and 20 NN
be_singlets <-
  RunUMAP(be_singlets, dims = 1:20, n.neighbors = 20)

cat("Saving results\n")
qs::qsave(
  be_singlets,
  paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap.qs")
)
