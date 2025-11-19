#!/usr/bin/env Rscript

# DESCRIPTION
# - Calculate EED KO score via gene module scoring with Seurat

# INPUT
# - Seurat object where cell cycle is regressed out
# - List of EED KO gene module genes, that are upregulated upon EED KO in K562
#   but robust to Interferon stimulation

# OUTPUT
# - Seurat object with EED KO score

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data

be_singlets <-
  qs::qread(
    paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides.qs")
  )

# Robust gene module of 84 genes that are upregulated upon EED KO in K562 but robust to Interferon stimulation,
# from intersection of Marian's and Christina's data
eed_ko_gene_module_robust <- read_lines(
  paste0(result_dir, "christina_marian_eed_ko_gene_module.txt")
)

## Calculate EED KO score

# Calculate Seurat module scores to see which cells are most similar to EED KO phenotype
be_singlets <-
  AddModuleScore(
    be_singlets,
    features = list(eed_ko_gene_module_robust = eed_ko_gene_module_robust),
    assay = "RNA",
    name = c("eed_ko_gene_module_robust_")
  )

# Normalize EED KO score by the control sample, so that it makes sense that Ctrl is around 0 and other samples are higher.
# Subtract the mean EED KO score of the Ctrl sample.
# This is inspired by Dixit et al.
mean_eed_ko_Ctrl_HLA_pos <- mean(subset(be_singlets, sample == "Ctrl-HLA-pos")$eed_ko_gene_module_robust_1)

be_singlets$eed_ko_score_normalized <- be_singlets$eed_ko_gene_module_robust_1 - mean_eed_ko_Ctrl_HLA_pos

qs::qsave(
  be_singlets,
  paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs")
)
