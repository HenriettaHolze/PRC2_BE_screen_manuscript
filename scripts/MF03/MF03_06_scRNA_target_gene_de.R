#!/usr/bin/env Rscript

# DESCRIPTION
# Calculate differential expression of the target genes for each guide against the Ctrl HLA+ sample.
# Do these for all samples.

# INPUT
# - Seurat object with log normalized gene expression
# - Guide annotation per cell

# OUTPUT
# Table with differential expression results.

library(Seurat)
library(tidyverse)
library(ltm)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data
be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )

guide_counts_single_integration <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_single_guide.tsv"
))
guide_counts_aggregated <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_aggregated.tsv"
))

## wilcoxon rank sum test of normalized gene expression 

# Define samples of interest
samples <- c(
  "V7-HLA-pos",
  "V14-HLA-pos",
  "TADA-HLA-pos",
  "CDA1-HLA-pos",
  "RAPO-HLA-pos",
  "CAS9-HLA-pos"
)

# genes of interest
target_genes <- c("EED", "EZH2", "SUZ12")

# cells of control sample
cell_ctrl <- colnames(be_singlets)[be_singlets$sample == "Ctrl-HLA-pos"]

# Use Seurat's FindMarkers function to perform wilcoxon rank sum test and get p-value and log2FC

# Iterate over sample
findmarkers_res_samples <- lapply(samples, function(selected_sample) {
  # get guides detected in sample
  guide_counts_single_integration_sample <- guide_counts_single_integration %>%
    filter(sample == selected_sample)
  guides <- as.list(unique(guide_counts_single_integration_sample$barcode))
  names(guides) <- guides
  
  cat(selected_sample, length(guides), "\n")
  
  # make list of guide to cells
  guides_cells_sample <- lapply(guides, function(b) {
    guide_counts_single_integration_sample %>% filter(barcode == b) %>% pull(cell_id)
  })
  
  # subset seurat object to sample + ctrl
  be_singlets_sample <- subset(be_singlets,
                               cells = c(as.vector(unlist(
                                 guides_cells_sample
                               )), cell_ctrl),
                               features = target_genes)
  
  # iterate over guides in sample
  findmarkers_res <- lapply(names(guides_cells_sample), function(guide) {
    # perform wilcoxon rank sum test
    de_res <- FindMarkers(
      object = subset(
        be_singlets_sample,
        cells = c(guides_cells_sample[[guide]], cell_ctrl)
      ),
      ident.1 = selected_sample,
      ident.2 = "Ctrl-HLA-pos",
      group.by = "sample",
      features = target_genes,
      logfc.threshold = 0
    )
    de_res$guide <- guide
    de_res$gene <- rownames(de_res)
    return(de_res)
  })
  
  findmarkers_res_df <- bind_rows(findmarkers_res)
  findmarkers_res_df$sample <- selected_sample
  
  return(findmarkers_res_df)
})

# combine into one dataframe. 
findmarkers_res_samples_df <- bind_rows(findmarkers_res_samples)

# the adjusted p-value is not meaningful anymore, so I'll remove it
findmarkers_res_samples_df$p_val_adj <- NULL

write_tsv(
  findmarkers_res_samples_df,
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_target_gene_de_single_guide_vs_ctrl_wilcoxon.tsv"
  )
)
