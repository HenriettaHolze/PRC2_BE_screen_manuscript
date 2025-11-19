#!/usr/bin/env Rscript

# DESCRIPTION
# - Demultiplex variant calling data into samples (i.e. merge with scRNA-seq data)
# - Remove cell barcodes overlapping between captures, due to ambiguity from pooling captures for 
#     nanopore sequencing
# NB: `cellsnp-lite` was run with `--minCOUNT 10`, i.e. only considering positions with at least 
#       10 UMI coverage (to exclude low confidence mapping regions, chimeric reads). 
#     `cellsnp-lite` was also run with `--minMAF 0.0001`, to only consider positions with an 
#       alternative allele (for speed and to reduce output to avoid non-informative output). 

# INPUT
# - cellsnp-lite variant calling results, run on sicelore consensus sequence 
#     (first run mode 2b, then mode 1a)
# - Filtered cell barcodes from Cell Ranger
# - Seurat object with demultiplexed samples

# OUTPUT
# - cellsnp-lite variant calling results, that overlap with scRNA-seq data

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/nanopore/"), recursive = TRUE)

## Load data
be_singlets <-
  qs::qread(
    paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs")
  )
be_singlets$cellid <- sub("P[12]_(.*)-[12]", "\\1", colnames(be_singlets))

# I only need the metadata
be_singlets_metadata <- be_singlets@meta.data
be_singlets <- NULL

# Filtered cell barcodes from Cell Ranger
barcodes_p1 <- read.table(
  paste0(result_dir, "MF03/scRNA/cellranger_results/MF03_P1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
)$V1
barcodes_p2 <- read.table(
  paste0(result_dir, "MF03/scRNA/cellranger_results/MF03_P2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
)$V1

# cellsnp variant data
cellsnp_res_dir <- paste0(result_dir, "MF03/nanopore/")
mat_alt <-
  Matrix::readMM(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.tag.AD.mtx'
    )
  )
mat_ref_alt <-
  Matrix::readMM(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.tag.DP.mtx'
    )
  )
mat_oth <-
  Matrix::readMM(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.tag.OTH.mtx'
    )
  )
vcf <-
  read.table(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.base.vcf.gz'
    ),
    comment.char = "",
    skip = 1,
    header = T
  )
cells <-
  read.table(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.samples.tsv'
    )
  )$V1

dim(mat_alt) == dim(mat_ref_alt)
dim(mat_ref_alt) == dim(mat_oth)

nrow(mat_alt) == nrow(vcf)
ncol(mat_alt) == length(cells)

colnames(mat_alt) <- cells
rownames(mat_alt) <- vcf$POS

colnames(mat_ref_alt) <- cells
rownames(mat_ref_alt) <- vcf$POS

colnames(mat_oth) <- cells
rownames(mat_oth) <- vcf$POS

vcf <- vcf %>%
  separate(INFO, into = c("AD", "DP", "OTH"), sep = ";") %>%
  mutate(AD = as.integer(gsub("AD=", "", AD)),
         DP = as.integer(gsub("DP=", "", DP)),
         OTH = as.integer(gsub("OTH=", "", OTH)))

vcf <- vcf %>%
  mutate(gene = ifelse(X.CHROM == 11, "EED", ifelse(
    X.CHROM == 7, "EZH2", ifelse(X.CHROM == 17, "SUZ12", NA)
  )))

# sum(is.na(vcf$gene))

### Filter cells that are demultiplexed into samples

# Demultiplex cells into samples, remove cell barcodes that are present in both captures
barcodes_overlap <- gsub("-1", "", barcodes_p1[barcodes_p1 %in% barcodes_p2])

mat_alt_samples <-
  mat_alt[, colnames(mat_alt) %in% be_singlets_metadata$cellid &
            !(colnames(mat_alt) %in% barcodes_overlap)]
mat_ref_alt_samples <-
  mat_ref_alt[, colnames(mat_ref_alt) %in% be_singlets_metadata$cellid &
                !(colnames(mat_ref_alt) %in% barcodes_overlap)]
mat_oth_samples <-
  mat_oth[, colnames(mat_oth) %in% be_singlets_metadata$cellid &
            !(colnames(mat_oth) %in% barcodes_overlap)]

# Write to file
Matrix::writeMM(mat_alt_samples, file = paste0(result_dir, "MF03/nanopore/mat_alt_samples_sicelore_cellsnp.mtx"), sep = " ")
Matrix::writeMM(mat_ref_alt_samples, file = paste0(result_dir, "MF03/nanopore/mat_ref_alt_samples_sicelore_cellsnp.mtx"), sep = " ")
Matrix::writeMM(mat_oth_samples, file = paste0(result_dir, "MF03/nanopore/mat_oth_samples_sicelore_cellsnp.mtx"), sep = " ")
write_lines(
  rownames(mat_ref_alt_samples),
  paste0(result_dir, "MF03/nanopore/positions_samples_sicelore_cellsnp.tsv")
)
write_lines(
  colnames(mat_ref_alt_samples),
  paste0(result_dir, "MF03/nanopore/cells_samples_sicelore_cellsnp.tsv")
)
