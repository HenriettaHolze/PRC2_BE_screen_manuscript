#!/usr/bin/env Rscript

# DESCRIPTION
# Identify heterozygous edits (or frequent sequencing errors) that are high frequency in K562 sample.
# The K562 sample was not sorted for HLA expression. Therefore, it is less likely to be clonally
# restricted due to selection pressure, or to have accuired spontaneous mutations.

# INPUT
# - cellSNP variant calling results, filtered for cells demultiplexed in scRNA-seq data
# - Seurat object with sample annotation

# OUTPUT
# Table with allele frequency of edits in K562 control sample (unsorted)

library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/nanopore/"), recursive = TRUE)

## Load data

cellsnp_res_dir <-
  paste0(result_dir, "MF03/nanopore/")
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
# read in reference and alternative UMI counts for cells that were demultiplexed into samples
mat_alt_samples <-
  Matrix::readMM(paste0(
    result_dir,
    "MF03/nanopore/mat_alt_samples_sicelore_cellsnp.mtx"
  ))
mat_ref_alt_samples <-
  Matrix::readMM(paste0(
    result_dir,
    "MF03/nanopore/mat_ref_alt_samples_sicelore_cellsnp.mtx"
  ))

cells <-
  read_lines(paste0(
    result_dir,
    "MF03/nanopore/cells_samples_sicelore_cellsnp.tsv"
  ))
positions <-
  read_lines(paste0(
    result_dir,
    "MF03/nanopore/positions_samples_sicelore_cellsnp.tsv"
  ))
rownames(mat_alt_samples) <- positions
rownames(mat_ref_alt_samples) <- positions
colnames(mat_alt_samples) <- cells
colnames(mat_ref_alt_samples) <- cells
cells <- NULL
positions <- NULL

# Seurat object with annotation which cell belongs to which sample
be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )
be_singlets$cellid <- sub("P[12]_(.*)-[12]", "\\1", colnames(be_singlets))

# I only need the metadata
be_singlets_metadata <- be_singlets@meta.data
be_singlets <- NULL

# subset to control sample
ctrl_cells <- be_singlets_metadata[be_singlets_metadata$sample == "K562", "cellid"]
ctrl_cells <- ctrl_cells[ctrl_cells %in% colnames(mat_ref_alt_samples)]

mat_ref_alt_samples_ctrl <- mat_ref_alt_samples[, ctrl_cells]
mat_alt_samples_ctrl <- mat_alt_samples[, ctrl_cells]

# get positions that have alternative allele detected in control sample
mask <- Matrix::rowSums(mat_alt_samples_ctrl > 0) > 0
mat_alt_samples_ctrl <- mat_alt_samples_ctrl[mask, ]
mat_ref_alt_samples_ctrl <- mat_ref_alt_samples_ctrl[mask, ]

# make long-format dataframes of alternative and total allele count per cell and position
alt_df <- as.data.frame(Matrix::t(mat_alt_samples_ctrl)) %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(cols = c(-cell_id),
               values_to = "alt_count",
               names_to = "POS")

ref_alt_df <- as.data.frame(Matrix::t(mat_ref_alt_samples_ctrl)) %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(cols = c(-cell_id),
               values_to = "ref_alt_count",
               names_to = "POS")

# combine alternative allele and total allele count in one dataframe
sample_alt_ref_df <-
  cbind(alt_df, ref_alt_df[, "ref_alt_count", drop = F])

sample_alt_ref_df <- sample_alt_ref_df %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(gene = ifelse(POS < 40000000, "SUZ12", ifelse(POS < 100000000, "EED", "EZH2")))

# Calculate the proportion of nanopore reads at each position that have the alternative allele.
# Only calculate this for positions with a minimum coverage, otherwise it will primarily show
# messed up (chimeric) reads.
# EED has the lowest coverage and a median coverage of 400 UMI across positions
ctrl_proportion_edited_transcripts <- sample_alt_ref_df %>%
  group_by(POS) %>%
  mutate(
    total_count_position = sum(ref_alt_count),
    alt_count_position = sum(alt_count),
    frac_alt_count = alt_count_position / total_count_position
  ) %>%
  arrange(desc(frac_alt_count)) %>%
  ungroup() %>%
  dplyr::select(POS,
                gene,
                frac_alt_count,
                total_count_position,
                alt_count_position) %>%
  filter(total_count_position > 100) %>%
  unique() %>%
  mutate(rank = row_number())

write_tsv(
  ctrl_proportion_edited_transcripts,
  paste0(result_dir, "MF03/nanopore/MF03_k562_heterozygosity.tsv")
)
