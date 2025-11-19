#!/usr/bin/env Rscript

# DESCRIPTION
# Calculate the alternative allele frequency for all edits for cells in guide clones. 
# Calculate both the number of cells with ONT transcript detected at a position and the total number of transcripts.
# Divide by the counts for just the alternative allele.
# So you get a "cell frequency" and "allele frequency"

# INPUT
# - cellSNP variant calling results, filtered for cells demultiplexed in scRNA-seq data
# - Filtered guide counts, guides aggregated by cell
# - annotation which cell belongs to which sample from Seurat object

library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")

# cellsnp data that was filtered for cells that can be merged with scRNA-seq data
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
gc()

# Guide UMI count table to assign cells to guide clones
guide_counts_aggregated <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_aggregated.tsv"
))

guide_counts_aggregated <- guide_counts_aggregated %>% mutate(barcode = gsub("_", "-", barcode))

sample_alt_ref_df_guides_anno <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation.tsv"
  )
)

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

####################################################################################################



# get guides annotated to cells for UMI threshold
guide_counts_filtered <- guide_counts %>%
  filter(bc.umi.count >= 3) %>%
  mutate(guide = gsub("_", "-", barcode)) %>%
  group_by(sample, guide) %>%
  mutate(n_cells_guide_detected = n()) %>%
  filter(n_cells_guide_detected >= 10) %>%
  dplyr::select(cell_id, cellid, guide, sample)


be_samples <- c("CDA1-HLA-pos",
                "TADA-HLA-pos",
                "RAPO-HLA-pos",
                "V7-HLA-pos",
                "V14-HLA-pos")

guide_counts_filtered_be <- guide_counts_aggregated %>%
  filter(sample %in% be_samples) %>%
  mutate(guide = barcode) %>%
  dplyr::select(cell_id, cellid, guide, sample)


genes <- c("EED", "EZH2", "SUZ12")

guide_be_combinations <- guide_counts_filtered_be %>%
  select(guide, sample) %>%
  unique()

guide_counts_filtered_be_ont <- guide_counts_filtered_be %>%
  filter(cellid %in% colnames(mat_alt_samples))

nrow(guide_be_combinations)

be_singlets_metadata_ont <- be_singlets_metadata %>%
  filter(cellid %in% colnames(mat_ref_alt_samples))

be_singlets_metadata_ont <- be_singlets_metadata_ont[rownames(be_singlets_metadata_ont) %in% guide_counts_filtered_be$cell_id, ]

mat_alt_samples_guides_be <- mat_alt_samples[, be_singlets_metadata_ont$cellid]
mat_ref_alt_samples_guides_be <- mat_ref_alt_samples[, be_singlets_metadata_ont$cellid]

# 203 guide clone x base editor combinations

n_cells_coverage_be_guide <- lapply(seq(nrow(guide_be_combinations)), function(i) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  selected_sample <- guide_be_combinations[i, ]$sample
  selected_guide <- guide_be_combinations[i, ]$guide
  cells <- guide_counts_filtered_be_ont %>% filter(sample == selected_sample, guide == selected_guide) %>%
    pull(cellid)
  n_cells_coverage <- Matrix::rowSums(mat_ref_alt_samples_guides_be[, cells] > 0)
  return(n_cells_coverage)
})
names(n_cells_coverage_be_guide) <- paste0(guide_be_combinations$sample,
                                           "_",
                                           guide_be_combinations$guide)
n_cells_coverage_be_guide_df <- bind_rows(n_cells_coverage_be_guide, .id = "sample_guide_clone")
n_cells_coverage_be_guide_df <- n_cells_coverage_be_guide_df %>%
  pivot_longer(
    cols = -sample_guide_clone,
    names_to = "POS",
    values_to = "n_cells_sample_guide_clone_coverage"
  ) %>%
  separate_wider_delim(
    cols = sample_guide_clone,
    delim = "_",
    names = c("sample", "guide_clone")
  )

n_cells_edit_be_guide <- lapply(seq(nrow(guide_be_combinations)), function(i) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  selected_sample <- guide_be_combinations[i, ]$sample
  selected_guide <- guide_be_combinations[i, ]$guide
  cells <- guide_counts_filtered_be_ont %>% filter(sample == selected_sample, guide == selected_guide) %>%
    pull(cellid)
  n_cells_edit <- Matrix::rowSums(mat_alt_samples_guides_be[, cells] > 0)
  return(n_cells_edit)
})
names(n_cells_edit_be_guide) <- paste0(guide_be_combinations$sample,
                                       "_",
                                       guide_be_combinations$guide)
n_cells_edit_be_guide_df <- bind_rows(n_cells_edit_be_guide, .id = "sample_guide_clone")
n_cells_edit_be_guide_df <- n_cells_edit_be_guide_df %>%
  pivot_longer(
    cols = -sample_guide_clone,
    names_to = "POS",
    values_to = "n_cells_sample_guide_clone_edit"
  ) %>%
  separate_wider_delim(
    cols = sample_guide_clone,
    delim = "_",
    names = c("sample", "guide_clone")
  )


n_umi_coverage_be_guide <- lapply(seq(nrow(guide_be_combinations)), function(i) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  selected_sample <- guide_be_combinations[i, ]$sample
  selected_guide <- guide_be_combinations[i, ]$guide
  cells <- guide_counts_filtered_be_ont %>% filter(sample == selected_sample, guide == selected_guide) %>%
    pull(cellid)
  n_umi_coverage <- Matrix::rowSums(mat_ref_alt_samples_guides_be[, cells])
  return(n_umi_coverage)
})
names(n_umi_coverage_be_guide) <- paste0(guide_be_combinations$sample,
                                         "_",
                                         guide_be_combinations$guide)
n_umi_coverage_be_guide_df <- bind_rows(n_umi_coverage_be_guide, .id = "sample_guide_clone")
n_umi_coverage_be_guide_df <- n_umi_coverage_be_guide_df %>%
  pivot_longer(
    cols = -sample_guide_clone,
    names_to = "POS",
    values_to = "n_umi_sample_guide_clone_coverage"
  ) %>%
  separate_wider_delim(
    cols = sample_guide_clone,
    delim = "_",
    names = c("sample", "guide_clone")
  )

n_umi_edit_be_guide <- lapply(seq(nrow(guide_be_combinations)), function(i) {
  if (i %% 100 == 0) {
    cat(i, "\n")
  }
  selected_sample <- guide_be_combinations[i, ]$sample
  selected_guide <- guide_be_combinations[i, ]$guide
  cells <- guide_counts_filtered_be_ont %>% filter(sample == selected_sample, guide == selected_guide) %>%
    pull(cellid)
  n_umi_edit <- Matrix::rowSums(mat_alt_samples_guides_be[, cells])
  return(n_umi_edit)
})
names(n_umi_edit_be_guide) <- paste0(guide_be_combinations$sample,
                                     "_",
                                     guide_be_combinations$guide)
n_umi_edit_be_guide_df <- bind_rows(n_umi_edit_be_guide, .id = "sample_guide_clone")
n_umi_edit_be_guide_df <- n_umi_edit_be_guide_df %>%
  pivot_longer(
    cols = -sample_guide_clone,
    names_to = "POS",
    values_to = "n_umi_sample_guide_clone_edit"
  ) %>%
  separate_wider_delim(
    cols = sample_guide_clone,
    delim = "_",
    names = c("sample", "guide_clone")
  )


summary_stats_samples_guides <- cbind(
  n_cells_coverage_be_guide_df,
  n_cells_edit_be_guide_df,
  n_umi_coverage_be_guide_df,
  n_umi_edit_be_guide_df
)

summary_stats_samples_guides_df <- summary_stats_samples_guides[, c(1, 2, 3, 4, 8, 12, 16)]

write_tsv(
  summary_stats_samples_guides_df,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_be_guide_clones_10_cells.tsv"
  )
)

summary_stats_samples_guides_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_be_guide_clones_10_cells.tsv"
  )
)

head(summary_stats_samples_guides_df)

summary_stats_samples_guides_df %>%
  filter(n_cells_sample_guide_clone_edit > 0) %>%
  nrow()

# filter positions to be detected in at least 3 cells, else it does not make sense to calculate allele frequency
summary_stats_samples_guides_df <- summary_stats_samples_guides_df %>%
  filter(n_cells_sample_guide_clone_edit >= 3) %>%
  mutate(
    cells_edit_vs_coverage = n_cells_sample_guide_clone_edit / n_cells_sample_guide_clone_coverage,
    umi_edit_vs_coverage = n_umi_sample_guide_clone_edit / n_umi_sample_guide_clone_coverage
  )

nrow(summary_stats_samples_guides_df)

sample_alt_ref_df_guides_anno <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation.tsv"
  )
)

sample_alt_ref_df_guides_anno_af <- merge(
  sample_alt_ref_df_guides_anno,
  summary_stats_samples_guides_df,
  by = c("sample", "guide_clone", "POS"),
  all.x = T
)

write_tsv(
  sample_alt_ref_df_guides_anno_af,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af.tsv"
  )
)
