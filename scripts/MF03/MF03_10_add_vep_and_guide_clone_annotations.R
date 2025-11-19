#!/usr/bin/env Rscript

# DESCRIPTION
# Add VEP and guide clone information to edit x guide table

# INPUT
# - VEP data of cellSNP results
# - Guide UMI count table, filtered and guides aggregated per cell
# - jaccard index table
# - Target genes differential expression results
# - EED KO score t-test results

# OUTPUT
# - Table enriched with VEP annotation of edits, and per-guide clone annotation of target gene
#     differential expression (only guide clones with single guide), EED KO score and co-occuring guides (jaccard)


library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/nanopore/"), recursive = TRUE)

## Load data
sample_alt_ref_df_guides <- read_tsv(paste0(
  result_dir,
  "MF03/nanopore/MF03_sicelore_detected_edits_guides.tsv"
))

# Load VEP results (web server) for edits detected by cellsnp-lite
vep_cellsnp <-
  read_tsv(paste0(
    result_dir,
    "MF03/nanopore/MF03_nanopore_sicelore_cellsnp_1a_vep.txt"
  )) %>%
  mutate(POS = gsub(".*-", "", Location))

# Guide UMI count table to assign cells to guide clones
guide_counts_aggregated <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_aggregated.tsv"
))

# t-test of eed KO score of guide clones vs Ctrl HLA+
ttest_samples_guides_aggregated_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_aggregated_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

# wilcoxon rank sum test of target genes of guide clones vs. Ctrl HLA+
findmarkers_res_samples_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_target_gene_de_single_guide_vs_ctrl_wilcoxon.tsv"
  )
)

jaccard_samples_df_abba <- read_tsv(paste0(result_dir, "MF03/scRNA_seq/MF03_guides_jaccard.tsv"))


## Format tables

# select the relevant transcripts
vep_cellsnp_formatted <- vep_cellsnp %>%
  filter(Feature %in% c("NM_001203247.2", "NM_003797.5", "NM_015355.4")) %>%
  mutate(
    aa_mutation = gsub(".*:p.", "", HGVSp),
    aa_position = as.numeric(Protein_position),
    aa_ref = gsub("^([A-z]*)[0-9].*", "\\1", aa_mutation),
    aa_alt = gsub(".*[0-9]([A-z]*)", "\\1", aa_mutation),
    aa_ref_one_letter = gsub("/.*", "", Amino_acids),
    aa_alt_one_letter = gsub(".*/", "", Amino_acids),
  ) %>%
  mutate(Consequence_single = gsub(",.*", "", Consequence)) %>%
  dplyr::select(
    POS,
    UPLOADED_ALLELE,
    aa_mutation,
    aa_position,
    Consequence_single,
    am_pathogenicity,
    Feature
  )

# add VEP info
sample_alt_ref_df_guides_anno <- merge(sample_alt_ref_df_guides,
                                       vep_cellsnp_formatted,
                                       by = "POS",
                                       all.x = T)

# clean up mess with cell id column names
sample_alt_ref_df_guides_anno <- sample_alt_ref_df_guides_anno %>%
  mutate(cell_barcode = cell_id) %>%
  mutate(cell_id = cell_id_detected_guide)
sample_alt_ref_df_guides_anno$cell_id_detected_guide <- NULL

# add EED KO score t-test
sample_alt_ref_df_guides_anno <- merge(
  sample_alt_ref_df_guides_anno,
  guide_counts_aggregated %>%
    mutate(guide_clone = barcode) %>%
    select(cell_id, guide_clone),
  by = "cell_id",
  all.x = T
) %>%
  mutate(guide_clone = gsub("_", "-", guide_clone))

sample_alt_ref_df_guides_anno <- merge(
  sample_alt_ref_df_guides_anno,
  ttest_samples_guides_aggregated_df,
  by.x = c("guide_clone", "sample"),
  by.y = c("guide", "sample"),
  all.x = T
) %>%
  mutate(n_cells_guide_clone_detected = n_cells)

# add target gene DEG
findmarkers_res_samples_df <- findmarkers_res_samples_df %>%
  mutate(guide = gsub("_", "-", guide)) %>%
  select(p_val, avg_log2FC, guide, sample, gene) %>%
  pivot_wider(
    id_cols = c(guide, sample),
    names_from = gene,
    values_from = c(p_val, avg_log2FC)
  )

sample_alt_ref_df_guides_anno <- merge(
  sample_alt_ref_df_guides_anno,
  findmarkers_res_samples_df,
  by.x = c("guide_clone", "sample"),
  by.y = c("guide", "sample"),
  all.x = T
)


# add information on which guide are detected together with a high jaccard index
jaccard_samples_df_abba_high <- jaccard_samples_df_abba %>%
  filter(jaccard_index > 0.5) %>%
  # sort alphabetically
  arrange(guide_b) %>%
  group_by(sample, guide_a) %>%
  summarize(guides_jaccard = paste(guide_b, collapse = ";")) %>%
  mutate(
    guides_jaccard = gsub("_", "-", guides_jaccard),
    guide_a =  gsub("_", "-", guide_a)
  )

sample_alt_ref_df_guides_anno <- merge(
  sample_alt_ref_df_guides_anno,
  jaccard_samples_df_abba_high,
  by.x = c("sample", "guide"),
  by.y = c("sample", "guide_a"),
  all.x = T
)

write_tsv(
  sample_alt_ref_df_guides_anno,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation.tsv"
  )
)
