#!/usr/bin/env Rscript

# DESCRIPTION
# Identify the editing window of base editors
# - for that, first filter the guides for >= 3 UMI detected
# - identify the edits detected around guide binding sites
# NB: this script takes >32GB memory

# INPUT
# - cellSNP variant calling results, filtered for cells demultiplexed in scRNA-seq data
# - Guide UMI count table
# - Guide metadata (cut site, strand)

# OUTPUT
# - Table with edits and detected guides per cell, and position of edit relative to detected guide

library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/nanopore/"), recursive = TRUE)

####################################################################################################
## Load data

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

# Guide UMI count table
guide_counts <- read_tsv(paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts.tsv"))

# Load guide metadata (cut position and strand).
guide_metadata <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12.tsv"
    )
  )
guide_metadata$sgRNA_ID <- gsub("_", "-", guide_metadata$sgRNA_ID)

# Load frequency of edits in K562 sample (heterozygosity or sequencing errors)
ctrl_alt_freq <- read_tsv(paste0(result_dir, "MF03/nanopore/MF03_k562_heterozygosity.tsv"))
ctrl_alt_freq <- ctrl_alt_freq %>%
  dplyr::rename(k562_alt_freq = frac_alt_count) %>%
  dplyr::select(POS, k562_alt_freq)


####################################################################################################
## Filter variant data to cells with guide detected and positions with alternative alleles

# get guides annotated to cells for UMI threshold
guide_counts_filtered <- guide_counts %>%
  filter(bc.umi.count >= 3) %>%
  mutate(guide = gsub("_", "-", barcode)) %>%
  dplyr::select(cell_id, cellid, guide, sample)

# subset matrices to cells, to be more efficient
cells_with_guide <- guide_counts_filtered$cellid[guide_counts_filtered$cellid %in% colnames(mat_alt_samples)]

mat_ref_alt_samples_guides <- mat_ref_alt_samples[, cells_with_guide]
mat_alt_samples_guides <- mat_alt_samples[, cells_with_guide]

# subset count matrices and vcf file to positions with with alternative allele detected
mask <- Matrix::rowSums(mat_alt_samples_guides) > 0
mat_ref_alt_samples_guides <- mat_ref_alt_samples_guides[mask, ]
mat_alt_samples_guides <- mat_alt_samples_guides[mask, ]
vcf_filtered <- vcf %>%
  filter(POS %in% rownames(mat_alt_samples_guides))

cat(dim(mat_ref_alt_samples), "\n")
# 25769 24753 
cat(dim(mat_ref_alt_samples_guides), "\n")
# 18326 24103
# 24103/24753 cells have a guide detected
# 18326/25769 positions where at least 1 cell has a guide and an alternative allele

# clean up
mask <- NULL
vcf <- NULL
mat_alt_samples <- NULL
mat_ref_alt_samples <- NULL
gc()

# make long-format dataframes of alternative and total allele count per cell and position
alt_df <- as.data.frame(as.matrix(Matrix::t(mat_alt_samples_guides))) %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(cols = c(-cell_id),
               values_to = "alt_count",
               names_to = "POS")

ref_alt_df <- as.data.frame(as.matrix(Matrix::t(mat_ref_alt_samples_guides))) %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(cols = c(-cell_id),
               values_to = "ref_alt_count",
               names_to = "POS")

# combine alternative allele and total allele count in one dataframe
sample_alt_ref_df <-
  cbind(alt_df, ref_alt_df[, "ref_alt_count", drop = F])

# Remove all cell x position rows without alternative allele to reduce size
sample_alt_ref_df <- sample_alt_ref_df %>%
  filter(alt_count != 0)

# annotate gene to position
sample_alt_ref_df <- sample_alt_ref_df %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(gene = ifelse(POS < 40000000, "SUZ12", ifelse(POS < 100000000, "EED", "EZH2")))

# Merge info on what's the reference and alternative allele
sample_alt_ref_df <- merge(sample_alt_ref_df, dplyr::select(vcf_filtered, c(POS, REF, ALT)), by = "POS")

# clean up environment
alt_df <- NULL
ref_alt_df <- NULL
mat_alt_samples <- NULL
mat_ref_alt_samples <- NULL
mat_alt_samples_guides <- NULL
mat_ref_alt_samples_guides <- NULL
vcf_filtered <- NULL
gc()


####################################################################################################
## Annotate the relative position of the edit to the guide binding site that's detected in the cell

# Add guide detection data into variant data
# First add metadata to guide count table
guide_counts_filtered <- merge(
  guide_counts_filtered,
  guide_metadata %>% dplyr::select(sgRNA_ID, Strand.of.sgRNA, sgRNA.Cut.Position..1.based.),
  by.x = "guide",
  by.y = "sgRNA_ID"
)

# Then merge variant data and guide detection data based on cell barcode
sample_alt_ref_df_guides <- merge(
  sample_alt_ref_df,
  guide_counts_filtered,
  by.x = "cell_id",
  by.y = "cellid",
  suffixes = c("", "_detected_guide"),
)

# Merge K562 allele frequency information
sample_alt_ref_df_guides <- merge(sample_alt_ref_df_guides,
                                  ctrl_alt_freq,
                                  by = "POS",
                                  all.x = T)

# Cut site is annotated as guide position 17 for guides on (-) strans and position 18 for guides
# on (+) strand (by [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public).
# I checked with UCSC for some guides, that position - cut position + 18 for guide on (+)
# and cut position - position + 17 for guide on (-)
# gives the position of the edit on the guide sequence (1-based).
sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>% mutate(
  position_on_guide = ifelse(
    Strand.of.sgRNA == "+",
    POS - sgRNA.Cut.Position..1.based. + 18,
    sgRNA.Cut.Position..1.based. - POS + 17
  )
)

# annotate what the edit is, given the strand of the guide
getCompl = list(
  "A" = "T",
  "T" = "A",
  "C" = "G",
  "G" = "C"
)

sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>%
  mutate(detected_edit = ifelse(
    Strand.of.sgRNA == "+",
    paste0(REF, " -> ", ALT),
    paste0(getCompl[REF], " -> ", getCompl[ALT])
  ))

sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>%
  mutate(editor = gsub("-HLA-pos", "", sample))

# Make pretty human-readable labels
sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>%
  # there is no position 0, instead change to -1
  mutate(
    position_on_guide_label = ifelse(position_on_guide <= 0, position_on_guide - 1, position_on_guide)
  ) %>%
  # last 3 positions are PAM sequence
  mutate(position_on_guide_label = ifelse(
    position_on_guide_label == 21,
    "P",
    ifelse(
      position_on_guide_label == 22,
      "A",
      ifelse(position_on_guide_label == 23, "M", position_on_guide_label)
    )
  ))

# annotate whether the edit matches the capacity of the base editor or not
sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>%
  mutate(correct_edits = ifelse(
    editor %in% c("CDA1", "RAPO"),
    detected_edit == "C -> T",
    ifelse(
      editor == "TADA",
      detected_edit == "A -> G",
      ifelse(
        editor %in% c("V7", "V14"),
        detected_edit %in% c("A -> G", "C -> T"),
        NA
      )
    )
  ))

# only annotate realtive position of edit to guide if they are actually on the same gene.
sample_alt_ref_df_guides <- sample_alt_ref_df_guides %>% 
  mutate(guide_target = gsub("-.*", "", guide)) %>%
  mutate(position_on_guide = ifelse(guide_target == gene, position_on_guide, NA),
         position_on_guide_label = ifelse(guide_target == gene, position_on_guide_label, NA))


# Save the table with information of edits in cells around the binding site of the detected guide
write_tsv(
  sample_alt_ref_df_guides,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides.tsv"
  )
)
