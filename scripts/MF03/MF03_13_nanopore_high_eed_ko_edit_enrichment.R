#!/usr/bin/env Rscript

# DESCRIPTION
# Identify mutations enriched in cells with high EED KO score, independent from guide detection.
# Do so, by comparing alternative allele frequency to K562 parental and K562 Ctrl HLA+ samples.
# This gets rid of potential sequencing errors, especially non-random sequencing errors.

library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")

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
gc()

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

gc()

# Load VEP results (web server) for edits detected by cellsnp-lite
vep_cellsnp <-
  read_tsv(
    "/dawson_genomics/Projects/PRC2_BE_screen/data/MF03_nanopore/cellsnp_vep/MF03_nanopore_sicelore_cellsnp_1a_vep.txt"
  ) %>%
  mutate(POS = gsub(".*-", "", Location))

####################################################################################################
# Process data
# - Filter cells of base editor samples for PRC2 loss of function by EED KO score > 0.06
# - ONT coverage per sample per position, as number of cells with transcript detected at the position
#     and UMI transcripts
# - Calculate the fraction of cells/transcripts with edit detected per position
# - Calculate log2 FC compared to the fraction in the K562 sample

## Select cells with high EED KO score (>0.06) or K562 or Ctrl HLA+ samples
cells_high_eed_ko <- be_singlets_metadata$cellid[be_singlets_metadata$eed_ko_score_normalized > 0.06]
cells_high_eed_ko <- unique(c(cells_high_eed_ko, be_singlets_metadata$cellid[be_singlets_metadata$sample %in% c("K562", "Ctrl-HLA-pos")]))
length(cells_high_eed_ko)

cells_high_eed_ko <- cells_high_eed_ko[cells_high_eed_ko %in% colnames(mat_alt_samples)]
length(cells_high_eed_ko)

mat_ref_alt_samples_eed_ko <- mat_ref_alt_samples[, cells_high_eed_ko]
mat_alt_samples_eed_ko <- mat_alt_samples[, cells_high_eed_ko]

table(be_singlets_metadata$sample)

be_singlets_metadata <- be_singlets_metadata %>%
  mutate(sample_k562_be = ifelse(
    sample %in% c("Ctrl-HLA-pos", "K562"),
    "K562-Ctrl-HLA-pos",
    sample
  ))

# calculate number of cells with nanopore coverage for each sample at each position
k562_be_samples <- c(
  "K562-Ctrl-HLA-pos",
  "V7-HLA-pos",
  "TADA-HLA-pos",
  "CDA1-HLA-pos",
  "V14-HLA-pos",
  "RAPO-HLA-pos"
)

n_cells_coverage_sample <- lapply(k562_be_samples, function(selected_sample) {
  cells <- be_singlets_metadata$cellid[be_singlets_metadata$sample_k562_be == selected_sample]
  cells <- cells[cells %in% colnames(mat_ref_alt_samples_eed_ko)]
  n_cells_coverage <- Matrix::rowSums(mat_ref_alt_samples_eed_ko[, cells] != 0)
  return(n_cells_coverage)
})
names(n_cells_coverage_sample) <- k562_be_samples
n_cells_coverage_sample_df <- bind_rows(n_cells_coverage_sample, .id = "sample")
n_cells_coverage_sample_df <- n_cells_coverage_sample_df %>%
  pivot_longer(cols = -sample,
               names_to = "POS",
               values_to = "n_cells_sample_coverage")

n_cells_edit_sample <- lapply(k562_be_samples, function(selected_sample) {
  cells <- be_singlets_metadata$cellid[be_singlets_metadata$sample_k562_be == selected_sample]
  cells <- cells[cells %in% colnames(mat_alt_samples_eed_ko)]
  n_cells_coverage <- Matrix::rowSums(mat_alt_samples_eed_ko[, cells] != 0)
  return(n_cells_coverage)
})
names(n_cells_edit_sample) <- k562_be_samples
n_cells_edit_sample_df <- bind_rows(n_cells_edit_sample, .id = "sample")
n_cells_edit_sample_df <- n_cells_edit_sample_df %>%
  pivot_longer(cols = -sample,
               names_to = "POS",
               values_to = "n_cells_sample_edit")

n_umi_coverage_sample <- lapply(k562_be_samples, function(selected_sample) {
  cells <- be_singlets_metadata$cellid[be_singlets_metadata$sample_k562_be == selected_sample]
  cells <- cells[cells %in% colnames(mat_ref_alt_samples_eed_ko)]
  n_umi_coverage <- Matrix::rowSums(mat_ref_alt_samples_eed_ko[, cells])
  return(n_umi_coverage)
})
names(n_umi_coverage_sample) <- k562_be_samples
n_umi_coverage_sample_df <- bind_rows(n_umi_coverage_sample, .id = "sample")
n_umi_coverage_sample_df <- n_umi_coverage_sample_df %>%
  pivot_longer(cols = -sample,
               names_to = "POS",
               values_to = "n_umi_sample_coverage")

n_umi_edit_sample <- lapply(k562_be_samples, function(selected_sample) {
  cells <- be_singlets_metadata$cellid[be_singlets_metadata$sample_k562_be == selected_sample]
  cells <- cells[cells %in% colnames(mat_alt_samples_eed_ko)]
  n_umi_coverage <- Matrix::rowSums(mat_alt_samples_eed_ko[, cells])
  return(n_umi_coverage)
})
names(n_umi_edit_sample) <- k562_be_samples
n_umi_edit_sample_df <- bind_rows(n_umi_edit_sample, .id = "sample")
n_umi_edit_sample_df <- n_umi_edit_sample_df %>%
  pivot_longer(cols = -sample,
               names_to = "POS",
               values_to = "n_umi_sample_edit")


summary_stats_samples <- cbind(
  n_cells_coverage_sample_df,
  n_cells_edit_sample_df,
  n_umi_coverage_sample_df,
  n_umi_edit_sample_df
)

summary_stats_samples <- summary_stats_samples[, c(1, 2, 3, 6, 9, 12)]

# summary_stats_samples %>%
#   filter(n_cells_sample_edit != 0) %>%
#   ggplot(aes(x = n_cells_sample_edit, n_umi_sample_edit)) +
#   geom_point(alpha = 0.2) +
#   theme_bw() +
#   theme(panel.grid = element_blank())

summary_stats_samples <- summary_stats_samples %>%
  mutate(
    fraction_cells_sample = n_cells_sample_edit / n_cells_sample_coverage,
    fraction_umi_sample = n_umi_sample_edit / n_umi_sample_coverage
  ) %>%
  mutate(
    fraction_cells_sample = replace_na(fraction_cells_sample, 0),
    fraction_umi_sample = replace_na(fraction_umi_sample, 0)
  )


# annotate gene to position
summary_stats_samples <- summary_stats_samples %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(gene = ifelse(POS < 40000000, "SUZ12", ifelse(POS < 100000000, "EED", "EZH2")))

# Merge info on what's the reference and alternative allele
summary_stats_samples <- merge(summary_stats_samples, dplyr::select(vcf, c(POS, REF, ALT)), by = "POS")


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


summary_stats_samples <- merge(summary_stats_samples,
                               vep_cellsnp_formatted,
                               by = "POS",
                               all.x = T)

base_change_possible_list <- c("A_G", "C_T", "G_A", "T_C")

summary_stats_samples <- summary_stats_samples %>%
  mutate(editor = gsub("-HLA-pos", "", sample)) %>%
  mutate(base_change = paste0(REF, "_", ALT)) %>%
  mutate(aa_mutation = ifelse(aa_mutation == "-", Consequence_single, aa_mutation)) %>%
  mutate(
    base_change_possible = ifelse(
      editor %in% c("V14", "V7"),
      base_change %in% base_change_possible_list,
      ifelse(
        editor %in% c("RAPO", "CDA1"),
        base_change %in% c("C_T", "G_A"),
        ifelse(editor == "TADA", base_change %in% c("A_G", "T_C"), T)
      )
    )
  )

summary_stats_k562 <- summary_stats_samples %>%
  filter(sample == "K562-Ctrl-HLA-pos") %>%
  dplyr::select(POS, fraction_cells_sample, fraction_umi_sample) %>%
  dplyr::rename(fraction_cells_k562 = fraction_cells_sample,
                fraction_umi_k562 = fraction_umi_sample)

summary_stats_samples <- merge(summary_stats_samples, summary_stats_k562, by = "POS")

summary_stats_samples %>%
  filter(n_cells_sample_edit > 0) %>%
  filter(fraction_cells_k562 > 0) %>%
  filter(sample == "K562-Ctrl-HLA-pos", n_cells_sample_coverage > 100) %>%
  group_by(gene) %>%
  summarize(min(fraction_cells_k562), min(fraction_umi_k562))

summary_stats_samples <- summary_stats_samples %>%
  filter(n_cells_sample_edit > 0) %>%
  mutate(
    log2_fraction_cells_sample_over_fraction_cells_k562 = log2((fraction_cells_sample + 0.001) / (fraction_cells_k562 + 0.001)),
    log2_fraction_umi_sample_over_fraction_umi_k562 = log2((fraction_umi_sample + 0.001) / (fraction_umi_k562 + 0.001))
  )


write_tsv(
  summary_stats_samples,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells.tsv"
  )
)
