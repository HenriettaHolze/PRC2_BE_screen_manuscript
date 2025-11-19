#!/usr/bin/env Rscript

library(tidyverse)
library(seqinr)
result_dir <- Sys.getenv("RESULT_DIR")


sample_alt_ref_df_guides_anno <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af.tsv"
  )
)

# only consider base editor samples
sample_alt_ref_df_guides_anno <- sample_alt_ref_df_guides_anno %>%
  filter(!sample %in% c("CAS9-HLA-pos", "Ctrl-HLA-pos"))

##############################
# Make filtering columns

colnames_guide_clone_edit <- c(
  'sample',
  'guide',
  'guide_clone',
  'POS',
  'gene',
  'REF',
  'ALT',
  'Strand.of.sgRNA',
  'sgRNA.Cut.Position..1.based.',
  'k562_alt_freq',
  'position_on_guide',
  'detected_edit',
  'editor',
  'position_on_guide_label',
  'correct_edits',
  'guide_target',
  'UPLOADED_ALLELE',
  'aa_mutation',
  'aa_mutation_1letter',
  'aa_position',
  'Consequence_single',
  'am_pathogenicity',
  'Feature',
  'diff_eed_ko_score_normalized',
  'p.value_eed_ko_score_normalized',
  'eed_mean',
  'ezh2_mean',
  'suz12_mean',
  'p.val.adj_eed_ko_score_normalized',
  'n_cells_guide_clone_detected',
  'p_val_EZH2',
  'p_val_SUZ12',
  'p_val_EED',
  'avg_log2FC_EZH2',
  'avg_log2FC_SUZ12',
  'avg_log2FC_EED',
  'guides_jaccard',
  "cells_edit_vs_coverage",
  "n_cells_sample_guide_clone_edit",
  "n_cells_sample_guide_clone_coverage",
  "umi_edit_vs_coverage",
  "n_umi_sample_guide_clone_edit",
  "n_umi_sample_guide_clone_coverage"
)

# make better readable aa mutation labels
sample_alt_ref_df_guides_anno <- sample_alt_ref_df_guides_anno %>%
  mutate(aa_mutation_1letter = gsub("%3D", "", aa_mutation)) %>%
  mutate(aa_mutation_1letter = gsub("ext.*", "", aa_mutation_1letter)) %>%
  mutate(aa_mutation_1letter = ifelse(is.na(aa_mutation_1letter), "-", aa_mutation_1letter))

# convert 3 letter code into 1 letter code
aa_mutation_1letter_list <- unique(sample_alt_ref_df_guides_anno$aa_mutation_1letter)
names(aa_mutation_1letter_list) <- aa_mutation_1letter_list

aa_mutation_1letter_list <- lapply(aa_mutation_1letter_list, function(aa_mutation_1letter) {
  aa1 = a(gsub("[0-9]+.*", "", aa_mutation_1letter))
  aa2 = a(gsub(".*[0-9]+", "", aa_mutation_1letter))
  aa_position = gsub("[^0-9]*([0-9]+)[^0-9]*", "\\1", aa_mutation_1letter)
  aa2 = ifelse(grepl("[0-9]Ter", aa_mutation_1letter), "*", aa2)
  aa1 = ifelse(grepl("Ter[0-9]", aa_mutation_1letter), "*", aa1)
  aa2 = ifelse(is.na(aa2), "", aa2)
  aa_mutation_1letter = paste0(aa1, aa_position, aa2)
  return(aa_mutation_1letter)
})

sample_alt_ref_df_guides_anno[90,] %>% as.data.frame()

aa_mutation_1letter_column <- as.vector(unlist(aa_mutation_1letter_list[sample_alt_ref_df_guides_anno$aa_mutation_1letter]))
# length(aa_mutation_1letter_column) == nrow(sample_alt_ref_df_guides_anno)
sample_alt_ref_df_guides_anno$aa_mutation_1letter <- aa_mutation_1letter_column
sample_alt_ref_df_guides_anno <- sample_alt_ref_df_guides_anno %>% 
  mutate(aa_mutation_1letter = ifelse(is.na(aa_position), NA, aa_mutation_1letter))

sample_alt_ref_df_guides_anno_filter_columns <- sample_alt_ref_df_guides_anno %>%
  # mark heterozygous edits. NA means edit is not detected in K562 at all.
  mutate(edit_not_k562_heterozygous = (k562_alt_freq < 0.1 |
                                         is.na(k562_alt_freq)))

# I need to remove heterozygous edits before identifying the most frequent edit per guide clone.

# # Removed edits
# sample_alt_ref_df_guides_anno_filter_columns %>%
#   filter(!edit_not_k562_heterozygous) %>%
#   select(gene, POS, aa_mutation) %>%
#   unique() %>%
#   as.data.frame()
# # gene      POS aa_mutation
# # 1   EED 86257581   Leu207%3D
# # 2 SUZ12 31937006           -

sample_alt_ref_df_guides_anno_filter_columns <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(edit_not_k562_heterozygous)

# make filtering columns for all kinds of different conditions that indicate functional mutations
sample_alt_ref_df_guides_anno_filter_columns <- sample_alt_ref_df_guides_anno_filter_columns %>%
  # edits must be detected in at least 3 cells per guide clone
  group_by(across(all_of(colnames_guide_clone_edit)), edit_not_k562_heterozygous) %>%
  summarize(n_cells_edit_guide_clone_combination = n()) %>%
  mutate(edit_3_cells_guide_clone = n_cells_edit_guide_clone_combination >= 3) %>%
  ungroup() %>%
  # check which edit is most frequent per guide clone (among all edits)
  group_by(guide_clone, guide, sample) %>%
  mutate(
    max_cells = max(n_cells_edit_guide_clone_combination, na.rm = T),
    edit_is_most_frequent = max_cells == n_cells_edit_guide_clone_combination,
    # allele frequency only calculated for positions with 3 cells with edit, otherwise not reliable.
    # Then set to NA.
    max_allele_frequency = max(umi_edit_vs_coverage, na.rm = T),
    edit_highest_allele_frequency = max_allele_frequency == umi_edit_vs_coverage,
    max_cell_frequency = max(cells_edit_vs_coverage, na.rm = T),
    edit_highest_cell_frequency = max_cell_frequency == cells_edit_vs_coverage,
  ) %>%
  ungroup() %>%
  # check for edits possibly induced by the guide and editor
  mutate(
    edit_on_guide_window = (
      !is.na(position_on_guide) &
        position_on_guide < 20 &
        position_on_guide > -5
    ),
    edit_matches_guide_and_be = (edit_on_guide_window &
                                   correct_edits == TRUE)
  ) %>%
  # guide clones with only a single guide, and remove guides with a jaccard index > 0.5 with another guide
  mutate(
    guide_clone_single_guide = !grepl(";", guide_clone) &
      !is.na(guide_clone) &
      is.na(guides_jaccard)
  ) %>%
  # check for higher EED KO score than in ctrl
  mutate(
    guide_high_eed_ko = (
      diff_eed_ko_score_normalized > 0.06 &
        p.value_eed_ko_score_normalized < 0.05
    )
  ) %>%
  # remove SUZ12 guides where other genes are downregulated (likely other guides present)
  mutate(guide_suz12_no_off_target_nmd = !(guide_target == "SUZ12" &
                                             ((p_val_EZH2 < 0.05 &
                                                 avg_log2FC_EZH2 < 0) |
                                                (p_val_EED < 0.05 &
                                                   avg_log2FC_EED < 0)
                                             ))) %>%
  # remove EED guides where other genes are downregulated (likely other guides present)
  mutate(guide_eed_no_off_target_nmd = !(guide_target == "EED" &
                                           ((p_val_EZH2 < 0.05 &
                                               avg_log2FC_EZH2 < 0) |
                                              (p_val_SUZ12 < 0.05 &
                                                 avg_log2FC_SUZ12 < 0)
                                           ))) %>%
  # remove EZH2 guides where other genes are downregulated (likely other guides present)
  mutate(guide_ezh2_no_off_target_nmd = !(guide_target == "EZH2" &
                                            ((p_val_EED < 0.05 &
                                                avg_log2FC_EED < 0) |
                                               (p_val_SUZ12 < 0.05 &
                                                  avg_log2FC_SUZ12 < 0)
                                            )))


write_tsv(
  sample_alt_ref_df_guides_anno_filter_columns,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters.tsv"
  )
)

### Save different versions of the file to open in excel

# same as before, but don't remove guide clones that have other target downregulated
sample_alt_ref_df_guides_anno_export <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(
    edit_not_k562_heterozygous &
      edit_3_cells_guide_clone &
      (edit_matches_guide_and_be | edit_is_most_frequent)
  )

write_tsv(
  sample_alt_ref_df_guides_anno_export,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters_match_guide_or_most_frequent.tsv"
  )
)

# filter mutations that we are confident that they are disrupting PRC2
sample_alt_ref_df_guides_anno_functional_no_nmd_filter <- sample_alt_ref_df_guides_anno_filter_columns %>%
  # true across all filtering columns
  filter(if_all(all_of(
    c(
      "edit_not_k562_heterozygous",
      "edit_matches_guide_and_be",
      "edit_3_cells_guide_clone",
      "edit_is_most_frequent",
      "guide_clone_single_guide",
      "guide_high_eed_ko"
    )
  ), ~ . == TRUE))

write_tsv(
  sample_alt_ref_df_guides_anno_functional_no_nmd_filter,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters_functional_no_nmd_filter.tsv"
  )
)

# filter the bare minimum, i.e. edits that I'm guide confident were induced by guides and the base editors
sample_alt_ref_df_guides_anno_induced <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(
    edit_not_k562_heterozygous,
    edit_matches_guide_and_be,
    edit_3_cells_guide_clone,
    !is.na(guide_clone)
  )

write_tsv(
  sample_alt_ref_df_guides_anno_induced,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters_induced.tsv"
  )
)


####################################################################################################
# Add additional information to enriched mutations: COSMIC and ClinVar mutations at aa positions

# subset to missense variants, because COSMIC is also subsetted to missense variants.
# and subset to guide clones with high eed ko score
sample_alt_ref_df_guides_anno_induced_missense <- sample_alt_ref_df_guides_anno_induced %>% 
  filter(Consequence_single == "missense_variant") %>%
  filter(guide_high_eed_ko)

sample_alt_ref_df_guides_anno_induced_missense %>% select(gene, aa_mutation) %>% unique() %>% nrow()
# 107 aa mutations

### COSMIC
# Cross-reference with COSMIC data: how many mutations are in COSMIC, how many are annotated as heterozygous vs. heterozygous.
# see preprocessing /dawson_genomics/Projects/PRC2_BE_screen/manuscript/scripts/cosmic/cosmic_zygosity_v100.R

cosmic_db <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense.tsv"
)

# summarize COSMIC data, how many individuals are recorded with a specific aa_mutation
# and all aa_mutations recorded at a specific aa_position
cosmic_db_summary_aa_mutations <- cosmic_db %>%
  group_by(INDIVIDUAL_ID, aa_mutation, GENE_SYMBOL, aa_position) %>%
  select(INDIVIDUAL_ID,
         aa_mutation,
         GENE_SYMBOL,
         aa_position,
         MUTATION_ZYGOSITY) %>%
  unique() %>%
  group_by(aa_mutation, GENE_SYMBOL, aa_position) %>%
  mutate(cosmic_n_individuals_aa_mutation = n()) %>%
  mutate(cosmic_heterozygous = grepl("het", paste0(unique(MUTATION_ZYGOSITY), collapse = ";"))) %>%
  group_by(GENE_SYMBOL, aa_position) %>%
  arrange(aa_mutation) %>%
  mutate(cosmic_aa_mutations_at_aa_position = paste0(unique(aa_mutation), collapse = ";")) %>%
  mutate(aa_mutation_1letter = aa_mutation) %>%
  select(
    GENE_SYMBOL,
    aa_mutation_1letter,
    aa_position,
    cosmic_n_individuals_aa_mutation,
    cosmic_aa_mutations_at_aa_position,
    cosmic_heterozygous
  ) %>%
  unique() %>%
  ungroup()

cosmic_db_haema <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense_haematopoietic_neoplasm.tsv"
)

# summarize COSMIC data, how many individuals are recorded with a specific aa_mutation
# and all aa_mutations recorded at a specific aa_position
cosmic_db_haema_summary_aa_mutations <- cosmic_db_haema %>%
  group_by(INDIVIDUAL_ID, aa_mutation, GENE_SYMBOL, aa_position) %>%
  select(INDIVIDUAL_ID,
         aa_mutation,
         GENE_SYMBOL,
         aa_position,
         MUTATION_ZYGOSITY) %>%
  unique() %>%
  group_by(aa_mutation, GENE_SYMBOL, aa_position) %>%
  mutate(cosmic_haematopoietic_n_individuals_aa_mutation = n()) %>%
  mutate(cosmic_haematopoietic_heterozygous = grepl("het", paste0(unique(MUTATION_ZYGOSITY), collapse = ";"))) %>%
  group_by(GENE_SYMBOL, aa_position) %>%
  arrange(aa_mutation) %>%
  mutate(cosmic_haematopoietic_aa_mutations_at_aa_position = paste0(unique(aa_mutation), collapse = ";")) %>%
  mutate(aa_mutation_1letter = aa_mutation) %>%
  select(
    GENE_SYMBOL,
    aa_mutation_1letter,
    aa_position,
    cosmic_haematopoietic_n_individuals_aa_mutation,
    cosmic_haematopoietic_aa_mutations_at_aa_position,
    cosmic_haematopoietic_heterozygous
  ) %>%
  unique() %>%
  ungroup()

sample_alt_ref_df_guides_anno_induced_missense_cosmic <- merge(
  sample_alt_ref_df_guides_anno_induced_missense,
  select(
    cosmic_db_summary_aa_mutations,
    GENE_SYMBOL,
    aa_mutation_1letter,
    cosmic_n_individuals_aa_mutation,
    cosmic_heterozygous
  ),
  by.x = c("gene", "aa_mutation_1letter"),
  by.y = c("GENE_SYMBOL", "aa_mutation_1letter"),
  all.x = T
)
sample_alt_ref_df_guides_anno_induced_missense_cosmic <- merge(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic,
  select(
    cosmic_db_summary_aa_mutations,
    GENE_SYMBOL,
    aa_position,
    cosmic_aa_mutations_at_aa_position
  ) %>% unique(),
  by.x = c("gene", "aa_position"),
  by.y = c("GENE_SYMBOL", "aa_position"),
  all.x = T
)

sample_alt_ref_df_guides_anno_induced_missense_cosmic <- merge(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic,
  select(
    cosmic_db_haema_summary_aa_mutations,
    GENE_SYMBOL,
    aa_mutation_1letter,
    cosmic_haematopoietic_n_individuals_aa_mutation,
    cosmic_haematopoietic_heterozygous
  ),
  by.x = c("gene", "aa_mutation_1letter"),
  by.y = c("GENE_SYMBOL", "aa_mutation_1letter"),
  all.x = T
)
sample_alt_ref_df_guides_anno_induced_missense_cosmic <- merge(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic,
  select(
    cosmic_db_haema_summary_aa_mutations,
    GENE_SYMBOL,
    aa_position,
    cosmic_haematopoietic_aa_mutations_at_aa_position
  ) %>% unique(),
  by.x = c("gene", "aa_position"),
  by.y = c("GENE_SYMBOL", "aa_position"),
  all.x = T
)


#### ClinVar
clinvar <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/clinvar/variant_summary_hg38_cosmic_EED_EZH2_SUZ12_likely_path_path_lof_diseases.tsv"
)

clinvar <- clinvar %>%
  group_by(GeneSymbol, aa_position) %>%
  arrange(aa_mutation_1letter) %>%
  mutate(clinvar_number_submitters = NumberSubmitters) %>%
  mutate(clinvar_aa_mutations_at_aa_position = paste0(unique(aa_mutation_1letter), collapse = ";")) %>%
  ungroup() %>%
  select(
    GeneSymbol,
    aa_position,
    aa_mutation_1letter,
    clinvar_aa_mutations_at_aa_position,
    clinvar_number_submitters
  ) %>%
  unique()

sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar <- merge(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic,
  select(
    clinvar,
    GeneSymbol,
    aa_mutation_1letter,
    clinvar_number_submitters
  ),
  by.x = c("gene", "aa_mutation_1letter"),
  by.y = c("GeneSymbol", "aa_mutation_1letter"),
  all.x = T
)
sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar <- merge(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar,
  select(
    clinvar,
    GeneSymbol,
    aa_position,
    clinvar_aa_mutations_at_aa_position
  ) %>% unique(),
  by.x = c("gene", "aa_position"),
  by.y = c("GeneSymbol", "aa_position"),
  all.x = T
)


# check the enriched mutations that have any evidence in either COSMIC or ClinVar
sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar %>%
  filter(
    !is.na(clinvar_aa_mutations_at_aa_position) |
      !is.na(cosmic_aa_mutations_at_aa_position)
  )

sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar %>%
  filter(
    !is.na(clinvar_aa_mutations_at_aa_position) |
      !is.na(cosmic_aa_mutations_at_aa_position)
  ) %>% select(gene, aa_mutation) %>% unique() %>% nrow()
# 48 aa mutations

write_tsv(
  sample_alt_ref_df_guides_anno_induced_missense_cosmic_clinvar,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters_induced_missense_eed_high_cosmic_clinvar.tsv"
  )
)
