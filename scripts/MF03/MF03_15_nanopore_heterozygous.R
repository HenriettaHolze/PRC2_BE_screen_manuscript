#!/usr/bin/env Rscript

# DESCRIPTION
# - identify mutations that are potentially heterozygous, given that both alternative and reference
#     allele are detected in the same cell.
# - For this analysis, filter mutations that do not co-occur with another mutation

library(tidyverse)
library(ggpubr)
library(seqinr)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

summary_stats_samples_be_mean_plot <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells_z_score.tsv"
  )
)

enriched_edits_cooccurrence <- read_tsv(paste0(
  result_dir,
  "MF03/nanopore/MF03_enriched_edit_cooccurrence.tsv"
))

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

be_singlets_metadata %>%
  filter(
    eed_ko_score_normalized > 0.06,
    sample %in% c("CDA1-HLA-pos", "TADA-HLA-pos", "RAPO-HLA-pos")
  ) %>%
  group_by(sample) %>%
  summarize(n_cells_high_eed = n())

# z-score threshold for a 2.5% false positives on upper end
z_score_threshold <- 1.96

# filter variants to enriched variants, that match a potential guide and the BE
# and are not synonymous
predicted_functional_variants <- summary_stats_samples_be_mean_plot %>%
  filter(
    z_score > z_score_threshold,
    !is.na(sgRNA_ID),
    Consequence_single != "synonymous_variant"
  )

length(unique(predicted_functional_variants$POS))
# total of 139 SNV

# edit_lists <-
#   lapply(be_samples, function(x) {
#     predicted_functional_variants %>% filter(editor == x) %>% pull(POS) %>% unique()
#   })
# names(edit_lists) <- be_samples
#
# library(UpSetR)
# UpSetR::upset(
#   fromList(edit_lists),
#   nsets = 8,
#   order.by = "freq",
#   nintersects = NA
# )


# identify potential heterozygous dominant edits:
# If the fraction of cells in which the edit is detected differs from the fraction of UMI with the edit
# detected, then that means that there are cells with transcripts with both wt and mutant at the position.
# If the cell has different edits on each respective allele (from the same guide), then that could explain
# the presence of both wt and mut at a single position in a single cell.
# Therefore, exclude edits that were co-detected within cells with another edit in close distance.
co_occurring_muts <- enriched_edits_cooccurrence %>%
  filter(
    V1 %in% predicted_functional_variants$POS &
      V2 %in% predicted_functional_variants$POS
  ) %>%
  filter(containment_V1 >= 0.25, abs(V1 - V2) <= 20) %>%
  select(editor, V1, containment_V1) %>% unique()

####################################################################################################

# eed has too low transcript detection to hope to find cells with both wt and mutant
# SUZ12 I don't trust because parts also with low coverage and 5' end has super many guides in library
# V7 and V14 have high guide integration, so "heterozygous" looking edits could be cause by guides on different genes
predicted_functional_variants_hetero <- predicted_functional_variants %>%
  # filter for EZH2 and single base editors
  filter(!grepl("V", editor), gene == "EZH2") %>%
  merge(
    co_occurring_muts,
    by.x = c("editor", "POS"),
    by.y = c("editor", "V1"),
    all.x = T
  ) %>%
  filter(is.na(containment_V1))

getProportionCellsHeterozygous <- function(pos, selected_sample) {
  cells_sample <- be_singlets_metadata %>%
    filter(grepl(selected_sample, sample),
           eed_ko_score_normalized > 0.06) %>%
    pull(cellid)
  cells_sample <- cells_sample[cells_sample %in% colnames(mat_alt_samples)]
  cells_mutation <- colnames(mat_alt_samples)[mat_alt_samples[pos, ] != 0]
  cells_mutation <- cells_mutation[cells_mutation %in% cells_sample]
  
  # get the cells with the mutation and also reference transcript
  mat_ref_pos <- mat_ref_alt_samples[pos, cells_mutation] - mat_alt_samples[pos, cells_mutation]
  cells_also_ref <- names(mat_ref_pos)[mat_ref_pos != 0]
  
  umi_alt <- sum(mat_alt_samples[pos, cells_mutation])
  umi_ref <- sum(mat_ref_pos)
  
  return(list(
    length(cells_mutation),
    length(cells_also_ref),
    umi_alt,
    umi_ref
  ))
}

for (i in seq(nrow(predicted_functional_variants_hetero))) {
  if (i %% 10 == 0) {
    cat(i, "\n")
  }
  row <- predicted_functional_variants_hetero[i, ]
  cells_hetero <- getProportionCellsHeterozygous(as.character(row$POS), row$editor)
  predicted_functional_variants_hetero[i, "cells_alt_detected"] <- cells_hetero[[1]]
  predicted_functional_variants_hetero[i, "cells_alt_and_ref_detected"] <- cells_hetero[[2]]
  predicted_functional_variants_hetero[i, "umi_alt"] <- cells_hetero[[3]]
  predicted_functional_variants_hetero[i, "umi_ref"] <- cells_hetero[[4]]
}


write_tsv(
  predicted_functional_variants_hetero,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_enriched_edits_ezh2_sbe_heterozygosity.tsv"
  )
)

custom_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

samples <- c(
  "K562",
  "Ctrl-HLA-pos",
  "CAS9-HLA-pos",
  "V7-HLA-pos",
  "TADA-HLA-pos",
  "CDA1-HLA-pos",
  "V14-HLA-pos",
  "RAPO-HLA-pos"
)
# Dark2 color palette except for replaced gray with blue
samples_color <-
  list(
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#00468b"
  )
names(samples_color) <- samples

####################################################################################################
# Add additional information to enriched mutations: co-occurring enriched mutations in 20bp window or 
# across genes, mean and sd EED KO score, and COSMIC and ClinVar mutations at aa positions

# Co-occurring mutations
# 20bp window
cooccurring_muts_summary_20bp <- enriched_edits_cooccurrence %>%
  filter(
    V1 %in% predicted_functional_variants$POS &
      V2 %in% predicted_functional_variants$POS
  ) %>%
  filter(containment_V1 >= 0.25, abs(V1 - V2) <= 20) %>%
  merge(unique(
    select(predicted_functional_variants, POS, aa_mutation_1letter)
  ), by.x = "V1", by.y = "POS") %>%
  merge(
    unique(
      select(predicted_functional_variants, POS, aa_mutation_1letter)
    ),
    by.x = "V2",
    by.y = "POS",
    suffixes = c("_V1", "_V2")
  ) %>%
  group_by(editor, V1) %>%
  arrange(V1, V2) %>%
  mutate(
    coocurring_genomic_mutations_20bp = paste0(V2, collapse = ";"),
    coocurring_aa_mutations_20bp = paste0(aa_mutation_1letter_V2, collapse = ";")
  ) %>%
  select(V1,
         editor,
         coocurring_genomic_mutations_20bp,
         coocurring_aa_mutations_20bp) %>%
  unique()

# across genes
cooccurring_muts_summary <- enriched_edits_cooccurrence %>%
  filter(
    V1 %in% predicted_functional_variants$POS &
      V2 %in% predicted_functional_variants$POS
  ) %>%
  filter(containment_V1 >= 0.25) %>%
  merge(unique(
    select(predicted_functional_variants, gene, POS, aa_mutation_1letter)
  ), by.x = "V1", by.y = "POS") %>%
  merge(
    unique(
      select(predicted_functional_variants, gene, POS, aa_mutation_1letter)
    ),
    by.x = "V2",
    by.y = "POS",
    suffixes = c("_V1", "_V2")
  ) %>%
  group_by(editor, V1) %>%
  arrange(V1, V2) %>%
  mutate(
    V2 = paste0(gene_V2, "_", V2),
    aa_mutation_1letter_V2 = paste0(gene_V2, "_", aa_mutation_1letter_V2)
  ) %>%
  mutate(
    coocurring_genomic_mutations = paste0(V2, collapse = ";"),
    coocurring_aa_mutations = paste0(aa_mutation_1letter_V2, collapse = ";")
  ) %>%
  select(V1,
         editor,
         coocurring_genomic_mutations,
         coocurring_aa_mutations) %>%
  unique()

predicted_functional_variants <- predicted_functional_variants %>%
  merge(
    cooccurring_muts_summary_20bp,
    by.x = c("editor", "POS"),
    by.y = c("editor", "V1"),
    all.x = T
  ) %>%
  merge(
    cooccurring_muts_summary,
    by.x = c("editor", "POS"),
    by.y = c("editor", "V1"),
    all.x = T
  )

### EED KO score

mat_ref_alt_samples_filtered <- mat_ref_alt_samples[rownames(mat_ref_alt_samples) %in% as.character(predicted_functional_variants$POS),]
mat_alt_samples_filtered <- mat_alt_samples[rownames(mat_alt_samples) %in% as.character(predicted_functional_variants$POS),]

getEedKoEdit <- function(pos, selected_sample) {
  cells_sample <- be_singlets_metadata %>%
    filter(grepl(selected_sample, sample),
           eed_ko_score_normalized > 0.06) %>%
    pull(cellid)
  cells_sample <- cells_sample[cells_sample %in% colnames(mat_alt_samples_filtered)]
  cells_mutation <- colnames(mat_alt_samples_filtered)[mat_alt_samples_filtered[pos, ] != 0]
  cells_mutation <- cells_mutation[cells_mutation %in% cells_sample]
  
  eed_score_cells <- be_singlets_metadata %>% 
    filter(cellid %in% cells_mutation) %>%
    pull(eed_ko_score_normalized)
  
  return(list(
    mean(eed_score_cells),
    sd(eed_score_cells)
  ))
}

for (i in seq(nrow(predicted_functional_variants))) {
  if (i %% 10 == 0) {
    cat(i, "\n")
  }
  row <- predicted_functional_variants[i, ]
  cells_eed_ko <- getEedKoEdit(as.character(row$POS), row$editor)
  predicted_functional_variants[i, "avg_eed_ko_score"] <- cells_eed_ko[[1]]
  predicted_functional_variants[i, "sd_eed_ko_score"] <- cells_eed_ko[[2]]
}


p_avg_eed <- predicted_functional_variants %>%
  ggplot(aes(x = avg_eed_ko_score, y = z_score, color = sample)) +
  geom_point() +
  geom_errorbarh(
    aes(
      xmin = avg_eed_ko_score - sd_eed_ko_score,
      xmax = avg_eed_ko_score + sd_eed_ko_score
    ),
    height = 0.05,
    alpha = 0.3
  ) +
  # geom_smooth(method = "lm") +
  labs(x = "Mean", y = "Z-score", title = "Edit enrichment z-score vs. mean EED KO score with standard deviation error bars") +
  facet_wrap( ~ gene, axes = "all") +
  scale_color_manual(values = samples_color) +
  coord_fixed(ratio = 0.04) +
  custom_theme

pdf(
  file.path(plot_dir, "MF03_enriched_edits_avg_eed_ko_vs_z_score.pdf"),
  height = 9,
  width = 8
)

p_avg_eed

dev.off()

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

# subset to missense variants, because COSMIC is also subsetted to missense variants.
predicted_functional_variants_missense <- predicted_functional_variants %>% filter(Consequence_single == "missense_variant")

predicted_functional_variants_missense_cosmic <- merge(
  predicted_functional_variants_missense,
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
predicted_functional_variants_missense_cosmic <- merge(
  predicted_functional_variants_missense_cosmic,
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

predicted_functional_variants_missense_cosmic <- merge(
  predicted_functional_variants_missense_cosmic,
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
predicted_functional_variants_missense_cosmic <- merge(
  predicted_functional_variants_missense_cosmic,
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

predicted_functional_variants_missense_cosmic_clinvar <- merge(
  predicted_functional_variants_missense_cosmic,
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
predicted_functional_variants_missense_cosmic_clinvar <- merge(
  predicted_functional_variants_missense_cosmic_clinvar,
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

predicted_functional_variants_missense_cosmic_clinvar %>%
  filter(!is.na(clinvar_aa_mutations_at_aa_position))

# check the enriched mutations that have any evidence in either COSMIC or ClinVar
predicted_functional_variants_missense_cosmic_clinvar %>%
  filter(
    !is.na(clinvar_aa_mutations_at_aa_position) |
      !is.na(cosmic_aa_mutations_at_aa_position)
  )

# save table where enriched mutations are combined with COSMIC and ClinVar data.
write_tsv(
  predicted_functional_variants_missense_cosmic_clinvar,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_enriched_edits_missense_cosmic_v100_clinvar.tsv"
  )
)






predicted_functional_variants_missense_cosmic %>%
  filter(cosmic_heterozygous) %>%
  arrange(gene, aa_mutation_1letter) %>%
  select(gene, aa_mutation_1letter, cosmic_n_individuals_aa_mutation) %>%
  unique()
# GENE_SYMBOL aa_mutation n_individuals
# 1         EED       L260F             1
# 2        EZH2       A677V             1
# 3        EZH2       C566Y             1
# 4       SUZ12        G65E             1
# 5       SUZ12        S37L             1

predicted_functional_variants_missense_cosmic %>%
  filter(!is.na(cosmic_n_individuals_aa_mutation)) %>%
  arrange(gene, aa_mutation_1letter) %>%
  select(gene, aa_mutation_1letter, cosmic_n_individuals_aa_mutation) %>%
  unique()
# GENE_SYMBOL aa_mutation n_individuals
# 1 EED         D257N                   1
# 2 EED         L260F                   2
# 3 EZH2        A564T                   1
# 4 EZH2        A576T                   1
# 5 EZH2        A622T                   1
# 6 EZH2        A677V                   4
# 7 EZH2        C560Y                   1
# 8 EZH2        C566R                   1
# 9 EZH2        C566Y                   3
# 10 EZH2        G623D                   1
# 11 EZH2        H129R                   1
# 12 EZH2        V621A                   1
# 13 EZH2        W624R                   1
# 14 SUZ12       E421K                   1
# 15 SUZ12       G65E                    1
# 16 SUZ12       S37L                    1
# 17 SUZ12       S45G                    1
# 18 SUZ12       S56F                    1
# 19 SUZ12       Y480C                   1

cosmic_db %>% select(GENE_SYMBOL, aa_mutation) %>% unique() %>% nrow()
# 968
predicted_functional_variants_missense %>% select(gene, aa_mutation_1letter) %>% unique() %>% nrow()
# 124
predicted_functional_variants_missense_cosmic %>% filter(!is.na(cosmic_n_individuals_aa_mutation)) %>% select(gene, aa_mutation) %>% unique() %>% nrow()
# 19
# 15% of all enriched missense mutations are in COSMIC

cosmic_db %>% select(GENE_SYMBOL, aa_mutation) %>%
  filter(GENE_SYMBOL == "EZH2") %>% unique() %>% nrow()
# 540
predicted_functional_variants_missense %>% select(gene, aa_mutation_1letter) %>%
  filter(gene == "EZH2") %>% unique() %>% nrow()
# 52
predicted_functional_variants_missense_cosmic %>% filter(!is.na(cosmic_n_individuals_aa_mutation)) %>% select(gene, aa_mutation) %>%
  filter(gene == "EZH2") %>% unique() %>% nrow()
# 11
# 21% of enriched EZH2 missense mutations are in COSMIC


# Test for enrichment of enriched mutated amino acids in COSMIC DB
position_cosmic <- cosmic_db %>% select(GENE_SYMBOL, aa_position) %>%
  filter(GENE_SYMBOL == "EZH2") %>% pull(aa_position) %>% unique()

position_enriched <- predicted_functional_variants_missense %>% select(gene, aa_position) %>%
  filter(gene == "EZH2") %>% pull(aa_position) %>% unique()

# Count how many of your mutations are in the database
a <- sum(position_enriched %in% position_cosmic)  # in both your list and database
b <- sum(!position_enriched %in% position_cosmic) # in your list but not in database

# Assume the background is all mutations in the database plus all your mutations not in database
# c = database mutations not in your list
c <- sum(!position_cosmic %in% position_enriched)
# d = all other mutations (background)
d <- 746 - length(unique(c(position_cosmic, position_enriched)))

# Build table
contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
rownames(contingency_table) <- c("Enriched_positions", "Other")
colnames(contingency_table) <- c("COSMIC_positions", "Other")

contingency_table
# COSMIC_positions Other
# Enriched_positions               20    25
# Other                           329   372

fisher.test(contingency_table)

# Fisher's Exact Test for Count Data
#
# data:  contingency_table
# p-value = 0.7608
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.4670551 1.7318666
# sample estimates:
# odds ratio
#  0.9046962


# Test for enrichment of enriched mutated amino acids in COSMIC DB - only for haematopoietic malignancies
position_cosmic <- cosmic_db_haema %>% select(GENE_SYMBOL, aa_position) %>%
  filter(GENE_SYMBOL == "EZH2") %>% pull(aa_position) %>% unique()

position_enriched <- predicted_functional_variants_missense %>% select(gene, aa_position) %>%
  filter(gene == "EZH2") %>% pull(aa_position) %>% unique()

# Count how many of your mutations are in the database
a <- sum(position_enriched %in% position_cosmic)  # in both your list and database
b <- sum(!position_enriched %in% position_cosmic) # in your list but not in database

# Assume the background is all mutations in the database plus all your mutations not in database
# c = database mutations not in your list
c <- sum(!position_cosmic %in% position_enriched)
# d = all other mutations (background)
d <- 746 - length(unique(c(position_cosmic, position_enriched)))

# Build table
contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
rownames(contingency_table) <- c("Enriched_positions", "Other")
colnames(contingency_table) <- c("COSMIC_positions", "Other")

contingency_table
# COSMIC_positions Other
# Enriched_positions                9    36
# Other                           117   584

fisher.test(contingency_table)

# Fisher's Exact Test for Count Data
#
# data:  contingency_table
# p-value = 0.5405
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5142848 2.7296229
# sample estimates:
# odds ratio
#   1.247468
