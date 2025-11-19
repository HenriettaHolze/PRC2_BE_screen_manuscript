#!/usr/bin/env Rscript

# DESCRIPTION
# Plot mutations enriched in cells with high EED KO score, independent from guide detection.
# Goal: identify edits that have a sufficiently high enrichment over the control sample
#
# Approach: calculate a z-score relative to the synonymous mutations. Almost all synonymous variants
#     will be noise (sequencing errors or passenger mutations). Alternative approach: z-score to the
#     "impossible" mutations(e.g. C->G mutations). But they differ a lot for different editors
#     (dual vs. single base editor).
# Question: how "detailed" does the background distribution of nois have to be?
# Q1: is the log2FC different for different genes?
# Q2: is the background noise different along the protein sequence?
# Q3: is the log2FC different for each base editor?
# Q4: is the log2FC different for missense vs synonymous variants?

# The real question is why a difference would be expected across editors and genes/position.
# For genes, e.g. EED has lower sequencing depth and more variants that are not detected at all in K562,
# so pseudo count has a higher effect.

library(tidyverse)
library(ggpubr)
library(seqinr)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

## Load data
summary_stats_samples <- read_tsv(paste0(
  result_dir,
  "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells.tsv"
))

be_samples <- c("V14", "V7", "TADA", "RAPO", "CDA1")
genes <- c("EED", "EZH2", "SUZ12")

summary_stats_samples_be <- summary_stats_samples %>%
  filter(!editor %in% c("K562", "Ctrl", "CAS9", "K562-Ctrl"))

# calculate z-score only on edits that show in the plot later, i.e. only on edits that pass the minimum criteria
# n_cells_sample_edit >= 10 means at least 10 cells in this sample have the edit detected. 
# fraction_cells_sample > 0.01 means that at least 1% of cells that have transcript coverage at that position
# have the edit detected (considers different sequencing depth of genes)
summary_stats_samples_be_filtered <- summary_stats_samples_be %>%
  filter(n_cells_sample_edit >= 10) %>%
  filter(fraction_cells_sample > 0.01)

# check number of synonymous variants per sample to use for z-score
summary_stats_samples_be_filtered %>%
  group_by(gene, editor, Consequence_single) %>%
  summarize(n()) %>%
  filter(Consequence_single == "synonymous_variant")

# Q1: is the background noise different along the protein sequence?
# ** see MF03_plots_nanopore_high_eed_z_score_enrichment.R **
# Q3: is the log2FC different for each base editor?
# ** see MF03_plots_nanopore_high_eed_z_score_enrichment.R **
# Q4: is the log2FC different for missense vs synonymous variants?
# ** see MF03_plots_nanopore_high_eed_z_score_enrichment.R **

# calculate the mean and standard deviation per gene and editor, irrespective of position
mean_sd_gene_editor <- summary_stats_samples_be_filtered %>%
  filter(Consequence_single == "synonymous_variant",!is.na(aa_position)) %>%
  group_by(editor, gene) %>%
  summarize(
    avg = mean(log2_fraction_umi_sample_over_fraction_umi_k562),
    sd = sd(log2_fraction_umi_sample_over_fraction_umi_k562)
  )

mean_sd_gene_editor_all_muts <- summary_stats_samples_be_filtered %>%
  filter(!is.na(aa_position)) %>%
  group_by(editor, gene) %>%
  summarize(
    avg_all_muts = mean(log2_fraction_umi_sample_over_fraction_umi_k562),
    sd_all_muts = sd(log2_fraction_umi_sample_over_fraction_umi_k562)
  )

# add the average and sd of synonymous mutations to the dataframe
summary_stats_samples_be_mean <- merge(
  summary_stats_samples_be_filtered,
  mean_sd_gene_editor,
  by = c("gene", "editor"),
  all.x = T
)

summary_stats_samples_be_mean <- merge(
  summary_stats_samples_be_mean,
  mean_sd_gene_editor_all_muts,
  by = c("gene", "editor"),
  all.x = T
)

# calculate Z-score with respect to synonymous mutations for both rolling mean and mean for
# gene x editor combination
summary_stats_samples_be_mean_plot <- summary_stats_samples_be_mean %>%
  filter(!is.na(aa_position)) %>%
  mutate(z_score = (log2_fraction_umi_sample_over_fraction_umi_k562 - avg) / sd) %>%
  mutate(
    z_score_all_muts = (
      log2_fraction_umi_sample_over_fraction_umi_k562 - avg_all_muts
    ) / sd_all_muts
  )

##########
# Match edits to possible guides, same way as I annotated edits to guides detected in the cell.
# Combine all edits to all guides and all editors and then filter for the pairs that match.

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

# Load guide metadata (cut position and strand).
guide_metadata <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12.tsv"
    )
  )
guide_metadata$sgRNA_ID <- gsub("_", "-", guide_metadata$sgRNA_ID)

pos_subset <- unique(summary_stats_samples_be_mean_plot$POS)

# annotate gene to position
vcf <- vcf %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(gene = ifelse(POS < 40000000, "SUZ12", ifelse(POS < 100000000, "EED", "EZH2"))) %>%
  filter(POS %in% pos_subset)
guide_metadata_vcf <- merge(guide_metadata, vcf, by = "gene", all = T)

# Cut site is annotated as guide position 17 for guides on (-) strans and position 18 for guides
# on (+) strand (by [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public).
# I checked with UCSC for some guides, that position - cut position + 18 for guide on (+)
# and cut position - position + 17 for guide on (-)
# gives the position of the edit on the guide sequence (1-based).
guide_metadata_vcf <- guide_metadata_vcf %>% mutate(
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

guide_metadata_vcf <- guide_metadata_vcf %>%
  mutate(detected_edit = ifelse(
    Strand.of.sgRNA == "+",
    paste0(REF, " -> ", ALT),
    paste0(getCompl[REF], " -> ", getCompl[ALT])
  )) %>%
  mutate(edit_on_guide_window = (
    !is.na(position_on_guide) &
      position_on_guide < 20 &
      position_on_guide > -5
  ))

guide_metadata_vcf <- guide_metadata_vcf %>%
  filter(edit_on_guide_window)

editors <- c("V7", "V14", "CDA1", "RAPO", "TADA")
genes <- c("EZH2", "SUZ12", "EED")
guide_metadata_vcf <- merge(
  guide_metadata_vcf,
  expand.grid(editor = editors, gene = genes),
  by = "gene",
  all = T
)

# annotate whether the edit matches the capacity of the base editor or not
guide_metadata_vcf <- guide_metadata_vcf %>%
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

allowed_edits_editor <- guide_metadata_vcf %>%
  filter(correct_edits) %>%
  select(POS, editor, sgRNA_ID) %>%
  arrange(POS, editor, sgRNA_ID) %>%
  group_by(POS, editor) %>%
  summarise(sgRNA_ID = paste(sgRNA_ID, collapse = ","))

##########

summary_stats_samples_be_mean_plot <- summary_stats_samples_be_mean_plot %>%
  mutate(aa_mutation = gsub("%3D", "", aa_mutation))

# convert 3 letter code into 1 letter code
summary_stats_samples_be_mean_plot <- summary_stats_samples_be_mean_plot %>%
  mutate(aa1 = a(gsub("[0-9]+.*", "", aa_mutation)), aa2 = a(gsub(".*[0-9]+", "", aa_mutation))) %>%
  mutate(aa2 = ifelse(grepl("Ter", aa_mutation), "*", aa2)) %>%
  mutate(aa2 = replace_na(aa2, "")) %>%
  mutate(aa_mutation_1letter = paste0(aa1, aa_position, aa2)) %>%
  select(-aa2, -aa1)

summary_stats_samples_be_mean_plot <- summary_stats_samples_be_mean_plot %>%
  select(-z_score_all_muts, -avg_all_muts, -sd_all_muts)

summary_stats_samples_be_mean_plot <- merge(
  summary_stats_samples_be_mean_plot,
  allowed_edits_editor,
  by = c("POS", "editor"),
  all.x = T
)



# z-score threshold for a 2.5% false positives on upper end
z_score_threshold <- 1.96


####################################################################################################
# Check if synonymous variants that are enriched and possible base edits are detected together with
# non-synonymous (missense, splice, stop) mutations in the same cell.
# If cooccurrence is close to 100%, then this reveals if on the same allele is a more likely functional mutation.
# Across short distance (20bp) they are likely caused by the same guide, across longer distances it is from different guides.
# Calculate jaccard index and the containment index (how many cells with the synonymous mutation also
# have the other mutation).

# cellsnp data that was filtered for cells that can be merged with scRNA-seq data
cellsnp_res_dir <-
  paste0(result_dir, "MF03/nanopore/")
# read in reference and alternative UMI counts for cells that were demultiplexed into samples
mat_alt_samples <-
  Matrix::readMM(paste0(
    result_dir,
    "MF03/nanopore/mat_alt_samples_sicelore_cellsnp.mtx"
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
colnames(mat_alt_samples) <- cells
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

# jaccard index of neighbouring mutations
position_close_pairs_samples <- lapply(be_samples, function(selected_sample) {
  # filter positions that pass filters, i.e. high enrichment, possible base change
  positions_pass <- summary_stats_samples_be_mean_plot %>%
    filter(
      z_score >= z_score_threshold,
      base_change_possible,
      !is.na(sgRNA_ID),
      editor == selected_sample
    ) %>%
    pull(POS) %>% unique()
  position_close_pairs <- combn(positions_pass, 2)
  return(position_close_pairs)
})
names(position_close_pairs_samples) <- be_samples

position_close_pairs_samples <- t(bind_rows(position_close_pairs_samples)) %>%
  as.data.frame() %>%
  rownames_to_column("editor") %>%
  mutate(editor = gsub("\\..*", "", editor))
nrow(position_close_pairs_samples)
position_close_pairs_samples$jaccard <- 0

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# prefilter matrix
mat_alt_samples_filtered <- mat_alt_samples[rownames(mat_alt_samples) %in% as.character(c(
  position_close_pairs_samples$V1,
  position_close_pairs_samples$V2
)), ]

for (i in seq(nrow(position_close_pairs_samples))) {
  row <- position_close_pairs_samples[i, ]
  selected_sample <- row[["editor"]]
  cells_sample <- be_singlets_metadata %>% 
    filter(grepl(selected_sample, sample), eed_ko_score_normalized > 0.06) %>% 
    pull(cellid)
  # get cells for the sample with either mutation
  cells_pos1 <- colnames(mat_alt_samples_filtered)[mat_alt_samples_filtered[as.character(row[["V1"]]), ] != 0]
  # subset to cells from sample
  cells_pos1 <- cells_pos1[cells_pos1 %in% cells_sample]
  cells_pos2 <- colnames(mat_alt_samples_filtered)[mat_alt_samples_filtered[as.character(row[["V2"]]), ] != 0]
  cells_pos2 <- cells_pos2[cells_pos2 %in% cells_sample]
  # calculate and save jaccard
  position_close_pairs_samples[i, "jaccard"] <- jaccard(cells_pos1, cells_pos2)
  position_close_pairs_samples[i, "containment_V1"] <- length(intersect(cells_pos1, cells_pos2)) / length(cells_pos1)
  position_close_pairs_samples[i, "containment_V2"] <- length(intersect(cells_pos2, cells_pos1)) / length(cells_pos2)
}

position_close_pairs_samples_consequence <- position_close_pairs_samples %>%
  merge(
    summary_stats_samples_be_mean_plot %>% select(editor, POS, Consequence_single, aa_position),
    by.x = c("editor", "V1"),
    by.y = c("editor", "POS")
  ) %>%
  merge(
    summary_stats_samples_be_mean_plot %>% select(editor, POS, Consequence_single, aa_position),
    by.x = c("editor", "V2"),
    by.y = c("editor", "POS"),
    suffixes = c("_V1", "_V2")
  )

df_swapped <- position_close_pairs_samples_consequence
df_swapped <- df_swapped[, c(
  "editor",
  "V1",
  "V2",
  "jaccard",
  "containment_V1",
  "containment_V2",
  "Consequence_single_V2",
  "Consequence_single_V1",
  "aa_position_V1",
  "aa_position_V2"
)]
names(df_swapped) <- names(df)  # rename columns back to original order

# Combine original + swapped
position_close_pairs_samples_consequence_abba <- rbind(position_close_pairs_samples_consequence)

# check if there are any enriched synonymous variants co-occurring and on the same codon
position_close_pairs_samples_consequence_abba %>%
  filter(
    Consequence_single_V1 == "synonymous_variant" &
      Consequence_single_V2 == "synonymous_variant" &
      abs(V1 - V2) <= 20 & aa_position_V1 == aa_position_V2
  )
# none

# check if there are synonymous mutations co-occurring with non-synonymous in close distance, i.e. same guide
synonymous_position_coocurring_w_non_synonymous_same_guide <- position_close_pairs_samples_consequence_abba %>%
  filter(
    Consequence_single_V1 == "synonymous_variant" &
      Consequence_single_V2 != "synonymous_variant" &
      abs(V1 - V2) <= 20
  ) %>%
  group_by(V1, editor) %>%
  slice_max(order_by = containment_V1,
            n = 1,
            with_ties = F) %>%
  filter(containment_V1 >= 0.25) %>% pull(V1)
length(synonymous_position_coocurring_w_non_synonymous_same_guide)

summary_stats_samples_be_mean_plot %>%
  filter(
    z_score >= z_score_threshold,
    base_change_possible,
    !is.na(sgRNA_ID),
    Consequence_single == "synonymous_variant"
  ) %>%
  nrow()
# 21/65 cases (32%) (mutation x BE combination) there is another mutation co-occurring from same guide

# check for synonymous variants cooccurring with non-synonymous across genes
synonymous_position_coocurring_w_non_synonymous <- position_close_pairs_samples_consequence_abba %>%
  filter(
    Consequence_single_V1 == "synonymous_variant" &
      Consequence_single_V2 != "synonymous_variant"
  ) %>%
  group_by(V1, editor) %>%
  slice_max(order_by = containment_V1,
            n = 1,
            with_ties = F) %>%
  filter(containment_V1 >= 0.25) %>% select(V1, editor) %>%
  mutate(synonymous_variant_coocurring_w_non_synonymous = T)
nrow(synonymous_position_coocurring_w_non_synonymous)
# 32/65 cases (49%) (mutation x BE combination) there is another mutation co-occurring

# add information to dataframe
summary_stats_samples_be_mean_plot <- merge(
  summary_stats_samples_be_mean_plot,
  synonymous_position_coocurring_w_non_synonymous,
  by.x = c("POS", "editor"),
  by.y = c("V1", "editor"),
  all.x = T
) %>%
  mutate(synonymous_variant_coocurring_w_non_synonymous = replace_na(synonymous_variant_coocurring_w_non_synonymous, FALSE))

####################################################################################################
# Save data
write_tsv(
  summary_stats_samples_be_mean_plot,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells_z_score.tsv"
  )
)

write_tsv(
  position_close_pairs_samples_consequence_abba,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_enriched_edit_cooccurrence.tsv"
  )
)

# 
# summary_stats_samples_be_mean_plot_previous <- read_tsv(
#   paste0(
#     result_dir,
#     "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells_z_score.tsv"
#   )
# )
# 
# colnames(summary_stats_samples_be_mean_plot)
# 
# custom_theme <- theme_bw() +
#   theme(
#     panel.grid = element_blank(),
#     legend.position = "bottom",
#     strip.background = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(linewidth = 0.5)
#   )
# 
# summary_stats_samples_be_mean_plot %>%
#   filter(z_score > 1) %>%
#   ggplot(aes(
#     x = aa_position,
#     y = n_umi_sample_coverage,
#     color = !is.na(sgRNA_ID)
#   )) +
#   geom_point(size = 0.5, alpha = 0.5) +
#   facet_grid(editor ~ gene, scales = "free_x", axes = "all") +
#   custom_theme
# 
# summary_stats_samples_be_mean_plot %>%
#   ggplot(aes(
#     x = z_score,
#     y = n_umi_sample_coverage,
#     color = !is.na(sgRNA_ID)
#   )) +
#   geom_point(size = 0.5, alpha = 0.2) +
#   facet_grid(editor ~ gene, axes = "all") +
#   custom_theme
# 
# 
# summary_stats_samples_be_mean_plot_previous_filtered <- summary_stats_samples_be_mean_plot_previous %>% filter(z_score > 1.96, editor == "V7", gene == "SUZ12")
# summary_stats_samples_be_mean_plot_filtered <- summary_stats_samples_be_mean_plot %>% filter(z_score > 1.96, editor == "V7", gene == "SUZ12") %>% as_tibble()
# 
# 
# nrow(summary_stats_samples_be_mean_plot_previous_filtered)
# nrow(summary_stats_samples_be_mean_plot_filtered)
# 
# 
# tmp <- merge(summary_stats_samples_be_mean_plot_previous_filtered, summary_stats_samples_be_mean_plot_filtered, by=c("POS", "editor", "sgRNA_ID"), suffixes = c("_previous", ""))
# 
# ggplot(tmp, aes(z_score, z_score_previous, color = !is.na(sgRNA_ID), size = aa_position)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_hline(yintercept = 1.96) +
#   geom_vline(xintercept = 1.96) +
#   ggrepel::geom_text_repel(aes(label = aa_mutation_1letter_previous), max.overlaps = Inf)
