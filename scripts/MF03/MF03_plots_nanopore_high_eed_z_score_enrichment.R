#!/usr/bin/env Rscript

# DESCRIPTION
# See MF03_14_nanopore_high_eed_z_score_enrichment.R

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

summary_stats_samples_be_mean_plot <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_high_eed_ko_cells_z_score.tsv"
  )
)

# z-score threshold for a 2.5% false positives on upper end
z_score_threshold <- 1.96

be_samples <- c("V14", "V7", "TADA", "RAPO", "CDA1")
genes <- c("EED", "EZH2", "SUZ12")

summary_stats_samples_be <- summary_stats_samples %>%
  filter(!editor %in% c("K562", "Ctrl", "CAS9", "K562-Ctrl"))

custom_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

variant_effect_colors = list(
  "synonymous_variant" = "#E69F00",
  "missense_variant" = "#56B4E9",
  "stop_gained" = "#009E73",
  "3_prime_UTR_variant" = "#F0E442",
  "5_prime_UTR_variant" = "#0072B2",
  "start_lost" = "#CC79A7",
  "splice_region_variant" = "#D55E00"
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

editor_colors <- samples_color
names(editor_colors) <- gsub("-HLA-pos", "", names(editor_colors))

# calculate z-score only on edits that show in the plot later, i.e. only on edits that pass the minimum criteria
# n_cells_sample_edit >= 10 means at least 10 cells in this sample have the edit detected. 
# fraction_cells_sample > 0.01 means that at least 1% of cells that have transcript coverage at that position
# have the edit detected (considers different sequencing depth of genes)
summary_stats_samples_be_filtered <- summary_stats_samples_be %>%
  filter(n_cells_sample_edit >= 10) %>%
  filter(fraction_cells_sample > 0.01)

# Q1: is the background noise different along the protein sequence?
p1 <- summary_stats_samples_be_filtered %>%
  filter(Consequence_single == "synonymous_variant") %>%
  ggplot(aes(x = gene, y = log2_fraction_umi_sample_over_fraction_umi_k562, fill = gene)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("EED", "EZH2"), c("EZH2", "SUZ12"), c("SUZ12", "EED")),
                     label = "p.format") +
  ggtitle("Enrichment of synonymous variants across genes") +
  custom_theme

# Q3: is the log2FC different for each base editor?
p2 <- summary_stats_samples_be_filtered %>%
  filter(Consequence_single == "synonymous_variant") %>%
  ggplot(aes(x = editor, y = log2_fraction_umi_sample_over_fraction_umi_k562, fill = editor)) +
  scale_fill_manual(values = editor_colors) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  # facet_wrap(~ gene) +
  stat_compare_means(
    method = "t.test",
    comparisons = combn(
      be_samples,
      2,
      FUN = function(x)
        x,
      simplify = FALSE
    ),
    label = "p.format"
  ) +
  ggtitle("Enrichment of synonymous variants across editors") +
  custom_theme

# Q4: is the log2FC different for missense vs synonymous variants?
p3 <- summary_stats_samples_be_filtered %>%
  filter(Consequence_single %in% c("synonymous_variant", "missense_variant")) %>%
  ggplot(
    aes(x = Consequence_single, y = log2_fraction_umi_sample_over_fraction_umi_k562, fill = Consequence_single)
  ) +
  scale_fill_manual(values = variant_effect_colors) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  stat_compare_means(method = "t.test",
                     comparisons = list(c(
                       "synonymous_variant", "missense_variant"
                     )),
                     label = "p.format") +
  ggtitle("Enrichment of synonymous vs missense variants") +
  custom_theme

pdf(
  file.path(
    plot_dir,
    "MF03_high_eed_ko_score_edit_enrichment_comparisons.pdf"
  ),
  height = 6,
  width = 6
)

p1
p2
p3

dev.off()

####################################################################################################

# compare the proportion of synonymous and of "impossible" edits between enriched and not enriched mutations.
pp3 <- summary_stats_samples_be_mean_plot %>%
  mutate(
    enriched = z_score >= z_score_threshold,
    synonymous = Consequence_single == "synonymous_variant"
  ) %>%
  ungroup() %>%
  select(enriched, POS, base_change_possible, synonymous) %>%
  unique() %>%
  ggplot(aes(x = enriched)) +
  geom_bar(aes(fill = synonymous), position = "fill") +
  ggtitle("Background distribution: synonymous variants, by editor and gene") +
  ylab("edits") +
  xlab("z-score >= 2") +
  custom_theme

pp4 <- summary_stats_samples_be_mean_plot %>%
  mutate(enriched = z_score >= z_score_threshold) %>%
  ungroup() %>%
  select(enriched, POS, base_change_possible) %>%
  unique() %>%
  ggplot(aes(x = enriched)) +
  geom_bar(aes(fill = !base_change_possible), position = "fill") +
  xlab("z-score >= 2") +
  ylab("edits") +
  guides(fill = guide_legend(title = "BE incompatible")) +
  custom_theme

pp5 <- summary_stats_samples_be_mean_plot %>%
  mutate(enriched = z_score >= z_score_threshold) %>%
  ungroup() %>%
  select(enriched, POS, sgRNA_ID) %>%
  unique() %>%
  ggplot(aes(x = enriched)) +
  geom_bar(aes(fill = is.na(sgRNA_ID)), position = "fill") +
  xlab("z-score >= 2") +
  ylab("edits") +
  guides(fill = guide_legend(title = "BE+library \nincompatible")) +
  custom_theme


pdf(
  file.path(
    plot_dir,
    "MF03_high_eed_ko_score_enriched_edits_comparisons.pdf"
  ),
  height = 4,
  width = 8
)

(pp3 | pp4 | pp5)

dev.off()


####################################################################################################
# To plot the whole protein sequence
data_dummy <-
  data.frame(list(
    aa_position = c(0, 441, 0, 746, 0, 739),
    gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12"),
    z_score = c(0, 0, 0, 0, 0, 0)
  ))

p6 <- summary_stats_samples_be_mean_plot %>%
  mutate(label = ifelse(
    z_score >= z_score_threshold,
    paste0(aa_mutation_1letter, " ", editor),
    ""
  )) %>%
  mutate(label = gsub("V14", "EC", label)) %>%
  ggplot(aes(x = aa_position, y = z_score)) +
  geom_blank(data = data_dummy) +
  geom_hline(yintercept = z_score_threshold,
             linetype = 2,
             alpha = 0.3) +
  # show edits below z-score cutoff as gray dots
  geom_point(
    data = . %>% dplyr::filter(label == ""),
    color = "gray",
    alpha = 0.2
  ) +
  # color edits above z-score cutoff by mutation consequence
  # transparent dots = base change not possible by editor
  geom_point(data = . %>% dplyr::filter(label != ""),
             aes(color = Consequence_single, alpha = !is.na(sgRNA_ID))) +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_size_manual(values = c(1.5, 2.5)) +
  scale_color_manual(values = variant_effect_colors) +
  facet_wrap(~ gene,
             scales = "free_x",
             axes = "all",
             ncol = 1) +
  ylab("Z-score") +
  ggtitle("Z-score of edit enrichment", 
          subtitle = "Background distribution: synonymous variants, by editor and gene") +
  custom_theme

p6_labels <- p6 +
  ggrepel::geom_text_repel(
    data = . %>% filter(!is.na(sgRNA_ID), label != ""),
    aes(label = label, color = Consequence_single),
    size = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  )


# do the same but don't show the variants that are not possible by BE and synonymous variants
p8 <- summary_stats_samples_be_mean_plot %>%
  filter(Consequence_single != "synonymous_variant") %>%
  filter(!is.na(sgRNA_ID)) %>%
  mutate(label = ifelse(
    z_score >= z_score_threshold,
    paste0(aa_mutation_1letter, " ", editor),
    ""
  )) %>%
  mutate(label = gsub("V14", "EC", label)) %>%
  ggplot(aes(x = aa_position, y = z_score)) +
  geom_blank(data = data_dummy) +
  geom_hline(yintercept = z_score_threshold,
             linetype = 2,
             alpha = 0.3) +
  # show edits below z-score cutoff as gray dots
  geom_point(
    data = . %>% dplyr::filter(label == ""),
    color = "gray",
    alpha = 0.2
  ) +
  # color edits above z-score cutoff by mutation consequence
  geom_point(data = . %>% dplyr::filter(label != ""),
             aes(color = Consequence_single)) +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_size_manual(values = c(1.5, 2.5)) +
  scale_color_manual(values = variant_effect_colors) +
  facet_wrap(~ gene,
             scales = "free_x",
             axes = "all",
             ncol = 1) +
  ylab("Z-score") +
  ggtitle("Z-score of edit enrichment", subtitle = "Background distribution: synonymous variants, by editor and gene \nFiltered for non-synonymous mutations and BE compatible edits") +
  custom_theme

p8_labels <- p8 +
  ggrepel::geom_text_repel(
    data = . %>% filter(!is.na(sgRNA_ID), label != ""),
    aes(label = label, color = Consequence_single),
    size = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  )


# same thing as above but color by editor. No need to mention editor in label.
p9 <- summary_stats_samples_be_mean_plot %>%
  filter(Consequence_single != "synonymous_variant") %>%
  filter(!is.na(sgRNA_ID)) %>%
  mutate(label = ifelse(z_score >= z_score_threshold, aa_mutation_1letter, "")) %>%
  mutate(label = gsub("V14", "EC", label)) %>%
  ggplot(aes(x = aa_position, y = z_score)) +
  geom_blank(data = data_dummy) +
  geom_hline(yintercept = z_score_threshold,
             linetype = 2,
             alpha = 0.3) +
  # show edits below z-score cutoff as gray dots
  geom_point(
    data = . %>% dplyr::filter(label == ""),
    color = "gray",
    alpha = 0.2
  ) +
  # color edits above z-score cutoff by mutation consequence
  geom_point(data = . %>% dplyr::filter(label != ""), aes(color = sample)) +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_size_manual(values = c(1.5, 2.5)) +
  scale_color_manual(values = samples_color) +
  facet_wrap(~ gene,
             scales = "free_x",
             axes = "all",
             ncol = 1) +
  ylab("Z-score") +
  ggtitle("Z-score of edit enrichment", subtitle = "Background distribution: synonymous variants, by editor and gene \nFiltered for non-synonymous mutations and BE compatible edits") +
  custom_theme

p9_labels <- p9 +
  ggrepel::geom_text_repel(
    data = . %>% filter(!is.na(sgRNA_ID), label != ""),
    aes(label = label, color = sample),
    size = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  )


# do the same but ONLY show the variants that are not possible by BE and synonymous variants
p10 <- summary_stats_samples_be_mean_plot %>%
  filter(Consequence_single == "synonymous_variant" |
           is.na(sgRNA_ID)) %>%
  mutate(color_category = ifelse(
    !base_change_possible,
    "BE incompatible edit",
    ifelse(
      is.na(sgRNA_ID),
      "BE+library incompatible edit",
      ifelse(
        synonymous_variant_coocurring_w_non_synonymous,
        "synonymous_variant co-occurring\n with non-synonymous",
        "synonymous_variant"
      )
    )
  )) %>%
  mutate(label = ifelse(
    z_score >= z_score_threshold,
    paste0(aa_mutation_1letter, " ", editor),
    ""
  )) %>%
  mutate(label = gsub("V14", "EC", label)) %>%
  ggplot(aes(x = aa_position, y = z_score)) +
  geom_blank(data = data_dummy) +
  geom_hline(yintercept = z_score_threshold,
             linetype = 2,
             alpha = 0.3) +
  # show edits below z-score cutoff as gray dots
  geom_point(
    data = . %>% dplyr::filter(label == ""),
    color = "gray",
    alpha = 0.2
  ) +
  geom_point(data = . %>% dplyr::filter(label != ""), aes(color = color_category)) +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_size_manual(values = c(1.5, 2.5)) +
  scale_color_manual(values = c("#E66100", "#5D3A9B", "#D4B87B", "#56B4F1")) +
  facet_wrap(~ gene,
             scales = "free_x",
             axes = "all",
             ncol = 1) +
  ylab("Z-score") +
  ggtitle("Z-score of edit enrichment", subtitle = "Background distribution: synonymous variants, by editor and gene \nFiltered for synonymous mutations and BE and library incompatible edits") +
  custom_theme

p10_labels <- p10 +
  ggrepel::geom_text_repel(
    data = . %>% filter(label != "" &
                          color_category == "synonymous_variant"),
    aes(label = label, color = color_category),
    size = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  )


pdf(
  file.path(
    plot_dir,
    "MF03_high_eed_ko_score_edit_enrichment_z_score.pdf"
  ),
  height = 9,
  width = 9
)

p6
p6_labels

p8
p8_labels

p9
p9_labels

p10
p10_labels

dev.off()

