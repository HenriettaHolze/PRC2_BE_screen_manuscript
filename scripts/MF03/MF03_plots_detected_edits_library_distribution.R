#!/usr/bin/env Rscript

# DESCRIPTION
# Plot distribution of edits that were induced by base editors across protein sequence.
# Do not filter for whether they are likely functional.
# Also plot overview over kinds of edits for each base editor.

library(tidyverse)
library(patchwork)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

####################################################################################################
# Load data
be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )
be_singlets_metadata <- be_singlets@meta.data
be_singlets_metadata$cell_id <- rownames(be_singlets_metadata)

sample_alt_ref_df_guides_anno <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation.tsv"
  )
)
head(sample_alt_ref_df_guides_anno)

summary_stats_samples_guides_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_edit_frequency_be_guides_10_cells.tsv"
  )
)
head(summary_stats_samples_guides_df)

# Load guide metadata (cut position and strand).
guide_metadata <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12.tsv"
    )
  )
guide_metadata$sgRNA_ID <- gsub("_", "-", guide_metadata$sgRNA_ID)

guide_counts_library <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_guide_library/bartab/results_18bp/counts/all_counts_combined.tsv"
    )
  )
colnames(guide_counts_library) <- c("sgRNA_ID", "count_bartab_18bp")
guide_counts_library$sgRNA_ID <- gsub("_", "-", guide_counts_library$sgRNA_ID)
guide_metadata <- merge(guide_metadata, guide_counts_library, by="sgRNA_ID", all=T)

guide_counts_samples <- read_tsv(paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts.tsv"))
guide_counts_ctrl <- guide_counts_samples %>% 
  mutate(sgRNA_ID = gsub("_", "-", barcode)) %>%
  filter(bc.umi.count >= 3) %>% 
  filter(sample == "Ctrl-HLA-pos") %>%
  group_by(sgRNA_ID) %>%
  summarize(n_cells_ctrl_hla_pos = n())

guide_metadata <- merge(guide_metadata, guide_counts_ctrl, by = "sgRNA_ID", all=T)

variant_effect_colors = list(
  "synonymous_variant" = "#E69F00",
  "missense_variant" = "#56B4E9",
  "stop_gained" = "#009E73",
  "3_prime_UTR_variant" = "#F0E442",
  "5_prime_UTR_variant" = "#0072B2",
  "start_lost" = "#CC79A7",
  "splice_region_variant" = "#D55E00"
)

####################################################################################################
# Filter edits to base editor samples, matching guide and base editor in cell, detected in at least
# 3 cells in sample matching that guide, not heterozygous in K562

cells_high_eed_ko_score <- dplyr::filter(be_singlets_metadata, eed_ko_score_normalized > 0.06) %>%
  pull(cell_id)

be_samples <- c("CDA1-HLA-pos",
                "TADA-HLA-pos",
                "RAPO-HLA-pos",
                "V7-HLA-pos",
                "V14-HLA-pos")

# this table only contains cells with a guide detected.
# Column "guide" contains the guide detected with >= 3 UMI
# Column "guide_clone" contains the information to which guide clone the cell belongs.
# So there can be multiple rows per cell, one for each guide that it contains, 
# as well as every mutation in the cell that matches the respective guide
edits_detected_plot <- sample_alt_ref_df_guides_anno %>%
  # filter only for base editor samples, not cas9 or controls
  filter(sample %in% be_samples) %>%
  # remove heterozygous positions
  filter(k562_alt_freq < 0.1 | is.na(k562_alt_freq)) %>%
  # filter edits that match the guide detected in the cell, and the BE
  filter(
    !is.na(position_on_guide) &
      position_on_guide < 20 &
      position_on_guide > -5 &
      correct_edits == TRUE
  ) %>%
  # count number of cells per BE/guide/edit combination
  group_by(guide,
           sample,
           POS,
           k562_alt_freq,
           gene,
           aa_mutation,
           Consequence_single,
           aa_position) %>%
  summarize(n_cells_guide_edit_detected = n()) %>%
  # only keep if at least 3 cells
  filter(n_cells_guide_edit_detected >= 3)

edits_detected_high_eed_plot <- sample_alt_ref_df_guides_anno %>%
  filter(sample %in% be_samples) %>%
  filter(k562_alt_freq < 0.1 | is.na(k562_alt_freq)) %>%
  filter(
    !is.na(position_on_guide) &
      position_on_guide < 20 &
      position_on_guide > -5 &
      correct_edits == TRUE
  ) %>%
  group_by(guide,
           sample,
           POS,
           k562_alt_freq,
           gene,
           aa_mutation,
           Consequence_single,
           aa_position) %>%
  filter(cell_id %in% cells_high_eed_ko_score) %>%
  summarize(n_cells_guide_edit_detected = n()) %>%
  filter(n_cells_guide_edit_detected >= 3)

# # all edits have a higher allele frequency in the cells with the guide than in the K562 sample
# edits_detected_plot <- edits_detected_plot %>%
#   mutate(allele_frequency = n_umi_sample_guide_edit / n_umi_sample_guide_coverage) %>%
#   mutate(cell_frequency = n_cells_sample_guide_edit / n_cells_sample_guide_coverage) %>%
#   mutate(k562_alt_freq = replace_na(k562_alt_freq, 0))
# 
# edits_detected_plot %>%
#   mutate(log2_allele_frequency_vs_k562 = log2((allele_frequency + 0.0001) / (k562_alt_freq + 0.0001))) %>%
#   select(sample,
#          gene,
#          aa_mutation,
#          allele_frequency,
#          cell_frequency,
#          k562_alt_freq,
#          log2_allele_frequency_vs_k562) %>%
#   ungroup() %>%
#   arrange(log2_allele_frequency_vs_k562) %>%
#   head()

####################################################################################################
# Plot kinds of edits for each editor

p1 <- edits_detected_plot %>%
  # do not count edits induced by multiple guides as independent edits
  ungroup() %>%
  select(sample, POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  group_by(sample, Consequence_single) %>%
  summarize(n_edits = n()) %>%
  ggplot(aes(x = sample, y = n_edits)) +
  geom_bar(aes(fill = Consequence_single), stat = "identity") +
  ylab("Edits") +
  scale_fill_manual(values = variant_effect_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ggtitle("Edits in cells matching guide detected in cell")

p1

p1_high_eed <- edits_detected_high_eed_plot %>%
  # do not count edits induced by multiple guides as independent edits
  ungroup() %>%
  select(sample, POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  group_by(sample, Consequence_single) %>%
  summarize(n_edits = n()) %>%
  ggplot(aes(x = sample, y = n_edits)) +
  geom_bar(aes(fill = Consequence_single), stat = "identity") +
  ylab("Edits") +
  scale_fill_manual(values = variant_effect_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ggtitle("Edits in EED KO high cells matching guide detected in cell")

p1_high_eed

####################################################################################################
# Plot distribution of guides in library across protein sequence

data_dummy <-
  data.frame(list(
    Target.Cut.Length_aa = c(0, 441, 0, 752, 0, 739),
    gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12"),
    count_bartab_18bp = c(0, 0, 0, 0, 0, 0)
  ))

p2 <- guide_metadata %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10) +
  xlab("Guide cut position on protein sequence") +
  ylab("Guides") +
  facet_wrap(~gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2)) +
  ggtitle("Histogram of guide cut positions in library")

p2

p3 <- guide_metadata %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(
    binwidth = 10,
    aes(weight = count_bartab_18bp / 10000)
  ) +
  xlab("Guide cut position on protein sequence") +
  ylab("Reads (10k)") +
  facet_wrap(~gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2)) +
  ggtitle("Histogram of guide cut positions in library, weighted by read counts")

p3

# plot guide counts in HLA+ Ctrl sample
p2_ctrl <- guide_metadata %>%
  filter(!is.na(n_cells_ctrl_hla_pos)) %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10) +
  xlab("Guide cut position on protein sequence") +
  ylab("Guides") +
  facet_wrap(~gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2)) +
  ggtitle("Histogram of guide cut positions of guides detected in Ctrl HLA+ sample")

p2_ctrl

p3_ctrl <- guide_metadata %>%
  filter(!is.na(n_cells_ctrl_hla_pos)) %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(
    binwidth = 10,
    aes(weight = n_cells_ctrl_hla_pos)
  ) +
  xlab("Guide cut position on protein sequence") +
  ylab("Cells") +
  facet_wrap(~gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2)) +
  ggtitle("Histogram of guide cut positions of guides detected in Ctrl HLA+ sample")

p3_ctrl



####################################################################################################
# Plot distribution of detected edits across base editor sample that match the guide in the cell
# sum arcoss guides and base editors

data_dummy <-
  data.frame(list(
    aa_position = c(0, 441, 0, 746, 0, 739),
    gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12")
  ))

# edits_detected_plot %>%
#   select(aa_position, aa_mutation, POS, gene, Consequence_single, guide) %>%
#   unique() %>%
#   filter(is.na(aa_position))

p4 <- edits_detected_plot %>%
  # do not count edits induced by multiple BE and guides as independent edits
  ungroup() %>%
  select(POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  ggplot(aes(x = aa_position)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10, aes(fill = Consequence_single)) +
  scale_fill_manual(values = variant_effect_colors) +
  xlab("Edit position on protein sequence") +
  ylab("Edits") +
  facet_wrap( ~ gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    legend.position = "bottom"
  ) +
  ggtitle("Histogram of edits detected across BEs \nthat match detected guide")

p4_high_eed <- edits_detected_high_eed_plot %>%
  # do not count edits induced by multiple BE and guides as independent edits
  ungroup() %>%
  select(POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  ggplot(aes(x = aa_position)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10, aes(fill = Consequence_single)) +
  scale_fill_manual(values = variant_effect_colors) +
  xlab("Edit position on protein sequence") +
  ylab("Edits") +
  facet_wrap( ~ gene, ncol = 1, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    legend.position = "bottom"
  ) +
  ggtitle("Histogram of edits detected across BEs \nthat match detected guide, cells EED KO score > 0.06")

p4
p4_high_eed

###############################################################

data_dummy <-
  data.frame(list(
    Target.Cut.Length_aa = c(0, 441, 0, 752, 0, 739),
    gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12"),
    count_bartab_18bp = c(0, 0, 0, 0, 0, 0)
  ))
p5 <- guide_metadata %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10) +
  xlab("Guide cut position on protein sequence") +
  ylab("Guides (library)") +
  facet_wrap(~gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2),
        aspect.ratio=0.4)

p6 <- guide_metadata %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(
    binwidth = 10,
    aes(weight = count_bartab_18bp / 10000)
  ) +
  xlab("Guide cut position on protein sequence") +
  ylab("Read count (10k) (library)") +
  facet_wrap(~gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2),
        aspect.ratio=0.4)


p5_ctrl <- guide_metadata %>%
  filter(!is.na(n_cells_ctrl_hla_pos)) %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10) +
  xlab("Guide cut position on protein sequence") +
  ylab("Guides (Ctrl HLA+)") +
  facet_wrap(~gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2),
        aspect.ratio=0.4)

p6_ctrl <- guide_metadata %>%
  filter(!is.na(n_cells_ctrl_hla_pos)) %>%
  mutate(Target.Cut.Length_aa = ceiling(Target.Cut.Length / 3)) %>%
  ggplot(aes(x = Target.Cut.Length_aa)) +
  geom_blank(data = data_dummy) +
  geom_histogram(
    binwidth = 10,
    aes(weight = n_cells_ctrl_hla_pos)
  ) +
  xlab("Guide cut position on protein sequence") +
  ylab("Cells (Ctrl HLA+)") +
  facet_wrap(~gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.2),
        aspect.ratio=0.4)


data_dummy <-
  data.frame(list(
    aa_position = c(0, 441, 0, 746, 0, 739),
    gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12")
  ))
p7 <- edits_detected_plot %>%
  # do not count edits induced by multiple BE and guides as independent edits
  ungroup() %>%
  select(POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  ggplot(aes(x = aa_position)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10, aes(fill = Consequence_single)) +
  scale_fill_manual(values = variant_effect_colors) +
  xlab("Edit position on protein sequence") +
  ylab("Edits") +
  facet_wrap( ~ gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    legend.position = "None",
    aspect.ratio=0.4)

p8 <- edits_detected_high_eed_plot %>%
  # do not count edits induced by multiple BE and guides as independent edits
  ungroup() %>%
  select(POS, gene, aa_mutation, Consequence_single, aa_position) %>%
  unique() %>%
  ggplot(aes(x = aa_position)) +
  geom_blank(data = data_dummy) +
  geom_histogram(binwidth = 10, aes(fill = Consequence_single)) +
  scale_fill_manual(values = variant_effect_colors) +
  xlab("Edit position on protein sequence") +
  ylab("Edits") +
  facet_wrap( ~ gene, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    legend.position = "bottom",
    aspect.ratio=0.4)

pdf(
  file.path(
    plot_dir,
    "MF03_detected_edits_and_library_protein_sequence.pdf"
  ),
  height = 5,
  width = 5
)

p1
p1_high_eed
p2
p3
p2_ctrl
p3_ctrl
p4
p4_high_eed

# log transform of bar plots is not a good idea
# https://github.com/tidyverse/ggplot2/issues/4751
# And cannot be applied in a reasonable way to a stacked bar plot
p2 + scale_y_continuous(trans = "log1p")
p3 + scale_y_continuous(trans = "log1p")
p2_ctrl + scale_y_continuous(trans = "log1p")
p3_ctrl + scale_y_continuous(trans = "log1p")

dev.off()


pdf(
  file.path(
    plot_dir,
    "MF03_detected_edits_and_library_protein_sequence_stacked.pdf"
  ),
  height = 12,
  width = 8
)
p5 / p6 / p5_ctrl / p6_ctrl / p7 / p8
# log transform of bar plots is not a good idea
# https://github.com/tidyverse/ggplot2/issues/4751
(p5 + scale_y_continuous(trans = "log1p")) / (p6 + scale_y_continuous(trans = "log1p")) / (p5_ctrl + scale_y_continuous(trans = "log1p")) / (p6_ctrl + scale_y_continuous(trans = "log1p")) / plot_spacer() / plot_spacer()

dev.off()
