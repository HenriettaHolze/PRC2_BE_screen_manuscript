#!/usr/bin/env Rscript

# DESCRIPTION

library(tidyverse)
library(ggpubr)
library(seqinr)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

####################################################################################################
# Plot the results of the heterozygosity analysis

predicted_functional_variants_hetero <- read_tsv(
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

p1 <- predicted_functional_variants_hetero %>%
  ggplot(
    aes(
      x = cells_alt_and_ref_detected / cells_alt_detected,
      y = cells_alt_detected,
      color = Consequence_single
    )
  ) +
  geom_point() +
  facet_grid(gene ~ sample, axes = "all") +
  xlim(c(0, 1)) +
  ylim(c(0, NA)) +
  ylab("cells with alt detected \n(~enrichment in screen)") +
  xlab("heterozygous cells vs all cells with alt detected") +
  ggrepel::geom_text_repel(
    aes(label = aa_mutation_1letter),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  ) +
  scale_color_manual(values = variant_effect_colors) +
  custom_theme

p2 <- predicted_functional_variants_hetero %>%
  ggplot(aes(
    x = umi_alt / (umi_alt + umi_ref),
    y = cells_alt_detected,
    color = Consequence_single
  )) +
  geom_point() +
  facet_grid(gene ~ sample, axes = "all") +
  xlim(c(0, 1)) +
  ylim(c(0, NA)) +
  ylab("cells with alt detected \n(~enrichment in screen)") +
  xlab("alt UMI vs total UMI in cells with alt detected") +
  ggrepel::geom_text_repel(
    aes(label = aa_mutation_1letter),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F
  ) +
  scale_color_manual(values = variant_effect_colors) +
  custom_theme


# collate the 2 above figures into 1.
# Extra function to get ggrepel working with size aesthetic
# https://github.com/slowkow/ggrepel/issues/14#issuecomment-624823565
my_pal <- function(range = c(1, 6)) {
  force(range)
  function(x)
    scales::rescale(x, to = range, from = c(0, 1))
}
p3 <- predicted_functional_variants_hetero %>%
  ggplot(
    aes(
      x = umi_alt / (umi_alt + umi_ref),
      y = cells_alt_and_ref_detected / cells_alt_detected,
      fill = Consequence_single,
      color = Consequence_single,
      size = cells_alt_detected
    )
  ) +
  geom_point(alpha = 0.7,
             pch = 21,
             stroke = NA) +
  facet_grid(gene ~ sample, axes = "all") +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  xlab("alt UMI vs total UMI in cells with alt detected") +
  ylab("heterozygous cells vs all \ncells with alt detected") +
  ggrepel::geom_text_repel(
    aes(label = aa_mutation_1letter, point.size = cells_alt_detected),
    size = 3,
    point.padding = 0,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = F,
    box.padding = 0.3
  ) +
  scale_color_manual(values = variant_effect_colors) +
  scale_fill_manual(values = variant_effect_colors) +
  guides(size = guide_legend(title = "cells with alt detected \n(~enrichment in screen)")) +
  continuous_scale(
    aesthetics = c("size", "point.size"),
    scale_name = "size",
    palette = my_pal(c(2, 8)),
    guide = guide_legend(override.aes = list(label = ""))
  ) +
  custom_theme +
  coord_fixed()

pdf(
  file.path(plot_dir, "MF03_heterozygosity.pdf"),
  height = 9,
  width = 8
)

(p1 + theme(legend.position = "None")) / (p2 + theme(legend.position = "None")) + p3

dev.off()

# # Plot cells with specific mutation on UMAP
# selected_sample <- "TADA"
# pos <- predicted_functional_variants_hetero %>%
#   filter(aa_mutation_1letter == "C560R") %>% pull(POS) %>% as.character()
#
# cells_sample <- be_singlets_metadata %>%
#   filter(grepl(selected_sample, sample),
#          eed_ko_score_normalized > 0.06) %>%
#   pull(cellid)
# cells_sample <- cells_sample[cells_sample %in% colnames(mat_alt_samples)]
# cells_mutation <- colnames(mat_alt_samples)[mat_alt_samples[pos, ] != 0]
# cells_mutation <- cells_mutation[cells_mutation %in% cells_sample]
#
# cell_ids_mutation <- be_singlets_metadata %>%
#   filter(cellid %in% cells_mutation) %>% rownames()
#
# DimPlot(be_singlets, cells.highlight = cell_ids_mutation)
