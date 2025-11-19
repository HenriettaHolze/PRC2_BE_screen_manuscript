#!/usr/bin/env Rscript

# DESCRIPTION
# Plot EED KO score vs the number of expressed genes per cell, to show that PRC2 KO induces increased expression

# INPUT
# Seurat object with EED KO score and metadata per cell

library(Seurat)
library(tidyverse)
library(qs)
library(ggpmisc)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)


## Load data

be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )

## Color scheme
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

p1 <- be_singlets@meta.data %>%
  ggplot(aes(x = nFeature_RNA, y = eed_ko_score_normalized, color = sample)) +
  geom_hline(yintercept = 0.06,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.1,
             stroke = 0,
             size = 1) +
  geom_smooth(
    formula = 'y ~ x',
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    color = "black",
    formula = y ~ x,
    parse = TRUE,
    size = 3
  ) +
  facet_wrap( ~ sample, axes = "all") +
  scale_color_manual(values = samples_color) +
  ggtitle("Number of expressed genes vs. EED KO score") +
  ylab("EED KO score") +
  xlab("Detected genes per cell") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  )

p2 <- be_singlets@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = eed_ko_score_normalized, color = sample)) +
  geom_hline(yintercept = 0.06,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.1,
             stroke = 0,
             size = 1) +
  geom_smooth(
    formula = 'y ~ x',
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    color = "black",
    formula = y ~ x,
    parse = TRUE,
    size = 3
  ) +
  facet_wrap( ~ sample, axes = "all") +
  scale_color_manual(values = samples_color) +
  ggtitle("Number of detected transcripts vs. EED KO score") +
  ylab("EED KO score") +
  xlab("Detected transcripts per cell") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  )


pdf(
  file.path(plot_dir, "MF03_eed_ko_score_vs_genes_detected.pdf"),
  height = 7.5,
  width = 7
)

p1
p2

dev.off()
