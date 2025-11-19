#!/usr/bin/env Rscript

# DESCRIPTION
# Plot the average gene expression against the EED KO score.

library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

## Load data
ttest_samples_guides_filtered_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_single_guide_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

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

# check what normal gene expression in the ctrl sample is
cell_ctrl <- colnames(be_singlets)[be_singlets$sample == "Ctrl-HLA-pos"]
ctrl_expression <- be_singlets[["RNA"]]@data[c("EED", "EZH2", "SUZ12"), cell_ctrl]

ctrl_means <- rowMeans(be_singlets[["RNA"]]@data[c("EED", "EZH2", "SUZ12"), cell_ctrl])

set.seed(1)
# Plot into 2x2 grid and place legend in bottom left corner


plot_df <- ttest_samples_guides_filtered_df %>%
  dplyr::select(sample,
                guide,
                diff_eed_ko_score_normalized,
                eed_mean,
                ezh2_mean,
                suz12_mean) %>%
  mutate(target_gene = paste0(gsub("-.*", "", guide), " guide")) %>%
  pivot_longer(
    cols = c(eed_mean, ezh2_mean, suz12_mean),
    names_to = "gene_expressed",
    values_to = "mean_normalized_expression"
  ) %>%
  mutate(gene_expressed = paste0(toupper(gsub(
    "_mean", "", gene_expressed
  )), " expression")) %>%
  sample_frac(1)
  

# Plot combined base editors and CAS9
# Combine guide targets
p1 <- plot_df %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_wrap( ~ gene_expressed, scales = "free_y", axes = "all", ncol = 2) +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    legend.spacing.y = unit(0, "lines"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

legend <- cowplot::get_legend(p1)

p_nolegend <- p1 + theme(legend.position = "none")
p1 <- p_nolegend + inset_element(legend, left = 0.6, bottom = 0.2, right = 1, top = 0.2)

# Plot each guide and each guide target
p2 <- plot_df %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_grid(target_gene ~ gene_expressed,
             scales = "free_y",
             axes = "all") +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  )


# Plot only base editors
p3 <- plot_df %>%
  filter(sample != "CAS9-HLA-pos") %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_wrap( ~ gene_expressed, scales = "free_y", axes = "all", ncol = 2) +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    legend.spacing.y = unit(0, "lines"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

legend <- cowplot::get_legend(p3)

p_nolegend <- p3 + theme(legend.position = "none")
p3 <- p_nolegend + inset_element(legend, left = 0.6, bottom = 0.2, right = 1, top = 0.2)


p4 <- plot_df %>%
  filter(sample != "CAS9-HLA-pos") %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_grid(target_gene ~ gene_expressed,
             scales = "free_y",
             axes = "all") +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  )


# Plot only CAS9
p5 <- plot_df %>%
  filter(sample == "CAS9-HLA-pos") %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_wrap( ~ gene_expressed, scales = "free_y", axes = "all", ncol = 2) +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    legend.spacing.y = unit(0, "lines"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

legend <- cowplot::get_legend(p5)

p_nolegend <- p5 + theme(legend.position = "none")
p5 <- p_nolegend + inset_element(legend, left = 0.6, bottom = 0.2, right = 1, top = 0.2)


p6 <- plot_df %>%
  filter(sample == "CAS9-HLA-pos") %>%
  ggplot(aes(x = diff_eed_ko_score_normalized, y = mean_normalized_expression, color = sample)) +
  geom_vline(xintercept = 0.06, linetype = 3, alpha = 0.5) +
  geom_hline(
    data = data.frame(
      gene_expressed = paste0(names(ctrl_means), " expression"),
      yintercept = as.numeric(ctrl_means)
    ),
    aes(yintercept = yintercept, linetype = "Average Ctrl HLA+"),
  ) +
  geom_point(aes(color = sample), alpha = 0.5, size = 1) +
  scale_color_manual(values = samples_color) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_grid(target_gene ~ gene_expressed,
             scales = "free_y",
             axes = "all") +
  ylab("Average normalized expression") +
  xlab("EED KO score") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  )



pdf(
  file.path(plot_dir, "MF03_guide_target_expression.pdf"),
  height = 6,
  width = 9
)

p1
p2
p3
p4
p5
p6

dev.off()
