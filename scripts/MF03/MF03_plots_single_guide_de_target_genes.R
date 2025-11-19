#!/usr/bin/env Rscript

# DESCRIPTION
# Plot differential expression of target genes for each guide against the Ctrl HLA+ sample
# across all samples.
# Only consider cells with a single guide integrated.

# INPUT
# Differential expression results of target genes (wilcoxon rank sum test, using Seurat's FindMarkers function)

library(Seurat)
library(tidyverse)
library(ltm)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)


## Load data
# wilcoxon rank sum test of normalized gene expression
findmarkers_res_samples_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_target_gene_de_single_guide_vs_ctrl_wilcoxon.tsv"
  )
)


# Define samples of interest
samples <- c(
  "V7-HLA-pos",
  "V14-HLA-pos",
  "TADA-HLA-pos",
  "CDA1-HLA-pos",
  "RAPO-HLA-pos",
  "CAS9-HLA-pos"
)

# genes of interest
target_genes <- c("EED", "EZH2", "SUZ12")

x_min_max <- c(min(findmarkers_res_samples_df$avg_log2FC), max(findmarkers_res_samples_df$avg_log2FC))
y_min_max <- c(min(-log10(findmarkers_res_samples_df$p_val)), max(-log10(findmarkers_res_samples_df$p_val)))

plot_list <- lapply(samples, function(selected_sample) {
  findmarkers_res_df <- findmarkers_res_samples_df %>% filter(sample == selected_sample)
  
  n_test <- length(unique(findmarkers_res_df$guide))
  
  # identify guides that induce reduced expression of their target
  # i.e. NMD guides or knock-down guides
  nmd_guides <- findmarkers_res_df %>%
    mutate(target = gsub("_.*", "", guide)) %>%
    filter(target == gene) %>%
    filter(p_val < 0.05, avg_log2FC < 0) %>%
    pull(guide)
  
  findmarkers_res_df <- findmarkers_res_df %>%
    mutate(label = ifelse(p_val < 0.05 & avg_log2FC < 0, guide, ""))
  
  findmarkers_res_df <- findmarkers_res_df %>%
    mutate(nmd_guide = guide %in% nmd_guides) %>%
    # where target and gene are the same, to highlight rectangle in those facets
    mutate(matches = gsub("_.*", "", guide) == gene) %>%
    mutate(gene = paste0(gene, " expression")) %>%
    arrange(desc(matches)) %>%
    group_by(gene) %>%
    mutate(rank = row_number()) %>%
    mutate(target = paste0(gsub("_.*", "", guide), " guide"))
  
  p1 <- findmarkers_res_df %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val))) +
    # highlight how NMD guides were selected
    geom_rect(
      aes(
        xmin = -Inf * matches * (rank == 1),
        xmax = 0,
        ymin = -log10(0.05),
        ymax = Inf
      ),
      fill = "#005AB5",
      alpha = 0.1
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      alpha = 0.5
    ) +
    geom_vline(xintercept = 0,
               linetype = 2,
               alpha = 0.5) +
    # highlight NMD guides across expression of all genes
    geom_point(
      aes(color = nmd_guide),
      alpha = 0.7,
      size = 1.5,
      stroke = 0
    ) +
    ggrepel::geom_text_repel(
      aes(label = label,
          color = nmd_guide),
      min.segment.length = 0,
      max.overlaps = Inf,
      size = 3
    ) +
    scale_color_manual(
      "Guide induces \nreduced expression \nof target",
      values = c("TRUE" = "#005AB5", "FALSE" = "#DC3220")
    ) +
    facet_grid(target ~ gene) +
    ylim(y_min_max) +
    ggtitle(paste0(selected_sample, " guide target expression vs Ctrl HLA+")) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p1)
})


pdf(
  file.path(plot_dir, "MF03_single_guide_de_target_genes.pdf"),
  height = 7,
  width = 9
)

plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]
plot_list[[5]]
plot_list[[6]]

dev.off()
