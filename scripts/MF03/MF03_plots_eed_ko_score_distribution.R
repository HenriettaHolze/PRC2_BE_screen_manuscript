#!/usr/bin/env Rscript

# DESCRIPTION
# Plot the distribution of the EED KO score and HLA-B per sample as a violin plot

library(tidyverse)
library(Seurat)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "") dir.create(plot_dir, recursive = TRUE)

## Load data
be_singlets <-
  qs::qread(paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
  ))

## color scheme
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


## Plot EED KO score and HLA-B distribution per sample
samples_plot_order <- be_singlets@meta.data %>%
  group_by(sample) %>%
  summarize(median(eed_ko_score_normalized)) %>%
  arrange(`median(eed_ko_score_normalized)`) %>%
  pull(sample)

Idents(object = be_singlets) <- "sample"
be_singlets$sample_factor <- factor(x = be_singlets$sample, levels = samples_plot_order)

p1 <- VlnPlot(
  be_singlets,
  features = "eed_ko_score_normalized",
  cols = samples_color,
  group.by = "sample_factor",
  pt.size = 0
) + 
  NoLegend() + 
  xlab(element_blank()) +
  ylab("EED KO score") +
  geom_hline(yintercept = 0.06, linetype = 2)

p2 <- VlnPlot(
  be_singlets,
  features = "HLA-B",
  cols = samples_color,
  group.by = "sample_factor",
  assay = "RNA",
  slot = "data",
  pt.size = 0
) + 
  NoLegend() + 
  xlab(element_blank())

# Without V14
p3 <- VlnPlot(
  subset(be_singlets, sample != "V14-HLA-pos"),
  features = "eed_ko_score_normalized",
  cols = samples_color,
  group.by = "sample_factor",
  pt.size = 0
) + 
  NoLegend() + 
  xlab(element_blank()) +
  ylab("EED KO score") +
  geom_hline(yintercept = 0.06, linetype = 2)

p4 <- VlnPlot(
  subset(be_singlets, sample != "V14-HLA-pos"),
  features = "HLA-B",
  cols = samples_color,
  group.by = "sample_factor",
  assay = "RNA",
  slot = "data",
  pt.size = 0
) + 
  NoLegend() + 
  xlab(element_blank())

# plot target gene expression across samples
gene_plots <- lapply(c("EZH2", "EED", "SUZ12"), function(gene) {
  VlnPlot(
    be_singlets,
    features = gene,
    cols = samples_color,
    group.by = "sample_factor",
    assay = "RNA",
    slot = "data",
    pt.size = 0
  ) +
    NoLegend() +
    xlab(element_blank())
})

# check if in BE samples with high EED KO, the target genes are downregulated.
be_samples <- c("V14", "V7", "TADA", "RAPO", "CDA1")

selected_cells <- be_singlets@meta.data %>%
  # keep BE sample cells with high EED KO and Ctrl sample cells
  filter(eed_ko_score_normalized > 0.06 | sample == "Ctrl-HLA-pos") %>%
  filter(!sample %in% c("CAS9-HLA-pos", "K562")) %>%
  rownames()
be_singlets_high_eed <- subset(be_singlets, cells = selected_cells)

# plot target gene expression across samples
gene_plots_high_eed <- lapply(c("EZH2", "EED", "SUZ12"), function(gene) {
  VlnPlot(
    be_singlets_high_eed,
    features = gene,
    cols = samples_color,
    group.by = "sample_factor",
    assay = "RNA",
    slot = "data",
    pt.size = 0
  ) +
    NoLegend() +
    xlab(element_blank())
})


deg_eed_high_target_genes <- lapply(be_samples, function(selected_sample) {
  df <- FindMarkers(
    be_singlets_high_eed,
    features = c("EZH2", "EED", "SUZ12"),
    ident.2 = "Ctrl-HLA-pos",
    ident.1 = paste0(selected_sample, "-HLA-pos"),
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1
  )
  df$gene <- rownames(df)
  return(df)
})
names(deg_eed_high_target_genes) <- paste0(be_samples, "-HLA-pos")

deg_eed_high_target_genes <- bind_rows(deg_eed_high_target_genes, .id = "sample")

# the adjusted p-value is not meaningful because the function adjusts it to not only the number of tested genes
deg_eed_high_target_genes$p_val_adj <- deg_eed_high_target_genes$p_val * nrow(deg_eed_high_target_genes)

p5 <- deg_eed_high_target_genes %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = sample)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -0.5, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 0.5, linetype = 2, alpha = 0.5) +
  facet_wrap(~gene, axes = "all", ncol = 2) +
  xlim(c(-1, 1)) +
  scale_color_manual(values = samples_color) +
  guides(color=guide_legend(nrow=2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

be_singlets_high_eed@meta.data %>%
  group_by(sample) %>%
  summarize(mean(nFeature_RNA))


pdf(
  file.path(
    plot_dir, 
    "MF03_eed_ko_score_distribution.pdf"
  ), 
  height = 4.5, 
  width = 3.5
)

p1
p2
p3
p4
gene_plots[[1]]
gene_plots[[2]]
gene_plots[[3]]
gene_plots_high_eed[[1]]
gene_plots_high_eed[[2]]
gene_plots_high_eed[[3]]
p5

dev.off()
