#!/usr/bin/env Rscript

# DESCRIPTION
# UMAP plots with samples highlighted

library(qs)
library(Seurat)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )

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
samples <-
  c(
    "K562",
    "Ctrl-HLA-pos",
    "CAS9-HLA-pos",
    "V7-HLA-pos",
    "TADA-HLA-pos",
    "CDA1-HLA-pos",
    "V14-HLA-pos",
    "RAPO-HLA-pos"
  )
names(samples_color) <- samples

p_list <- list()
for (sample in samples) {
  cells_highlight = colnames(be_singlets)[be_singlets$sample == sample]
  p_list[[sample]] <-
    DimPlot(be_singlets,
            cells.highlight = cells_highlight,
            cols.highlight = samples_color[sample],
            sizes.highlight = 0.5) +
    coord_equal() +
    ggtitle(sample) +
    theme_void() +
    NoLegend() +
    NoAxes() +
    theme(plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ))
}

ggpubr::ggarrange(plotlist = p_list,
                  ncol = 3,
                  nrow = 3)



pdf(file.path(plot_dir, "MF03_umap_samples.pdf"),
    height = 12,
    width = 12)

ggpubr::ggarrange(plotlist = p_list,
                  ncol = 3,
                  nrow = 3)

dev.off()

png(
  file.path(plot_dir, "MF03_umap_samples.png"),
  width = 24,
  height = 24,
  units = "cm",
  res = 300
)

ggpubr::ggarrange(plotlist = p_list,
                  ncol = 3,
                  nrow = 3)

dev.off()
