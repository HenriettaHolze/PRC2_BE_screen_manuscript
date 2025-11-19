#!/usr/bin/env Rscript

# DESCRIPTION
# UMAP colored by EED KO score (binary and continuous)

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

cells_high_eed_ko <- colnames(be_singlets)[be_singlets$eed_ko_score_normalized > 0.06]
p1 <- DimPlot(
  be_singlets,
  cells.highlight = list("EED KO score > 0.06" = cells_high_eed_ko),
  sizes.highlight = 0.2,
  cols.highlight = c("#a90505")
) +
  coord_equal() +
  NoAxes()

cells_low_eed_ko <- colnames(be_singlets)[be_singlets$eed_ko_score_normalized < 0.06]
p2 <- DimPlot(
  be_singlets,
  cells.highlight = list("EED KO score < 0.06" = cells_low_eed_ko),
  sizes.highlight = 0.2,
  cols.highlight = c("#a90505")
) +
  coord_equal() +
  NoAxes()

p3 <- FeaturePlot(be_singlets, features = "eed_ko_score_normalized", order = TRUE) +
  scale_color_viridis_c("EED KO score", direction = -1, option = "magma") +
  coord_equal() +
  NoAxes() +
  ggtitle(element_blank())

pdf(
  file.path(plot_dir, "MF03_umap_eed_ko_cutoff.pdf"),
  height = 4.5,
  width = 12
)

p1 | p2
p3

dev.off()
