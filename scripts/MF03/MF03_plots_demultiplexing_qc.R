#!/usr/bin/env Rscript

## DESCRIPTION
## QC plot of sample demultiplexing: number of cells per sample and doublets, negatives

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
dir.create(plot_dir, recursive = TRUE)


be <- qs::qread(paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed.qs"))

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

samples_color["Doublet"] <- "black"
samples_color["Negative"] <- "gray"

p1 <- be@meta.data %>%
  mutate(is_sample = sample %in% samples) %>%
  mutate(sample = factor(sample, levels = names(samples_color), ordered = T)) %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  facet_grid(~ is_sample, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = samples_color) +
  ylab("# cells") +
  xlab("") +
  ggtitle("Cells per sample after LMO demultiplexing") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    strip.text = element_blank(),
    legend.position = "None",
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.spacing.x = unit(1, "lines")
  )

pdf(
  file.path(plot_dir, "MF03_demultiplexing_qc.pdf"),
  height = 4,
  width = 4
)

p1

dev.off()
