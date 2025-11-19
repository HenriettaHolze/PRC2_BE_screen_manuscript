#!/usr/bin/env Rscript

## DESCRIPTION
## - plot number of cells per sample with each target gene detected in ONT and 10X data

library(Seurat)
library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "") dir.create(plot_dir, recursive = TRUE)


## Load data
cellsnp_res_dir <-
  "/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sc_nanopore_variant_calling/"
vcf <-
  read.table(
    paste0(
      cellsnp_res_dir,
      'sicelore_cellsnp_1a_SUZ12_EED_EZH2/cellSNP.base.vcf.gz'
    ),
    comment.char = "",
    skip = 1,
    header = T
  )
# read in reference and alternative UMI counts for cells that were demultiplexed into samples
mat_alt_samples <-
  Matrix::readMM(paste0(
    result_dir,
    "MF03/nanopore/mat_alt_samples_sicelore_cellsnp.mtx"
  ))
mat_ref_alt_samples <-
  Matrix::readMM(paste0(
    result_dir,
    "MF03/nanopore/mat_ref_alt_samples_sicelore_cellsnp.mtx"
  ))

cells <-
  read_lines(paste0(
    result_dir,
    "MF03/nanopore/cells_samples_sicelore_cellsnp.tsv"
  ))
positions <-
  read_lines(paste0(
    result_dir,
    "MF03/nanopore/positions_samples_sicelore_cellsnp.tsv"
  ))
rownames(mat_alt_samples) <- positions
rownames(mat_ref_alt_samples) <- positions
colnames(mat_alt_samples) <- cells
colnames(mat_ref_alt_samples) <- cells
cells <- NULL
positions <- NULL

# Seurat object with annotation which cell belongs to which sample
be_singlets <-
  qs::qread(
    paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs")
  )
be_singlets$cellid <- sub("P[12]_(.*)-[12]", "\\1", colnames(be_singlets))

be_singlets_metadata <- be_singlets@meta.data
n_cells_sample <- as.list(table(be_singlets_metadata$sample))

vcf <- vcf %>%
  mutate(gene = ifelse(X.CHROM == 11, "EED", ifelse(
    X.CHROM == 7, "EZH2", ifelse(X.CHROM == 17, "SUZ12", NA)
  )))

## count cells with transcript

samples <- as.list(unique(be_singlets_metadata$sample))
names(samples) <- samples

genes <- list("EED", "SUZ12", "EZH2")
names(genes) <- genes

# iterate over genes
ont_cells <- lapply(genes, function(selected_gene) {
  # iterate over samples
  lapply(samples, function(selected_sample) {
    cat(selected_gene, selected_sample, "\n")
    positions <- vcf %>% filter(gene == selected_gene) %>% pull(POS)
    positions <- rownames(mat_ref_alt_samples)[rownames(mat_ref_alt_samples) %in% positions]
    cells <- be_singlets_metadata %>%
      filter(sample == selected_sample,
             cellid %in% colnames(mat_ref_alt_samples)) %>%
      pull(cellid)
    
    # count number of cells with number of transcripts detected
    max_transcripts_per_cell <- apply(mat_ref_alt_samples[positions, cells], 2, max)
    n_transcripts_cell <- table(max_transcripts_per_cell)
  
    return(n_transcripts_cell)
  })
})

# make df with number of cells with number of transcripts detected
ont_n_transcripts_cell <- lapply(genes, function(selected_gene) {
  lapply(samples, function(selected_sample) {
    # tmp <- as.data.frame(ont_cells[[selected_gene]][[selected_sample]]$n_transcripts_cell)
    tmp <- as.data.frame(ont_cells[[selected_gene]][[selected_sample]])
    tmp$gene <- selected_gene
    tmp$sample <- selected_sample
    return(tmp)
  })
})

ont_n_transcripts_cell_df <- bind_rows(lapply(ont_n_transcripts_cell, function(x) bind_rows(x)))

# the above takes some time to run, so I'll save the intermediate file
write_tsv(
  ont_n_transcripts_cell_df,
  paste0(
    result_dir,
    "MF03/nanopore/MF03_ont_transcript_detection.tsv"
  )
)

# the above takes some time to run, so I'll save the intermediate file
ont_n_transcripts_cell_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_ont_transcript_detection.tsv"
  )
)


# plot

p1 <- ont_n_transcripts_cell_df %>%
  filter(max_transcripts_per_cell != 0) %>%
  mutate(max_transcripts_per_cell = as.integer(as.character(max_transcripts_per_cell))) %>%
  ggplot(aes(x = gene, y = Freq)) +
  ggplot2::geom_hline(
    data = data.frame(
      sample = names(n_cells_sample),
      yintercept = as.numeric(n_cells_sample)
    ),
    aes(yintercept = yintercept, linetype = "# cells in sample"),
  ) +
  geom_bar(aes(fill = max_transcripts_per_cell), stat = "identity") +
  facet_wrap(~sample, axes = "all", scales = "free_y") +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  ylab("# cells") +
  scale_fill_viridis_c("ONT transcripts \nper cell", option = "B", direction = -1, end = 0.95) +
  ggtitle("ONT transcripts detected per cell") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# get number of cells with ONT transcript detected per sample and gene
ont_cells_df <- ont_n_transcripts_cell_df %>%
  filter(max_transcripts_per_cell != 0) %>%
  group_by(sample, gene) %>%
  summarize(n_cells = sum(Freq))

# add the correct number of cells without transcript for the respective gene in the samle
tmp <- t(bind_rows(n_cells_sample, .id = "sample")) %>%
  as.data.frame() %>%
  rownames_to_column("sample")
colnames(tmp) <- c("sample", "total_cells")

ont_n_transcripts_cell_df <- rbind(
  ont_n_transcripts_cell_df %>% filter(max_transcripts_per_cell != 0),
  merge(ont_cells_df, tmp, by = c("sample")) %>%
    mutate(Freq = total_cells - n_cells) %>%
    mutate(max_transcripts_per_cell = 0) %>%
    select(sample, gene, Freq, max_transcripts_per_cell)
)

# plot number of transcripts per cells as a violin plot
p2 <- ont_n_transcripts_cell_df %>%
  uncount(Freq) %>%
  ggplot(aes(x = gene, y = max_transcripts_per_cell)) +
  geom_violin(scale = "width") +
  facet_wrap(~sample, axes = "all") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  ) +
  ylab("Maximum transcript coverage per cell")

# plot number of cells with transcript detected
p3 <- ont_cells_df %>%
  ggplot(aes(x = gene, y = n_cells)) +
  geom_bar(stat = "identity") +
  facet_wrap(~sample, axes = "all", scales = "free_y") +
  ggplot2::geom_hline(
    data = data.frame(
      sample = names(n_cells_sample),
      yintercept = as.numeric(n_cells_sample)
    ),
    aes(yintercept = yintercept, linetype = "# cells in sample"),
  ) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  ylab("# cells") +
  ggtitle("Cells with ONT transcript detected") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )


# get number of cells with 10X transcript detected
expression_10x <- be_singlets@assays$RNA@counts[c("EED", "EZH2", "SUZ12"),]
# make pretty dataframe
tenx_cells_df <- merge(be_singlets_metadata, t(as.data.frame(expression_10x != 0)), by="row.names") %>%
  dplyr::select(sample, EED, SUZ12, EZH2) %>%
  pivot_longer(cols = -sample, names_to = "gene", values_to = "gene_detected") %>%
  group_by(sample, gene) %>%
  summarize(n_cells = sum(gene_detected))

p4 <- tenx_cells_df %>%
  ggplot(aes(x = gene, y = n_cells)) +
  geom_bar(stat = "identity") +
  facet_wrap(~sample, axes = "all", scales = "free_y") +
  ggplot2::geom_hline(
    data = data.frame(
      sample = names(n_cells_sample),
      yintercept = as.numeric(n_cells_sample)
    ),
    aes(yintercept = yintercept, linetype = "# cells in sample"),
  ) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("# cells") +
  ggtitle("Cells with 10X transcript detected") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# combine 10X and ONT transcript detection in 1 figure
ont_cells_df$method <- "ONT"
tenx_cells_df$method <- "10X"

p5 <- rbind(tenx_cells_df, ont_cells_df) %>%
  ggplot(aes(x = gene, y = n_cells, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~sample, axes = "all", scales = "free_y") +
  ggplot2::geom_hline(
    data = data.frame(
      sample = names(n_cells_sample),
      yintercept = as.numeric(n_cells_sample)
    ),
    aes(yintercept = yintercept, linetype = "# cells in sample"),
  ) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("# cells") +
  ggtitle("Cells with transcript detected") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )



pdf(
  file.path(
    plot_dir, 
    "MF03_transcript_detection.pdf"
  ), 
  height = 6, 
  width = 7
)

p1
p2
p3
p4
p5

dev.off()
