#!/usr/bin/env Rscript

# DESCRIPTION
# Figures on guide detection in single-cell data.
# Purpose: give an idea for what fraction of cells we were able to recover guides,
# how we chose the minimum UMI threshold,
# and how bad the multiple guide integration problem is

# INPUT
# - Seurat object with demultiplexed samples (to get number of cells per sample with scRNA-seq data)
# - Guide counts with samples annotated

library(Seurat)
library(tidyverse)
library(qs)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
dir.create(plot_dir, recursive = TRUE)

## Read data

# Demultiplexed cells
be_singlets <-
  qs::qread(paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap.qs"
  ))
n_cells_sample <- as.list(table(be_singlets$sample))

counts <- read_tsv(paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts.tsv"))

# for plotting, unique cell ID column is expected to be cellid
counts <- counts %>% mutate(cellid = cell_id)


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


### Number of UMI supporting most abundant guide per cell

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


plotUmiPerBarcodeSample <- function (counts,
                                     fraction = TRUE,
                                     n_cells_sample = NULL)
{
  max.umi.per.cell <- counts %>%
    dplyr::group_by(.data$cellid, .data$sample) %>%
    dplyr::summarise(max = max(.data$bc.umi.count)) %>%
    group_by(sample) %>%
    mutate(n_sample = n())
  max.umi.per.cell <-
    max.umi.per.cell %>%
    dplyr::group_by(.data$sample, .data$n_sample, .data$max) %>%
    dplyr::count() %>%
    dplyr::mutate(frac = .data$n / .data$n_sample)
  if (fraction) {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = .data$max, y = .data$frac))
  }
  else {
    p <- max.umi.per.cell %>%
      ggplot2::ggplot(aes(x = .data$max, y = .data$n))
    if (!is.null(n_cells_sample)) {
      p <- p +
        ggplot2::geom_hline(
          data = data.frame(
            sample = names(n_cells_sample),
            yintercept = as.numeric(n_cells_sample)
          ),
          aes(yintercept = yintercept, linetype = "# cells in sample"),
        ) +
        ggplot2::scale_linetype_manual(name = "", values = c(2))
    }
  }
  p <- p + ggplot2::theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
    ggplot2::geom_bar(aes(fill = sample), stat = "identity", width = 1, show.legend = F) +
    ggplot2::xlab("# UMI") +
    ggplot2::ggtitle("Number of UMI supporting the most frequent guide per cell") +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::theme_bw() +
    facet_wrap( ~ .data$sample, scales = "free_y", axes = "all", ncol = 4) +
    ggplot2::theme(panel.grid = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(), 
                   legend.position = "None")
  if (fraction) {
    p <- p + ylab("Fraction of cells")
  }
  else {
    p <- p + ylab("# cells")
  }
  return(p)
}

p1 <- plotUmiPerBarcodeSample(counts, fraction = F, n_cells_sample = n_cells_sample) + 
  xlim(0, 75) +
  scale_fill_manual(values = samples_color)

p2 <- plotUmiPerBarcodeSample(counts, fraction = T) + 
  xlim(0, 75) +
  scale_fill_manual(values = samples_color)


### Number of cells with single guide detected for different UMI thresholds

getCellsSingleBarcode <- function(counts, umi_threshold = 1) {
  counts %>%
    filter(bc.umi.count >= umi_threshold) %>%
    group_by(cellid) %>%
    summarize(n_cells = n()) %>%
    filter(n_cells == 1) %>%
    pull(cellid)
}

plotCellsSingleBarcode <- function(counts, n_cells_sample = NULL) {
  samples <- as.list(unique(counts$sample))
  names(samples) <- samples
  df <- lapply(samples, function(selected_sample) {
    counts_sample <- counts %>% filter(sample == selected_sample)
    threshold_cells <- as.list(sapply(seq(1, 10), function(i)
      length(
        getCellsSingleBarcode(counts = counts_sample, umi_threshold = i)
      )))
    names(threshold_cells) <- seq(1, 10)
    return(threshold_cells)
  })
  
  df_long <- bind_rows(df, .id = "sample") %>%
    pivot_longer(cols = -sample,
                 names_to = "UMI threshold",
                 values_to = "n_cells") %>%
    mutate(`UMI threshold` = as.integer(`UMI threshold`))
  
  p <- df_long %>%
    ggplot(aes(x = `UMI threshold`, y = n_cells)) +
    geom_bar(aes(fill = sample), stat = "identity", width = 0.75, show.legend = F) +
    facet_wrap( ~ sample, axes = "all", ncol = 4) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::ggtitle(paste(
      "Number of cells with single guide annotated at different UMI thresholds"
    )) +
    ggplot2::theme(panel.grid = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(), 
                   legend.position = "None") +
    ylab("# cells")
  
  if (!is.null(n_cells_sample)) {
    p <- p +
      geom_hline(
        data = data.frame(
          sample = names(n_cells_sample),
          yintercept = as.numeric(n_cells_sample)
        ),
        aes(yintercept = yintercept, linetype = "# cells in sample"),
      ) +
      scale_linetype_manual(name = "", values = c(2))
  }
  
  return(p)
}

p3 <- plotCellsSingleBarcode(counts, n_cells_sample) +
  scale_fill_manual(values = samples_color)

### Number of guides detected per cell

plotBarcodesPerCellSample <-
  function (counts,
            fraction = TRUE,
            aggregated = FALSE,
            notDetected = "",
            sep = ";",
            n_cells_sample = NULL) {
    if (class(counts)[1] == "Seurat") {
      counts <- counts@meta.data
      aggregated <- TRUE
    }
    else if (class(counts)[1] == "SingleCellExperiment") {
      counts <- as.data.frame(counts@colData)
      aggregated <- TRUE
    }
    if (aggregated) {
      lineagePerCell.dist.df <- counts %>%
        tibble::rownames_to_column("cellid") %>%
        dplyr::select(.data$cellid, .data$barcode) %>%
        dplyr::mutate(barcode = ifelse(.data$barcode ==
                                         notDetected, NA, .data$barcode)) %>%
        dplyr::mutate(number_of_lineage_barcodes = stringr::str_count(.data$barcode, pattern = sep) + 1) %>%
        dplyr::mutate(number_of_lineage_barcodes = ifelse(
          is.na(.data$number_of_lineage_barcodes),
          0,
          .data$number_of_lineage_barcodes
        )) %>%
        dplyr::select(-.data$barcode)
    }
    else {
      lineagePerCell.dist.df <- counts %>%
        dplyr::select(.data$cellid, .data$barcode, .data$sample) %>%
        dplyr::group_by(.data$cellid, .data$sample) %>%
        dplyr::tally(., name = "number_of_lineage_barcodes")
    }
    lineagePerCell.dist.df <- lineagePerCell.dist.df %>%
      # calculate the fraction of all cells in sample
      merge(data.frame(list(
        n_sample = unlist(n_cells_sample),
        sample = names(n_cells_sample)
      )), by = "sample") %>%
      group_by(sample, n_sample) %>%
      dplyr::count(.data$number_of_lineage_barcodes) %>%
      dplyr::mutate(frac = .data$n / n_sample)
    
    if (fraction) {
      p <- lineagePerCell.dist.df %>%
        ggplot2::ggplot(aes(
          x = .data$number_of_lineage_barcodes,
          y = .data$frac
        ))
    }
    else {
      p <- lineagePerCell.dist.df %>%
        ggplot2::ggplot(aes(x = .data$number_of_lineage_barcodes, y = .data$n))
      
      if (!is.null(n_cells_sample)) {
        p <- p +
          geom_hline(
            data = data.frame(
              sample = names(n_cells_sample),
              yintercept = as.numeric(n_cells_sample)
            ),
            aes(yintercept = yintercept, linetype = "# cells in sample"),
          ) +
          scale_linetype_manual(name = "", values = c(2))
      }
      
    }
    p <- p + ggplot2::theme(
      axis.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold")
    ) + ggplot2::geom_bar(aes(fill = sample), stat = "identity", width = 0.75, show.legend = F) +
      ggplot2::xlab("# guides") +
      ggplot2::ggtitle(paste("Number of guides detected per cell")) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = integer_breaks()) +
      facet_wrap( ~ .data$sample, axes = "all", ncol = 4) +
      ggplot2::theme(panel.grid = element_blank(),
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(), 
                     legend.position = "None")
    
    if (fraction) {
      p <- p + ylab("Fraction of cells in sample")
    }
    else {
      p <- p + ylab("# cells")
    }
    return(p)
  }

p4 <- plotBarcodesPerCellSample(counts, fraction = F, n_cells_sample = n_cells_sample) +
  scale_fill_manual(values = samples_color)

# After application of minimum UMI threshold
counts_minimum_umi <- counts %>%
  filter(bc.umi.count >= 3)

# # proportion of cells with a guide detected that have a single guide detected
# counts_minimum_umi %>%
#   group_by(sample, cellid) %>%
#   summarize(n_guides = n()) %>%
#   group_by(sample) %>%
#   summarize(n_cells_with_guide = n(), sum(n_guides == 1) / n_cells_with_guide)

p5 <- plotBarcodesPerCellSample(counts_minimum_umi,
                          fraction = F,
                          n_cells_sample = n_cells_sample) +
  scale_fill_manual(values = samples_color) +
  ggtitle("Number of guides detected per cell with ≥3 UMI")

p6 <- plotBarcodesPerCellSample(counts_minimum_umi,
                                fraction = T,
                                n_cells_sample = n_cells_sample) +
  scale_fill_manual(values = samples_color) +
  ggtitle("Number of guides detected per cell with ≥3 UMI")


p7 <- counts_minimum_umi %>%
  group_by(sample, barcode) %>%
  mutate(n_cells = n()) %>%
  filter(n_cells >= 10) %>%
  group_by(sample) %>%
  summarize(n_guides = as.numeric(n_distinct(barcode))) %>%
  ggplot(aes(x = sample, y = n_guides, fill = sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = samples_color) +
  ylab("Guides") +
  ggtitle("Number of guides detected in ≥10 cells") +
  theme_bw() +
  ggplot2::theme(panel.grid = element_blank(),
                 strip.background = element_blank(),
                 panel.border = element_blank(),
                 axis.line = element_line(), 
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed(ratio = 0.05)



pdf(
  file.path(
    plot_dir, 
    "MF03_guides_qc.pdf"
  ), 
  height = 5, 
  width = 9
)

p1
p2
p3
p4
p5
p6
p7

dev.off()
