#!/usr/bin/env Rscript

# DESCRIPTION
# Create QC metrics for Nanopore data:
# - how many cells have a ONT transcript detected, for each gene and sample
# - how many target gene transcripts per cell, split by capture
# - what proportion of ONT transcripts is from target genes, for 10X and ONT
# - Coverage of the transcript, showing internal primer position

# .libPaths("/home/hholze/R/x86_64-pc-linux-gnu-library/4.4")
# Sys.setenv(RESULT_DIR="/scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/manuscript/results/")
# Sys.setenv(PLOT_DIR="/dawson_genomics/Projects/PRC2_BE_screen/manuscript/figures/")

library(tidyverse)
library(Seurat)
library(patchwork)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

############################################################################################
# Plot UMI counts in 10X vs ONT data
### Combine nanopore and 10X data in a Seurat object

# load ONT gene counts, same data that's used for variant calling
nano_data_genes <- data.table::fread(paste0(result_dir, "MF03/nanopore/sicelore/output/04b.matrices/sicelore_genematrix.txt"))
rownames(nano_data_genes) <- nano_data_genes$geneId
nano_data_genes$geneId <- NULL

# turn into Seurat object
be_nano <- CreateSeuratObject(counts = nano_data_genes, assay = "RNA_nano_genes", min.cells = 0, min.features = 0)

# Only keep cells that pass 10X scRNA-seq QC.  
# The other cells are not of use since they are not demultiplexed into samples. 
be_sc <-
  qs::qread(
    paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs")
  )

be_sc@meta.data <- be_sc@meta.data %>%
  mutate(cellid = gsub("P[12]_([ACTG]*)-1", "\\1", rownames(be_sc@meta.data)))
be_sc@meta.data$cellid_pool <- rownames(be_sc@meta.data)

be_nano <- be_nano[, colnames(be_nano) %in% be_sc$cellid]

# Remove cell barcodes that are present in both captures. 
# Captures were pooled for ONT
barcodes_p1 <-
  read.table(
    paste0(result_dir, "MF03/scRNA_seq/cellranger_results/MF03_P1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
  )$V1
barcodes_p2 <-
  read.table(
    paste0(result_dir, "MF03/scRNA_seq/cellranger_results/MF03_P2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
  )$V1
barcodes_overlap <- gsub("-1", "", barcodes_p1[barcodes_p1 %in% barcodes_p2])
be_nano <- be_nano[,!(colnames(be_nano) %in% barcodes_overlap)]

# convert the cell barcode into the cell ID that includes the pool information to 
# merge ONT data with 10X data
be_nano@meta.data$cellid <- rownames(be_nano@meta.data)

metadata <-
  merge(
    be_nano@meta.data,
    be_sc@meta.data,
    by = "cellid",
    all.x = T,
    suffixes = c("", "_10x")
  )
rownames(metadata) <- metadata$cellid
metadata <- metadata[colnames(be_nano),]

be_nano <- RenameCells(be_nano, new.names = metadata$cellid_pool)

# add cells that are detected in 10X data and missing in ONT data 
# so we can merge the data into 1 object
cells_10x_only <- colnames(be_sc)[!colnames(be_sc) %in% colnames(be_nano)]

mtx_cells_10x_only <- Matrix::Matrix(nrow = nrow(be_nano), ncol = length(cells_10x_only), data = 0, sparse = T)
colnames(mtx_cells_10x_only) <- cells_10x_only
rownames(mtx_cells_10x_only) <- rownames(be_nano)

# add cells to matrix
nano_data_genes_combine_sc <- Matrix::cbind2(x = be_nano[["RNA_nano_genes"]]@counts, mtx_cells_10x_only)

# get cells in same order
nano_data_genes_combine_sc <- nano_data_genes_combine_sc[,colnames(be_sc)]

# add assay to 10X seurat object
be_sc[["RNA_nano_genes"]] <- CreateAssayObject(counts = nano_data_genes_combine_sc)

# Get the number of ONT target gene transcripts per cell, split by pool and target gene
nano_target_umi_counts <- be_sc@assays$RNA_nano_genes@counts[c("EZH2", "SUZ12", "EED"),]
nano_target_umi_counts <- cbind(t(as.matrix(nano_target_umi_counts)), be_sc@meta.data)
nano_target_umi_counts <- nano_target_umi_counts %>% 
  pivot_longer(cols = c("EZH2", "SUZ12", "EED"), names_to = "gene", values_to = "ONT_umi_count")

# Same for 10X
tenx_target_umi_counts <- be_sc@assays$RNA@counts[c("EZH2", "SUZ12", "EED"),]
tenx_target_umi_counts <- cbind(t(as.matrix(tenx_target_umi_counts)), be_sc@meta.data)
tenx_target_umi_counts <- tenx_target_umi_counts %>% 
  pivot_longer(cols = c("EZH2", "SUZ12", "EED"), names_to = "gene", values_to = "10X_umi_count")

# compare number of cells with target gene transcript detected between 10X and ONT
combined_target_umi_counts <- merge(nano_target_umi_counts, select(tenx_target_umi_counts, c(cellid_pool, gene, `10X_umi_count`)), by = c("gene", "cellid_pool"))

combined_target_umi_counts <- combined_target_umi_counts %>%
  pivot_longer(cols = c(`10X_umi_count`, ONT_umi_count), names_to = "technology", values_to = "umi_count") %>%
  mutate(technology = gsub("_umi_count", "", technology))

# get total number of demultiplexed cells per sample
n_cells_sample <- as.list(table(be_sc$sample))

# Get the proportion of transcripts that are for target genes
combined_target_umi_counts_summarized <- combined_target_umi_counts %>%
  group_by(gene, technology) %>%
  summarize(sum_umi_count = sum(umi_count))

sum_all_genes_10x_umi_count <- sum(be_sc@assays$RNA@counts)
sum_all_genes_ont_umi_count <- sum(be_sc@assays$RNA_nano_genes@counts)

combined_target_umi_counts_summarized <- merge(
  combined_target_umi_counts_summarized,
  data.frame(list(
    technology = c("10X", "ONT"),
    total_umi_count = c(sum_all_genes_10x_umi_count, sum_all_genes_ont_umi_count)
  )),
  by = "technology",
  all = T
) %>%
  mutate(fraction_gene_umi_count = sum_umi_count / total_umi_count) 

# Calculate remaining fraction per technology
df_remaining <- combined_target_umi_counts_summarized %>%
  group_by(technology) %>%
  summarise(remaining = 1 - sum(fraction_gene_umi_count)) %>%
  mutate(gene = "other genes") %>%
  rename(fraction_gene_umi_count = remaining)

# Combine original fractions + remaining
combined_target_umi_counts_summarized <- bind_rows(combined_target_umi_counts_summarized, df_remaining)
combined_target_umi_counts_summarized$gene_factor <- factor(combined_target_umi_counts_summarized$gene, levels = c("other genes", "EED", "EZH2", "SUZ12"), ordered = T)
combined_target_umi_counts_summarized <- combined_target_umi_counts_summarized %>%
  arrange(desc(gene_factor))

colors_genes <- scales::hue_pal()(3)
names(colors_genes) <- c("EED", "EZH2", "SUZ12")
colors_genes["other genes"] <- "gray"


############################################################################################
# Plots

# plot number of ONT target gene transcripts per cell, split by pool and target gene
p1 <- nano_target_umi_counts %>%
  ggplot(aes(x = pool, y = ONT_umi_count, fill = pool)) +
  geom_violin(scale = "width") +
  geom_boxplot(outliers = F, width = 0.1) +
  facet_wrap(~gene, axes = "all") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "None"
  ) +
  ylab("ONT UMI per cell")

# compare the average UMI count per cell between pools as bar plot
p2 <- nano_target_umi_counts %>%
  group_by(gene, pool) %>%
  summarize(mean_ONT_umi_count = mean(ONT_umi_count)) %>%
  ggplot(aes(x = gene, y = mean_ONT_umi_count, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Average ONT UMI count per cell") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "None"
  ) 

p3 <- nano_target_umi_counts %>%
  ggplot(aes(x = gene, y = ONT_umi_count)) +
  geom_violin(scale = "width") +
  facet_wrap(~sample, axes = "all") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "None"
  ) +
  ylab("ONT UMI per cell")


p4 <- combined_target_umi_counts %>%
  group_by(sample, gene, technology) %>%
  summarize(n_cells_transcript = sum(umi_count > 0)) %>%
  ggplot(aes(x = gene, y = n_cells_transcript, fill = technology)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggplot2::geom_hline(
    data = data.frame(
      sample = names(n_cells_sample),
      yintercept = as.numeric(n_cells_sample)
    ),
    aes(yintercept = yintercept, linetype = "# cells in sample"),
  ) +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_wrap(~sample, axes = "all", scales = "free_y") +
  ggtitle("Cells with transcript detected") +
  ylab("Cells") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

p5 <- combined_target_umi_counts %>%
  group_by(sample, gene, technology) %>%
  summarize(n_cells_transcript = sum(umi_count > 0), n_cells_sample = n(), frac_cells_transcript = n_cells_transcript / n_cells_sample) %>%
  ggplot(aes(x = gene, y = frac_cells_transcript, fill = technology)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggplot2::scale_linetype_manual(name = "", values = c(2)) +
  facet_wrap(~sample, axes = "all") +
  ggtitle("Cells with transcript detected") +
  ylab("Proportion of cells in sample") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# # plot UMI counts of each cell in 10X vs. ONT for the target genes
# p5 <- combined_target_umi_counts %>%
#   pivot_wider(names_from = technology, values_from = umi_count) %>%
#   ggplot(aes(x = ONT, y = `10X`)) +
#   geom_point(alpha = 0.1) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~gene, axes = "all", scales = "free", ncol = 2) +
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(linewidth = 0.5)
#   )

# Plot the proportion of transcripts that are for target genes

p6 <- combined_target_umi_counts_summarized %>%
  ggplot(aes(x = technology, y = fraction_gene_umi_count, fill = gene_factor)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors_genes) +
  facet_wrap(~technology, scales = "free_x", axes = "all") +
  ylab("% UMI") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
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
    "MF03_plots_nanopore_qc.pdf"
  ), 
  height = 6, 
  width = 7
)

p1 / p2
p3
p4
p5
p6

dev.off()


sum_target_genes_ont_umi_count <- sum(nano_target_umi_counts$ONT_umi_count)
sum_all_genes_ont_umi_count <- sum(be_sc@assays$RNA_nano_genes@counts)

sum_target_genes_10x_umi_count <- sum(tenx_target_umi_counts$`10X_umi_count`)
sum_all_genes_10x_umi_count <- sum(be_sc@assays$RNA@counts)

cat("UMI in ONT data that map to target genes", round(sum_target_genes_ont_umi_count / sum_all_genes_ont_umi_count * 100, 2), "%\n")
cat("UMI in 10X data that map to target genes", round(sum_target_genes_10x_umi_count / sum_all_genes_10x_umi_count * 100, 2), "%\n")
# UMI in ONT data that map to target genes 25.1 %
# UMI in 10X data that map to target genes 0.02 %


############################################################################################
# Coverage of the transcript

# cellsnp data that was filtered for cells that can be merged with scRNA-seq data
cellsnp_res_dir <-
  paste0(result_dir, "MF03/nanopore/")
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

vcf <- vcf %>%
  separate(INFO, into = c("AD", "DP", "OTH"), sep = ";") %>%
  mutate(
    AD = as.integer(gsub("AD=", "", AD)),
    DP = as.integer(gsub("DP=", "", DP)),
    OTH = as.integer(gsub("OTH=", "", OTH))
  )
# matrix with counts for cells that are demultiplexed into cells
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
rownames(mat_ref_alt_samples) <- positions

coverage_demultiplexed_cells <- Matrix::rowSums(mat_ref_alt_samples)

vcf$coverage_cells <- coverage_demultiplexed_cells[as.character(vcf$POS)]

plotCdsCoverage <- function(exons, gene_strand, vcf, gene) {
  exon_list <-
    str_split(gsub("chr[0-9]*:", "", str_split(exons, "\n")[[1]]), "-")
  exon_nc_positions <-
    unlist(sapply(exon_list, function(x)
      seq(as.integer(x[1]), as.integer(x[2]))))
  
  if (!all(exon_nc_positions %in% vcf$POS)) {
    cat("Not all positions of transcript are in vcf file")
  }
  if (gene_strand == "sense") {
    vcf_cds <- data.frame(POS = exon_nc_positions) %>%
      mutate(CDS = dense_rank(POS)) %>%
      select(POS, CDS)
  } else if (gene_strand == "antisense") {
    vcf_cds <- data.frame(POS = exon_nc_positions) %>%
      mutate(CDS = dense_rank(desc(POS))) %>%
      select(POS, CDS)
  }
  
  vcf <- merge(vcf, vcf_cds, by = "POS", all = T)
  
  coverage_df <- vcf %>%
    filter(!is.na(CDS)) %>%
    dplyr::select(CDS, DP, coverage_cells) %>%
    mutate(gene = gene)
  
  p_cov_total <- vcf %>%
    filter(!is.na(CDS)) %>%
    # DP column is ref + alt UMI count
    ggplot(aes(x = CDS, y = DP)) +
    geom_bar(stat = "identity",
             fill = colors_genes[[gene]],
             # no gaps between bars
             width = 1) +
    xlim(min(vcf$CDS, na.rm = T), max(vcf$CDS, na.rm = T)) +
    ylab("UMI") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "None"
    )
  
  p_cov_cells <- vcf %>%
    filter(!is.na(CDS)) %>%
    # coverage_cells column is ref + alt UMI count for demultiplexed cells
    ggplot(aes(x = CDS, y = coverage_cells)) +
    geom_bar(stat = "identity",
             fill = colors_genes[[gene]],
             # no gaps between bars
             width = 1) +
    xlim(min(vcf$CDS, na.rm = T), max(vcf$CDS, na.rm = T)) +
    ylab("UMI") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "None"
    )
  return(list(p_cov_total, p_cov_cells, coverage_df))
}

### Gene-specific parameters
# For EZH2, you can see in the coverage plot, that not NM_004456.5 but rather NM_001203247.2 is expresseed.
# There is a gap in one exon chr7:148826454-148826632 and the exon boundary of NM_001203247.2 fits better (chr7:148826469-148826632).
# I use the same transcript for VEP effect prediction, protein position etc.
# Of course it's bad if we use a transcript that is "too short" if there are additional exons expressed in a subset of cells. 
# But it is worse to use a transcript that isn't actually expressed with additional exons. Because then the coverage of that additional
# exon is low quality wrongly mapped transcripts with likely sequencing errors. 

# Get ranges of all exons of target gene to get coding sequence. 
# Exon sequences are downloaded from UCSC table browser, Dec 2013 GRCh38/hg38

# SUZ12   NM_015355.4    chr17:31,937,007-32,001,038
# EED     NM_003797.5    chr11:86,244,753-86,278,810
# EZH2    NM_001203247.2 chr7:148,807,383-148884291

# NOT using this one
# EZH2    NM_004456.5    chr7:148,807,383-148,884,291

# `grep "NM_015355.4" 2024-04-03_SUZ12_exons.fasta | grep -o "chr17:[0-9]*-[0-9]*"`  
# Double check UTR and CDS boundaries in UCSC browser. 

gene <- "SUZ12"

exons <- "chr17:31937247-31937520
chr17:31940286-31940332
chr17:31940422-31940486
chr17:31947617-31947685
chr17:31966147-31966196
chr17:31973146-31973231
chr17:31975482-31975713
chr17:31976521-31976614
chr17:31982999-31983104
chr17:31988320-31988497
chr17:31993242-31993333
chr17:31993865-31994008
chr17:31994564-31994721
chr17:31995564-31995762
chr17:31996798-31996877
chr17:31998658-31999003"

UTR_start <- 31937007
UTR_end <- 32001038
gene_range <- c(17, 31937007, 32001038)
target_transcript <- "NM_015355.4"
gene_strand <- "sense"

pp1 <- plotCdsCoverage(exons, gene_strand, vcf, gene)
pp1[[1]] <- pp1[[1]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))
pp1[[2]] <- pp1[[2]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))

#####

gene <- "EED"

exons <- "chr11:86245230-86245343
chr11:86250296-86250448
chr11:86252148-86252240
chr11:86255222-86255287
chr11:86256387-86256512
chr11:86257515-86257596
chr11:86264172-86264263
chr11:86266083-86266216
chr11:86268456-86268561
chr11:86276980-86277138
chr11:86277918-86277991
chr11:86278399-86278525"

UTR_start <- 86244753
UTR_end <- 86278810
gene_range <- c(11, 86244753, 86278810)
target_transcript <- "NM_003797.5"
gene_strand <- "sense"

pp2 <- plotCdsCoverage(exons, gene_strand, vcf, gene)
pp2[[1]] <- pp2[[1]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))
pp2[[2]] <- pp2[[2]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))

#####

# EZH2 is on antisense strand. 
gene <- "EZH2"

# exons <- "chr7:148847182-148847298
# chr7:148846470-148846598
# chr7:148832634-148832750
# chr7:148829728-148829848
# chr7:148828740-148828880
# chr7:148827164-148827266
# chr7:148826454-148826632
# chr7:148819596-148819687
# chr7:148817877-148818117
# chr7:148817222-148817391
# chr7:148816684-148816778
# chr7:148815506-148815546
# chr7:148814914-148815039
# chr7:148813959-148814137
# chr7:148811625-148811720
# chr7:148810333-148810414
# chr7:148809310-148809390
# chr7:148809071-148809155
# chr7:148807646-148807706"

exons <- "chr7:148847182-148847298
chr7:148846470-148846598
chr7:148832634-148832750
chr7:148829728-148829848
chr7:148828740-148828880
chr7:148827164-148827266
chr7:148826469-148826632
chr7:148819596-148819687
chr7:148817877-148818117
chr7:148817222-148817391
chr7:148816684-148816778
chr7:148815506-148815546
chr7:148814914-148815039
chr7:148813959-148814137
chr7:148811625-148811720
chr7:148810333-148810414
chr7:148809310-148809390
chr7:148809071-148809155
chr7:148807646-148807706"

UTR_start <- 148807383
UTR_end <- 148884291
gene_range <- c(7, 148807383, 148884291)
# target_transcript <- "NM_004456.5"
target_transcript <- "NM_001203247.2"
gene_strand <- "antisense"

pp3 <- plotCdsCoverage(exons, gene_strand, vcf, gene)
pp3[[1]] <- pp3[[1]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))
pp3[[2]] <- pp3[[2]] +
  ggtitle(paste0("CDS ", gene, " (", target_transcript, ")"))


pdf(
  file.path(
    plot_dir, 
    "MF03_plots_nanopore_cds_coverage.pdf"
  ), 
  height = 3, 
  width = 20
)

# Plot total coverage in ONT data
pp1[[1]] | pp2[[1]] | pp3[[1]]
# Plot coverage only for cells that were demultiplexed into samples
pp1[[2]] | pp2[[2]] | pp3[[2]]

dev.off()

write_tsv(
  bind_rows(pp1[[3]], pp2[[3]], pp3[[3]]),
  paste0(
    result_dir,
    "MF03/nanopore/MF03_nanopore_cds_coverage.tsv"
  )
)
