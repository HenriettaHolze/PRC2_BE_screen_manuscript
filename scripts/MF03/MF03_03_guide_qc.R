#!/usr/bin/env Rscript

# DESCRIPTION
# - Demultiplex guide counts into samples according to scRNA-seq demultiplexing results, save as table.
# - add as assay to Seurat object
# - Filter guides for at least 3 UMI per cell and aggregate barcodes detected in a cell, 
#     remove barcode combinations detected in less than 10 cells per sample, save as table.
# - Remove cells with multiple guides in cell detected, save as table.
# - Store filtered guides as metadata slot in Seurat assay.
# - Calculate jaccard index between guides detected in a sample, (guides with >= 10 cells)
# NB: Guides are sometimes referred to as barcodes, because I use the BARtab/bartools packages to
# analyze the data which was developed for analysis of cellular barcoding data.

# INPUT
# - Seurat object with demultiplexed cells
# - Guide count tables from BARtab results

# OUTPUT
# - Guide counts per cell with sample annotation
# - Seurat object with guide count as assay and filtered guide annotation as metadata slot
# - Filtered guide counts per cell with sample annotation
# - Jaccard index between guides detected in a sample

library(Seurat)
library(tidyverse)
library(qs)
library(data.table)
result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data

# Demultiplexes cells
be_singlets <-
  qs::qread(
    paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap.qs")
  )

# BARtab results
counts_p1 <- utils::read.delim(paste0(result_dir, "MF03/scRNA_seq/bartab_guides/bartab_results/counts/MF03-1-CRISPR.counts.tsv"), header = TRUE)
counts_p1 <- counts_p1[, c("cell", "gene", "count")]
colnames(counts_p1) <- c("cellid", "barcode", "bc.umi.count")

counts_p2 <- utils::read.delim(paste0(result_dir, "MF03/scRNA_seq/bartab_guides/bartab_results/counts/MF03-2-CRISPR.counts.tsv"), header = TRUE)
counts_p2 <- counts_p2[, c("cell", "gene", "count")]
colnames(counts_p2) <- c("cellid", "barcode", "bc.umi.count")


## Demultiplex BARtab counts

# By merging them with the Seurat object
cell_barcode_table <- be_singlets@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_barcode = gsub("P._(.*)-.", "\\1", cell_id)) %>%
  dplyr::select(cell_id, cell_barcode, pool, hash.ID, sample)


# Split BARtab counts into samples
counts_p1 <-
  merge(
    counts_p1,
    filter(cell_barcode_table, pool == "P1"),
    by.x = "cellid",
    by.y = "cell_barcode"
  )
counts_p2 <-
  merge(
    counts_p2,
    filter(cell_barcode_table, pool == "P2"),
    by.x = "cellid",
    by.y = "cell_barcode"
  )

counts <- rbind(counts_p1, counts_p2)

# Save to file for later use
write_tsv(
  counts,
  paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts.tsv")
)

# aggregate all guides detected per cell with >= 3 UMI
# and keep barcode combinations present in at least 10 cells per sample
counts_aggregated <- counts %>%
  filter(bc.umi.count >= 3) %>%
  dplyr::arrange(.data$barcode)

counts_aggregated <- as.data.frame(data.table::dcast(
  data.table::setDT(counts_aggregated),
  cell_id ~ .,
  value.var = c("barcode", "bc.umi.count", "sample", "hash.ID", "cellid"),
  fun.aggregate = function(x) {
    paste(unique(x), collapse = ";")
  }
))

# keep barcode combinations that are present in at least 10 cells per sample
# with >=3 UMI
counts_aggregated <- counts_aggregated %>%
  group_by(barcode, sample) %>% 
  mutate(n_cells = n()) %>%
  filter(n_cells >= 10) %>%
  dplyr::select(-n_cells)

write_tsv(
  counts_aggregated,
  paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts_aggregated.tsv")
)

# remove cells that have multiple guides with >= 3 UMI
counts_single_guide <- counts_aggregated %>%
  filter(!grepl(";", barcode))

write_tsv(
  counts_single_guide,
  paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts_single_guide.tsv")
)

### Add guide counts to Seurat object

# Add the unfiltered UMI count matrix as an assay
counts_matrix <- counts %>%
  dplyr::select(barcode, cell_id, bc.umi.count) %>%
  pivot_wider(
    id_cols = barcode,
    names_from = cell_id,
    values_from = bc.umi.count,
    values_fill = 0
  )

# cells not in guide counts
cells_missing <-
  colnames(be_singlets)[!colnames(be_singlets) %in% colnames(counts_matrix)]

counts_matrix <-
  cbind(counts_matrix, setNames(lapply(cells_missing, function(x)
    x = 0), cells_missing))

rownames(counts_matrix) <- counts_matrix$barcode

counts_matrix <-
  counts_matrix[, colnames(be_singlets)]

be_singlets[["guides"]] <-
  CreateAssayObject(data = as.matrix(counts_matrix))

# Add the filtered guide annotation as a metadata slot
counts_single_guide_metadata <- counts_single_guide %>% 
  column_to_rownames("cell_id") %>%
  dplyr::rename(guide = barcode) %>%
  dplyr::select(-cellid, -hash.ID, -sample)

be_singlets <- AddMetaData(be_singlets, metadata = counts_single_guide_metadata)

qs::qsave(
  be_singlets,
  paste0(result_dir, "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides.qs")
)


## Calculate Jaccard index between guides
# When we remove cells with multiple guides detected, we get rid of a lot of ambiguity.
# Removing cells instead of guides has the advantage, that we keep cells with only 1 guide integrated
# if only one subclone of a guide has other guides integrated too.
# However, if we have poor guide detection for specific guides, e.g. due to poor PCR amplification, 
# then we can still have lots of cells where only 1 guide is detected, although on a guide level we 
# can see that there are multiple guides.


jaccard_df_list = list()
for (selected_sample in unique(counts$sample)) {
  
  # filter for minimal UMI threshold
  guide_list <- counts %>%
    filter(sample == selected_sample) %>%
    filter(bc.umi.count >= 3) %>%
    group_by(barcode) %>%
    summarise(cells = list(unique(cell_id))) %>%
    deframe()
  
  # we don't focus on guides with less than 10 cells anyways
  # But if a guide is in 10 cells and has another guide in 9 out of them, we want to know
  guide_present <- counts %>%
    filter(sample == selected_sample) %>%
    filter(bc.umi.count >= 3) %>%
    group_by(barcode, sample) %>%
    mutate(n_cells = n()) %>%
    filter(n_cells >= 10) %>%
    pull(barcode) %>% unique()
  
  # Generate all pairwise combinations of guide, where at least one guide is present in at least 10 cells
  guide_pairs <- list()
  counter <- 1
  for (i in guide_present) {
    for (j in names(guide_list)) {
      if (j %in% guide_present & j < i)
        next
      if (i != j) {
        guide_pairs[[counter]] <- c(i, j)
        counter <- counter + 1
      }
    }
  }
  
  cat(selected_sample, length(guide_pairs), "\n")
  
  if (length(guide_pairs) == 0) next
  
  # calculate 3 different measures of overlap and containment for each pair
  jaccard_results <- lapply(guide_pairs, function(pair) {
    set1 <- guide_list[[pair[1]]]
    set2 <- guide_list[[pair[2]]]
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    
    jaccard <- ifelse(union == 0, NA, intersection / union)
    dice <- ifelse(union == 0, NA, (2 * intersection) / (length(set1) + length(set2)))
    overlap <- ifelse(union == 0, NA, intersection / min(length(set1), length(set2)))
    
    data.frame(
      guide_a = pair[1],
      guide_b = pair[2],
      guide_a_cells = length(set1),
      guide_b_cells = length(set2),
      intersection = intersection,
      jaccard_index = jaccard,
      dice_coefficient = dice,
      overlap_coefficient = overlap
    )
  })
  
  # Combine results into a data frame
  jaccard_df <- do.call(rbind, jaccard_results) %>%
    arrange(desc(jaccard_index)) %>%
    mutate(sample = selected_sample)
  
  jaccard_df_list[[selected_sample]] <- jaccard_df
}

# combine all samples
jaccard_samples_df <- do.call(rbind, jaccard_df_list)

# reorder columns
jaccard_samples_df <- jaccard_samples_df %>%
  select(sample, guide_a, guide_b, guide_a_cells, guide_b_cells, intersection, jaccard_index, dice_coefficient, overlap_coefficient)

# make if A-B and B-A in the table, for easier merging with other tables and filtering in excel
jaccard_samples_df_ba <- jaccard_samples_df %>%
  dplyr::rename(
    guide_a = guide_b,
    guide_b = guide_a,
    guide_a_cells = guide_b_cells,
    guide_b_cells = guide_a_cells
  )
jaccard_samples_df_abba <- rbind(jaccard_samples_df, jaccard_samples_df_ba) %>%
  arrange(desc(jaccard_index), sample, intersection) %>%
  filter(jaccard_index > 0)

# save file
write_tsv(
  jaccard_samples_df_abba,
  paste0(result_dir, "MF03/scRNA_seq/MF03_guides_jaccard.tsv")
)
