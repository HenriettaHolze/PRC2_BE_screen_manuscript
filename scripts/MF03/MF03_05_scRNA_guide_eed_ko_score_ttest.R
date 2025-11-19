#!/usr/bin/env Rscript

# DESCRIPTION
# - Identify guides detected in a sample as detected with >=3 UMI and detected in at least 10 cells
# - Perform t-test for each guide-BE-combination against the HLA+ Ctrl sample for the EED KO score
# - On the go, get the average EED, SUZ12, EZH2 expression per guide
# - Do the same t-test once for cells with guides annotated above UMI threshold and once for cells
#     with only a single guide above that threshold.

# INPUT
# - Seurat object with guide counts as assay and EED KO score in metadata
# - Guide count table unfiltered and filtered

# OUTPUT
# - List of all detected guides per sample and their EED KO score

library(Seurat)
library(tidyverse)
library(ltm)
library(qs)

result_dir <- Sys.getenv("RESULT_DIR")
dir.create(file.path(result_dir, "MF03/scRNA_seq/"), recursive = TRUE)

## Load data
be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )

guide_counts <- read_tsv(paste0(result_dir, "MF03/scRNA_seq/MF03_guide_counts.tsv"))
guide_counts_single_integration <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_single_guide.tsv"
))
guide_counts_aggregated <- read_tsv(paste0(
  result_dir,
  "MF03/scRNA_seq/MF03_guide_counts_aggregated.tsv"
))

## t-test

# Define samples of interest
samples <- c(
  "V7-HLA-pos",
  "V14-HLA-pos",
  "TADA-HLA-pos",
  "CDA1-HLA-pos",
  "RAPO-HLA-pos",
  "CAS9-HLA-pos",
  "Ctrl-HLA-pos"
)

# Get Ctrl HLA+ sample phenotype as the reference to compare to
# The score was normalized to the Ctrl HLA+ sample, so by definition
# the average EED KO score of the Ctrl HLA+ sample is 0
ctrl_pheno_eed_ko_score_normalized <- be_singlets@meta.data %>%
  filter(sample == "Ctrl-HLA-pos") %>%
  pull(eed_ko_score_normalized)

# Function to perform t-test between ctrl sample and cells with guides
ttestSamplesGuides <- function(be_singlets,
                               guide_counts_filtered,
                               ctrl_pheno_eed_ko_score_normalized) {
  samples <- unique(be_singlets$sample)
  # iterate over samples
  ttest_samples_guides_list <-
    sapply(
      samples,
      FUN = function(selected_sample) {
        # subset seurat object to sample
        expr <- FetchData(object = be_singlets, vars = "sample")
        be_singlets_sample <- be_singlets[, which(x = expr == selected_sample)]
        
        # get guides detected in minimum number of cells with minimum UMI
        sample_detected_guides <- guide_counts_filtered %>%
          filter(sample == selected_sample) %>%
          pull(barcode) %>%
          unique()
        
        cat(paste0(
          "Sample ",
          selected_sample,
          " has ",
          length(sample_detected_guides),
          " guides detected\n"
        ))
        
        # get metadata for sample
        sample_metadata <-
          be_singlets@meta.data[be_singlets@meta.data$sample == selected_sample, ]
        
        # iterate over guides in that sample
        sapply(
          sample_detected_guides,
          FUN = function(guide) {
            cell_with_guide <- guide_counts_filtered %>%
              filter(sample == selected_sample, barcode == guide) %>%
              pull(cell_id)
            
            # need to make sure that it's not expression level of PRC2 core members that drives HLA- condition
            # calculate mean normalized expression of target genes
            eed_mean <-
              mean(be_singlets_sample[["RNA"]]@data["EED", cell_with_guide])
            ezh2_mean <-
              mean(be_singlets_sample[["RNA"]]@data["EZH2", cell_with_guide])
            suz12_mean <-
              mean(be_singlets_sample[["RNA"]]@data["SUZ12", cell_with_guide])
            
            # get phenotype vector of cells in sample with guide detected
            sample_guide_pheno_eed_ko_score_normalized <-
              sample_metadata[cell_with_guide, ] %>%
              pull(eed_ko_score_normalized)
            
            # perform t-test, to test whether means of K562 Ctrl and guides differ
            # by default paired = FALSE, var.equal = FALSE
            ttest_res_eed_ko_score_normalized <-
              t.test(
                ctrl_pheno_eed_ko_score_normalized,
                sample_guide_pheno_eed_ko_score_normalized
              )
            
            # get the difference of ctrl minus guide mean score (effect size) and p-value from t-test results
            return(
              list(
                # this is the same as the average EED KO score for that sample, because the mean EED KO score
                # of the Ctrl HLA+ sample is 0
                diff_eed_ko_score_normalized = as.numeric(
                  ttest_res_eed_ko_score_normalized$estimate[2] - ttest_res_eed_ko_score_normalized$estimate[1]
                ),
                p.value_eed_ko_score_normalized = ttest_res_eed_ko_score_normalized$p.value,
                
                n_cells = length(sample_guide_pheno_eed_ko_score_normalized),
                eed_mean = eed_mean,
                ezh2_mean = ezh2_mean,
                suz12_mean = suz12_mean
              )
            )
            
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
  
  # Make one big dataframe with all guides and all samples.
  ttest_samples_guides_df <-
    bind_rows(purrr::map(ttest_samples_guides_list, bind_rows, .id = "guide"),
              .id = "sample")
  
  ttest_samples_guides_df <- ttest_samples_guides_df %>%
    mutate(guide = gsub("_", "-", guide))
  
  # Adjust p-value by multiplying with number of tests
  ttest_samples_guides_df <- ttest_samples_guides_df %>%
    mutate(p.val.adj_eed_ko_score_normalized = p.value_eed_ko_score_normalized * nrow(ttest_samples_guides_df))
  
}

# Set thresholds what is considered a "detected guide" in a sample
umi_threshold <- 3
n_cells_threshold <- 10

guide_counts_filtered <- guide_counts %>%
  filter(bc.umi.count >= umi_threshold) %>%
  group_by(barcode, sample) %>%
  mutate(n_cells = n()) %>%
  filter(n_cells >= n_cells_threshold) %>%
  dplyr::select(-n_cells)

ttest_samples_guides_df <- ttestSamplesGuides(be_singlets,
                                              guide_counts_filtered,
                                              ctrl_pheno_eed_ko_score_normalized)

# Arrange by average EED KO score
ttest_samples_guides_df <- ttest_samples_guides_df %>%
  arrange(desc(diff_eed_ko_score_normalized))

write_tsv(
  ttest_samples_guides_df,
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)


# Do the same thing, but this time with the guide annotation without cells that have multiple guides detected
ttest_samples_guides_filtered_df <- ttestSamplesGuides(
  be_singlets,
  guide_counts_single_integration,
  ctrl_pheno_eed_ko_score_normalized
)

# Arrange by average EED KO score
ttest_samples_guides_filtered_df <- ttest_samples_guides_filtered_df %>%
  arrange(desc(diff_eed_ko_score_normalized))

write_tsv(
  ttest_samples_guides_filtered_df,
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_single_guide_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)


# Do the same thing, but this time with the guide annotation with guides aggregated per cell
ttest_samples_guides_aggregated_df <- ttestSamplesGuides(
  be_singlets,
  guide_counts_aggregated,
  ctrl_pheno_eed_ko_score_normalized
)

# Arrange by average EED KO score
ttest_samples_guides_aggregated_df <- ttest_samples_guides_aggregated_df %>%
  arrange(desc(diff_eed_ko_score_normalized))

write_tsv(
  ttest_samples_guides_aggregated_df,
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_aggregated_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

