#!/usr/bin/env Rscript

# DESCRIPTION
# Intersection of detected guides between samples
# Options: all cells with guide or only cells with single guide, filtered for EED KO score or not
# For filtereing EED KO score, take guides that have a significantly higher EED KO score than Ctrl HLA+ sample.
# Calculate adjusted p-value on the fly, depending on which samples to show in the plot.

# INPUT
# - table with guides in samples and EED KO score t-test results, filtered for at least 10 cells 
#     with the guide in the sample with 3 UMI
# - same as above but removed cells with multiple guides integrated before selecting guides with 
#     at least 10 cells per sample and doing t-test

library(UpSetR)
library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "") dir.create(plot_dir, recursive = TRUE)

## Load data
ttest_samples_guides_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

ttest_samples_guides_filtered_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_single_guide_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

# plotting function
plotGuideUpset <- function(ttest_samples_guides_df,
         samples,
         filter_eed_ko = F) {
  
  ttest_samples_guides_df <- ttest_samples_guides_df %>%
    filter(sample %in% samples)
  if (filter_eed_ko) {
    n_tests <- nrow(ttest_samples_guides_df)
    
    ttest_samples_guides_df <- ttest_samples_guides_df %>%
      filter(
        diff_eed_ko_score_normalized > 0.06,
        p.value_eed_ko_score_normalized * n_tests < 0.05
      )
  }
  
  guide_lists <-
    lapply(samples, function(x) {
      ttest_samples_guides_df %>% filter(sample == x) %>% pull(guide) %>% unique()
    })
  names(guide_lists) <- samples
  
  UpSetR::upset(
    fromList(guide_lists),
    nsets = 8,
    order.by = "freq",
    nintersects = NA
  )
}

# select only base editor samples to show
samples_be <- c("V7-HLA-pos",
                "TADA-HLA-pos",
                "CDA1-HLA-pos",
                "V14-HLA-pos",
                "RAPO-HLA-pos")


# do for all samples and including cells with multiple guides detected
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_df,
  samples = samples_be,
  filter_eed_ko = F
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p1 <- gridExtra::grid.arrange(grobs = list(vp), top = "Including multiple guides per cell")

# do for all samples and only cells with a single guides detected
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_filtered_df,
  samples = samples_be,
  filter_eed_ko = F
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p2 <- gridExtra::grid.arrange(grobs = list(vp), top = "Single guide per cell")

# do for all samples including cells with multiple guides detected 
# and only guides with EED KO score sign. higher than ctrl HLA+
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_df,
  samples = samples_be,
  filter_eed_ko = T
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p3 <- gridExtra::grid.arrange(grobs = list(vp), top = "Including multiple guides per cell, high EED KO score")

# do for all samples and only cells with a single guides detected 
# and only guides with EED KO score sign. higher than ctrl HLA+
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_filtered_df,
  samples = samples_be,
  filter_eed_ko = T
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p4 <- gridExtra::grid.arrange(grobs = list(vp), top = "Single guide per cell, high EED KO score")

####################################
####################################

# The whole thing again w/o V14
samples_be <- c("V7-HLA-pos", "TADA-HLA-pos", "CDA1-HLA-pos", "RAPO-HLA-pos")

# do for samples excl V14 and including cells with multiple guides detected
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_df,
  samples = samples_be,
  filter_eed_ko = F
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p5 <- gridExtra::grid.arrange(grobs = list(vp), top = "Including multiple guides per cell")

# do for samples excl V14 and only cells with a single guides detected
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_filtered_df,
  samples = samples_be,
  filter_eed_ko = F
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p6 <- gridExtra::grid.arrange(grobs = list(vp), top = "Single guide per cell per cell")

# do for samples excl V14 including cells with multiple guides detected 
# and only guides with EED KO score sign. higher than ctrl HLA+
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_df,
  samples = samples_be,
  filter_eed_ko = T
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p7 <- gridExtra::grid.arrange(grobs = list(vp), top = "Including multiple guides per cell, high EED KO score")


# do for samples excl V14 and only cells with a single guides detected 
# and only guides with EED KO score sign. higher than ctrl HLA+
plotGuideUpset(
  ttest_samples_guides_df = ttest_samples_guides_filtered_df,
  samples = samples_be,
  filter_eed_ko = T
)
grid::grid.edit('arrange', name = "movies")
vp <- grid::grid.grab()
p8 <- gridExtra::grid.arrange(grobs = list(vp), top = "Single guide per cell, high EED KO score")



pdf(
  file.path(
    plot_dir, 
    "MF03_guide_overlap.pdf"
  ), 
  height = 4, 
  width = 5
)

plot(p1)
plot(p2)
plot(p3)
plot(p4)
plot(p5)
plot(p6)
plot(p7)
plot(p8)

dev.off()
