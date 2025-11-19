#!/usr/bin/env Rscript

# DESCRIPTION
# MHC1 library QC
# Goal: Check quality, i.e. skewednes of the library.


library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")


guide_counts <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_guide_library/bartab/results_18bp/counts/all_counts_combined.tsv"
    )
  )
colnames(guide_counts) <- c("sgRNA_ID", "count")

guide_library <-
  read_tsv(
    paste0(
      result_dir,
      "MHC1_base_editing_library_EZH2_EED_SUZ12.tsv"
    )
  )
nrow(guide_library)

# All guides are detected in the sequenced library. 
guide_counts %>%
  summarize(n())


# Skewness: 90th / 10th percentile. Would be good if it was <2. 
CalculateSkewness <- function(counts, thresholds = c(0.9, 0.1)) {
  # return(3 * (mean(counts) - median(counts)) / sd(counts))
  return(as.vector(unlist(
    quantile(x = counts, probs = thresholds[1]) / quantile(x = counts, probs = thresholds[2])
  )))
}

thresholds_list <- list(c(0.9, 0.1), c(0.95, 0.05), c(0.99, 0.01))
skewness <- sapply(thresholds_list, function(thresholds) {
  CalculateSkewness(guide_counts %>% pull(count),
                    thresholds)
})
names(skewness) <- c("90/10", "95/5", "99/1")
skewness

# 90/10      95/5      99/1 
# 2.925484  4.930939 18.665291



p1 <- guide_counts %>%
  ggplot(mapping = aes(x = reorder(sgRNA_ID, -count), y = count)) +
  geom_bar(stat = "identity") +
  xlab("guides") +
  ylab("counts") +
  ggtitle("Distribution of counts in MHC1_BE_library", subtitle = paste0("Skewness 90th/10th percentile: ", round(skewness["90/10"], 2))) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

pdf(
  file.path(
    plot_dir, 
    "MHC1_library_skewness_18bp.pdf"
  ), 
  height = 3, 
  width = 5
)
p1
dev.off()
