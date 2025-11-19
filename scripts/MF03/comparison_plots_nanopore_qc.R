
library(tidyverse)
library(Seurat)
library(patchwork)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

cds_coverage_mf01 <- read_tsv(
  paste0(
    result_dir,
    "MF01/nanopore/MF01_nanopore_cds_coverage.tsv"
  )
)
cds_coverage_mf01$experiment <- "MF01"

cds_coverage_mf03 <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_nanopore_cds_coverage.tsv"
  )
)
cds_coverage_mf03$experiment <- "MF03"

# combine experiments
cds_coverage <- bind_rows(cds_coverage_mf01, cds_coverage_mf03)

# plot total UMI count in experiments
p1 <- cds_coverage %>%
  ggplot(aes(x = CDS, y = DP / 10000)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free") +
  ylab("UMI (10k)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

p1b <- cds_coverage %>%
  ggplot(aes(x = CDS, y = DP / 10000)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free_x") +
  ylab("UMI (10k)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# plot proportion of maximum coverage per gene per experiment, showing the "evenness"
p2 <- cds_coverage %>%
  group_by(gene, experiment) %>%
  mutate(max_coverage = max(DP, na.rm = T)) %>%
  ungroup() %>%
  mutate(frac_max_cov = DP / max_coverage * 100) %>%
  ggplot(aes(x = CDS, y = frac_max_cov)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free_x") +
  ylab("% max coverage gene") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# plot proportion of maximum coverage per experiment, showing the "evenness"
p3 <- cds_coverage %>%
  group_by(experiment) %>%
  mutate(max_coverage = max(DP, na.rm = T)) %>%
  ungroup() %>%
  mutate(frac_max_cov = DP / max_coverage * 100) %>%
  ggplot(aes(x = CDS, y = frac_max_cov)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free_x") +
  ylab("% max coverage") +
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
    "comparison_MF01_MF03_plots_nanopore_cds_coverage.pdf"
  ), 
  height = 7, 
  width = 8
)

p1 / p1b / p2 / p3 

dev.off()


cds_coverage_3prime <- cds_coverage %>%
  group_by(gene, experiment) %>%
  filter(!is.na(DP)) %>%
  slice_max(order_by = CDS, n = 1) %>%
  mutate(cov_3prime = DP) %>%
  select(cov_3prime, gene, experiment)

p4 <- merge(cds_coverage, cds_coverage_3prime, by = c("gene", "experiment")) %>%
  mutate(cov_norm_to_3prime = DP / cov_3prime) %>%
  ggplot(aes(x = CDS, y = cov_norm_to_3prime)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free") +
  ylab("coverage normalized to \nmost 3' position") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

# cds_coverage %>%
#   dplyr::filter(gene == "EZH2", experiment == "MF01") %>%
#   filter(!is.na(DP)) %>%
#   arrange(desc(CDS))%>%
#   head(100) %>%
#   tail()

cds_coverage_3prime <- cds_coverage %>%
  group_by(gene, experiment) %>%
  filter(!is.na(DP)) %>%
  slice_max(order_by = CDS, n = 50) %>%
  slice_max(order_by = coverage_cells, n = 1) %>%
  mutate(cov_3prime = DP) %>%
  select(cov_3prime, gene, experiment)

p5 <- merge(cds_coverage, cds_coverage_3prime, by = c("gene", "experiment")) %>%
  mutate(cov_norm_to_3prime = DP / cov_3prime) %>%
  ggplot(aes(x = CDS, y = cov_norm_to_3prime)) +
  geom_freqpoly(aes(color = experiment), stat="identity") +
  facet_wrap(~gene, ncol = 3, axes = "all", scales = "free_x") +
  ylab("coverage normalized \nto 3' primer") +
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
    "comparison_MF01_MF03_plots_nanopore_cds_coverage_norm_3prime.pdf"
  ), 
  height = 3.5, 
  width = 8
)

p4 / p5

dev.off()
