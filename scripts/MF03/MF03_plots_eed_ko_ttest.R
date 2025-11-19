#!/usr/bin/env Rscript

# DESCRIPTION
# Plot the results of the t-test, comparing the EED KO score of each guide against the HLA+ ctrl sample.
# For cells with only single guide detected or including cells with multiple guides detected.


library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)

## Load data
ttest_samples_guides_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)
ttest_samples_guides_single_guide_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_single_guide_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)
ttest_samples_guides_aggregated_df <- read_tsv(
  paste0(
    result_dir,
    "MF03/scRNA_seq/MF03_ranked_guides_aggregated_vs_ctrl_eed_ko_score_ttest.tsv"
  )
)

be_singlets <-
  qs::qread(
    paste0(
      result_dir,
      "MF03/scRNA_seq/MF03_be_demultiplexed_cc_umap_guides_eed_ko_score.qs"
    )
  )

## color scheme
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

eed_ko_ctrl_diff <- be_singlets$eed_ko_score_normalized[be_singlets$sample == "Ctrl-HLA-pos"]


## Plot results of EED KO score t-test

# shuffle so that points don't hide each other by sample
set.seed(1)
ttest_samples_guides_df_shuffled <- ttest_samples_guides_single_guide_df %>%
  slice_sample(prop = 1)

ylim_max <- max(-log10(
  ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized[ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized != 0]
))
ylim_max <- ylim_max + 0.1 * ylim_max

xlim_min_max <-
  c(
    min(
      mean(eed_ko_ctrl_diff) - sd(eed_ko_ctrl_diff),
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ),
    max(
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ) + 0.05
  )

n_tests <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

p_be <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  arrange(p.value_eed_ko_score_normalized) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(
    x = diff_eed_ko_score_normalized,
    # p-value adjustment on the fly
    y = -log10(p.value_eed_ko_score_normalized * n_tests),
    color = sample,
  )) +
  geom_point(
    alpha = 0.5,
    shape = 19,
    stroke = 0,
    size = 2
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  ylab("-log10 adjusted p-value") +
  xlab("Average EED KO score per guide") +
  ggrepel::geom_text_repel(
    # label the top10 guides by p-value
    data = . %>% mutate(label = ifelse(rank <= 10, guide, "")),
    aes(label = label),
    max.overlaps = 100,
    min.segment.length = 0,
    nudge_x = 0.02,
    force_pull = 0.1,
    force = 1,
    show_guide  = FALSE,
    size = 3
  ) +
  coord_cartesian(xlim = xlim_min_max, ylim = c(0, ylim_max)) +
  ggtitle("EED KO score guides compared to Ctrl HLA+", subtitle = "Only single guide integration") +
  scale_color_manual(values = samples_color) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    panel.border = element_blank(),
    axis.line = element_line()
  )

p_ctrl_cas9 <- be_singlets@meta.data %>%
  filter(sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = eed_ko_score_normalized)) +
  stat_density(
    aes(color = sample, fill = sample),
    geom = "area",
    position = "identity",
    alpha = 0.5
  ) +
  coord_cartesian(xlim = xlim_min_max) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  scale_color_manual(values = unlist(samples_color)) +
  scale_fill_manual(values = unlist(samples_color)) +
  xlab("EED KO score per cell")

p1 <- ggpubr::ggarrange(
  p_be,
  p_ctrl_cas9,
  align = "v",
  heights = c(2, 1),
  ncol = 1,
  nrow = 2,
  common.legend = F
)


#######################################
# Do the same but including cells with multiple guides

ttest_samples_guides_df_shuffled <- ttest_samples_guides_df %>%
  slice_sample(prop = 1)

ylim_max <- max(-log10(
  ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized[ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized != 0]
))
ylim_max <- ylim_max + 0.1 * ylim_max

xlim_min_max <-
  c(
    min(
      mean(eed_ko_ctrl_diff) - sd(eed_ko_ctrl_diff),
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ),
    max(
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ) + 0.05
  )

n_tests <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

p_be <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  arrange(p.value_eed_ko_score_normalized) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(
    x = diff_eed_ko_score_normalized,
    # p-value adjustment on the fly
    y = -log10(p.value_eed_ko_score_normalized * n_tests),
    color = sample,
  )) +
  geom_point(
    alpha = 0.5,
    shape = 19,
    stroke = 0,
    size = 2
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  ylab("-log10 adjusted p-value") +
  xlab("Average EED KO score per guide") +
  ggrepel::geom_text_repel(
    # label the top10 guides by p-value
    data = . %>% mutate(label = ifelse(rank <= 10, guide, "")),
    aes(label = label),
    max.overlaps = 100,
    min.segment.length = 0,
    nudge_x = 0.02,
    force_pull = 0.1,
    force = 1,
    show_guide  = FALSE,
    size = 3
  ) +
  coord_cartesian(xlim = xlim_min_max, ylim = c(0, ylim_max)) +
  ggtitle("EED KO score guides compared to Ctrl HLA+", subtitle = "Including multiple guide integration") +
  scale_color_manual(values = samples_color) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    panel.border = element_blank(),
    axis.line = element_line()
  )

p_ctrl_cas9 <- be_singlets@meta.data %>%
  filter(sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = eed_ko_score_normalized)) +
  stat_density(
    aes(color = sample, fill = sample),
    geom = "area",
    position = "identity",
    alpha = 0.5
  ) +
  coord_cartesian(xlim = xlim_min_max) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  scale_color_manual(values = unlist(samples_color)) +
  scale_fill_manual(values = unlist(samples_color)) +
  xlab("EED KO score per cell")


p2 <- ggpubr::ggarrange(
  p_be,
  p_ctrl_cas9,
  align = "v",
  heights = c(2, 1),
  ncol = 1,
  nrow = 2,
  common.legend = F
)

pdf(
  file.path(plot_dir, "MF03_eed_ko_score_ttest.pdf"),
  height = 7,
  width = 6
)

p1
p2

dev.off()



#######################################
# Do the same for guide combinations (aggregated guides per cell)

ttest_samples_guides_df_shuffled <- ttest_samples_guides_aggregated_df %>%
  slice_sample(prop = 1)

ylim_max <- max(-log10(
  ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized[ttest_samples_guides_df_shuffled$p.value_eed_ko_score_normalized != 0]
))
ylim_max <- ylim_max + 0.1 * ylim_max

xlim_min_max <-
  c(
    min(
      mean(eed_ko_ctrl_diff) - sd(eed_ko_ctrl_diff),
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ),
    max(
      ttest_samples_guides_df_shuffled$diff_eed_ko_score_normalized
    ) + 0.05
  )

n_tests <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

p_be <- ttest_samples_guides_df_shuffled %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  mutate(single_guide = !grepl(";", guide)) %>%
  # color the tope 10 guides with a single guide, not guide combination
  arrange(desc(single_guide), p.value_eed_ko_score_normalized) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(
    x = diff_eed_ko_score_normalized,
    # p-value adjustment on the fly
    y = -log10(p.value_eed_ko_score_normalized * n_tests),
    color = sample,
  )) +
  # make guide combinations a lower opacity
  geom_point(
    aes(alpha = single_guide),
    shape = 19,
    stroke = 0,
    size = 1.8
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  ylab("-log10 adjusted p-value") +
  xlab("Average EED KO score per guide clone") +
  ggrepel::geom_text_repel(
    # label the top10 guides by p-value
    # but only single guides
    data = . %>% mutate(label = ifelse(rank <= 10, guide, "")),
    aes(label = label),
    max.overlaps = 100,
    min.segment.length = 0,
    nudge_x = 0.02,
    force_pull = 0.1,
    force = 1,
    show_guide  = FALSE,
    size = 3
  ) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.2)) +
  guides(alpha = "none") +
  coord_cartesian(xlim = xlim_min_max, ylim = c(0, ylim_max)) +
  ggtitle("EED KO score guides compared to Ctrl HLA+", subtitle = "Guide combinations/clones") +
  scale_color_manual(values = samples_color) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    panel.border = element_blank(),
    axis.line = element_line()
  )

p_ctrl_cas9 <- be_singlets@meta.data %>%
  filter(sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = eed_ko_score_normalized)) +
  stat_density(
    aes(color = sample, fill = sample),
    geom = "area",
    position = "identity",
    alpha = 0.5
  ) +
  coord_cartesian(xlim = xlim_min_max) +
  geom_vline(xintercept = 0.06, linetype = "dotted") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  scale_color_manual(values = unlist(samples_color)) +
  scale_fill_manual(values = unlist(samples_color)) +
  xlab("EED KO score per cell")


p3 <- ggpubr::ggarrange(
  p_be,
  p_ctrl_cas9,
  align = "v",
  heights = c(2, 1),
  ncol = 1,
  nrow = 2,
  common.legend = F
)

pdf(
  file.path(plot_dir, "MF03_eed_ko_score_ttest.pdf"),
  height = 7,
  width = 6
)

p1
p2
p3

dev.off()



####################################################################################################
# plot number of guides per sample that pass the EED KO score t-test filtering

n_tests <- ttest_samples_guides_df %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

p4 <- ttest_samples_guides_df %>%
  filter(diff_eed_ko_score_normalized > 0.06,
         p.value_eed_ko_score_normalized * n_tests < 0.05) %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = sample)) +
  geom_bar(aes(fill = sample)) +
  ylab("# guides") +
  scale_fill_manual(values = samples_color) +
  ggtitle(
    "Number of guides per sample with EED KO score >0.06 and \nsignificantly higher than Ctrl HLA+ sample (pval.adj<0.05) \nIncluding multiple guide integration"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    plot.title = element_text(size = 8),
    panel.border = element_blank(),
    axis.line = element_line()
  )

n_tests <- ttest_samples_guides_single_guide_df %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

# plot number of guides per sample that pass the EED KO score t-test filtering
# with only single guide in cells
p5 <- ttest_samples_guides_single_guide_df %>%
  filter(diff_eed_ko_score_normalized > 0.06,
         p.value_eed_ko_score_normalized * n_tests < 0.05) %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = sample)) +
  geom_bar(aes(fill = sample)) +
  ylab("# guides") +
  scale_fill_manual(values = samples_color) +
  ggtitle(
    "Number of guides per sample with EED KO score >0.06 and \nsignificantly higher than Ctrl HLA+ sample (pval.adj<0.05) \nOnly single guide integration"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    plot.title = element_text(size = 8),
    panel.border = element_blank(),
    axis.line = element_line()
  )

n_tests <- ttest_samples_guides_aggregated_df %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>% nrow()

p6 <- ttest_samples_guides_aggregated_df %>%
  mutate(single_guide = !grepl(";", guide)) %>%
  filter(diff_eed_ko_score_normalized > 0.06,
         p.value_eed_ko_score_normalized * n_tests < 0.05) %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  ggplot(aes(x = sample)) +
  geom_bar(aes(fill = sample, alpha = single_guide)) +
  ylab("# guides") +
  scale_fill_manual(values = samples_color) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  ggtitle(
    "Number of guides per sample with EED KO score >0.06 and \nsignificantly higher than Ctrl HLA+ sample (pval.adj<0.05) \nShaded by single vs. multiple guide integration"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    plot.title = element_text(size = 8),
    panel.border = element_blank(),
    axis.line = element_line()
  )

p7 <- ttest_samples_guides_aggregated_df %>%
  mutate(single_guide = !grepl(";", guide)) %>%
  mutate(n_guides_combination = str_count(guide, ";") + 1) %>%
  filter(diff_eed_ko_score_normalized > 0.06,
         p.value_eed_ko_score_normalized * n_tests < 0.05) %>%
  filter(!sample %in% c("Ctrl-HLA-pos", "CAS9-HLA-pos")) %>%
  group_by(sample, n_guides_combination) %>%
  summarize(n_guides = n()) %>%
  ggplot(aes(x = sample, y = n_guides)) +
  geom_bar(aes(fill = sample, alpha = n_guides_combination), stat = "identity") +
  ylab("# guides") +
  scale_fill_manual(values = samples_color) +
  scale_alpha("Guides in \nclone", range = c(1, 0.1)) +
  ggtitle(
    "Number of guides per sample with EED KO score >0.06 and \nsignificantly higher than Ctrl HLA+ sample (pval.adj<0.05) \nShaded by guides in clone"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    plot.title = element_text(size = 6),
    # Legend title size
    legend.title = element_text(size = 6),
    # Legend label size
    legend.text = element_text(size = 5),
    # Size of legend boxes
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_blank(),
    axis.line = element_line()
  )



pdf(
  file.path(plot_dir, "MF03_eed_ko_score_ttest_n_guides.pdf"),
  height = 4,
  width = 3.5
)

p4
p5
p6
p7

dev.off()
