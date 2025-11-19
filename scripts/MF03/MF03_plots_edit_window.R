#!/usr/bin/env Rscript

# DESCRIPTION
# Plot the edits detected with respect to the position on the guide that was also detected in the cell.
# This indicates the edit window of a base editor.
# TODO
# Make a version with cells with a single guide, because sometimes there are overlapping guides 
# in the same cell, especially in V7 and V14 and SUZ12


library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "") dir.create(plot_dir, recursive = TRUE)


## Load data
sample_alt_ref_df_guides <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides.tsv"
  )
)

## Count kinds of edits at guide positions

# Filter positions that are around guide binding site
n_edits_sample_position_on_guide <- sample_alt_ref_df_guides %>%
  filter(sample != "K562",
         sample != "Ctrl-HLA-pos",
         sample != "CAS9-HLA-pos") %>%
  # filter edits to be remotely close to the guide
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  ungroup()

cat(nrow(n_edits_sample_position_on_guide))
# 11099

n_edits_sample_position_on_guide <- n_edits_sample_position_on_guide %>%
  filter(k562_alt_freq < 0.1 | is.na(k562_alt_freq))

cat(nrow(n_edits_sample_position_on_guide))
# 11079

# count the edits at each position relative to the guide cut site
# for each guide x base editor combination, record which edits have been detected
n_edits_sample_position_on_guide <-
  n_edits_sample_position_on_guide %>%
  dplyr::select(
    sample,
    editor,
    guide,
    REF,
    ALT,
    cell_id,
    POS,
    gene,
    Strand.of.sgRNA,
    position_on_guide,
    position_on_guide_label,
    detected_edit,
    correct_edits
  ) %>%
  unique() %>%
  group_by(
    sample,
    editor,
    Strand.of.sgRNA,
    gene,
    detected_edit,
    position_on_guide,
    position_on_guide_label,
    correct_edits
  ) %>%
  summarize(n_cells = n())

# n_edits_sample_position_on_guide <- n_edits_sample_position_on_guide %>%
#   filter(n_cells >= 3)

allowed_changes <- c("A -> G", "C -> T")
seen_changes <-
  unique(n_edits_sample_position_on_guide$detected_edit)

# order the edits, so that the allowed changes come first, and then all other edits
n_edits_sample_position_on_guide <-
  n_edits_sample_position_on_guide %>%
  mutate(detected_edit = factor(detected_edit, levels = rev(c(
    allowed_changes, seen_changes[!seen_changes %in% allowed_changes]
  )), ordered = T))

axis_label <- c(seq(-20, -1, 1), seq(1, 20, 1), "P", "A", "M", seq(24, 30, 1))

p1 <- n_edits_sample_position_on_guide %>%
  ungroup() %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = detected_edit
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ editor, ncol = 1, scales = "free_y", axes = "all") +
  scale_fill_brewer(palette = "Set3") +
  xlab("Position on guide (1-based)") +
  ylab("# cells") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

p2 <- n_edits_sample_position_on_guide %>%
  ungroup() %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  mutate(edit_color = ifelse(correct_edits, as.character(detected_edit), "off target")) %>%
  mutate(edit_color = factor(edit_color, levels = c("off target", "C -> T", "A -> G"), ordered= T)) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = edit_color
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ editor, ncol = 1, scales = "free_y", axes = "all") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Position on guide (1-based)") +
  ylab("# cells") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

pdf(
  file.path(
    plot_dir, 
    "MF03_nanopore_edit_window.pdf"
  ), 
  height = 9, 
  width = 6
)

p1
p2

dev.off()


p3 <- n_edits_sample_position_on_guide %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = detected_edit
  )) +
  geom_bar(stat = "identity") +
  facet_grid(editor ~ gene + Strand.of.sgRNA,
             scales = "free_y",
             axes = "all") +
  scale_fill_brewer(palette = "Set3") +
  xlab("Position on guide (1-based)") +
  ylab("# cells") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

p4 <- n_edits_sample_position_on_guide %>%
  ungroup() %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  mutate(edit_color = ifelse(correct_edits, as.character(detected_edit), "off target")) %>%
  mutate(edit_color = factor(edit_color, levels = c("off target", "C -> T", "A -> G"), ordered= T)) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = edit_color
  )) +
  geom_bar(stat = "identity") +
  facet_grid(editor ~ gene + Strand.of.sgRNA,
             scales = "free_y",
             axes = "all") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Position on guide (1-based)") +
  ylab("# cells") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

pdf(
  file.path(
    plot_dir, 
    "MF03_nanopore_edit_window_facet_all.pdf"
  ), 
  height = 7, 
  width = 20
)

p3
p4

dev.off()

####################################################################################################
# Number of edits instead of number of cells

# Filter positions that are around guide binding site
n_edits_sample_position_on_guide <- sample_alt_ref_df_guides %>%
  filter(sample != "K562",
         sample != "Ctrl-HLA-pos",
         sample != "CAS9-HLA-pos") %>%
  # filter edits to be remotely close to the guide
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  ungroup()

cat(nrow(n_edits_sample_position_on_guide))
# 11099

n_edits_sample_position_on_guide <- n_edits_sample_position_on_guide %>%
  filter(k562_alt_freq < 0.1 | is.na(k562_alt_freq))

cat(nrow(n_edits_sample_position_on_guide))
# 11079

# count the edits at each position relative to the guide cut site
# for each guide x base editor combination, record which edits have been detected
n_edits_sample_position_on_guide <-
  n_edits_sample_position_on_guide %>%
  dplyr::select(
    sample,
    editor,
    guide,
    REF,
    ALT,
    POS,
    gene,
    Strand.of.sgRNA,
    position_on_guide,
    position_on_guide_label,
    detected_edit,
    correct_edits
  ) %>%
  unique() %>%
  group_by(
    sample,
    editor,
    Strand.of.sgRNA,
    gene,
    detected_edit,
    position_on_guide,
    position_on_guide_label,
    correct_edits
  ) %>%
  summarize(n_cells = n())

n_edits_sample_position_on_guide <- n_edits_sample_position_on_guide %>%
  filter(n_cells >= 3)

allowed_changes <- c("A -> G", "C -> T")
seen_changes <-
  unique(n_edits_sample_position_on_guide$detected_edit)

# order the edits, so that the allowed changes come first, and then all other edits
n_edits_sample_position_on_guide <-
  n_edits_sample_position_on_guide %>%
  mutate(detected_edit = factor(detected_edit, levels = rev(c(
    allowed_changes, seen_changes[!seen_changes %in% allowed_changes]
  )), ordered = T))

axis_label <- c(seq(-20, -1, 1), seq(1, 20, 1), "P", "A", "M", seq(24, 30, 1))

p5 <- n_edits_sample_position_on_guide %>%
  ungroup() %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = detected_edit
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ editor, ncol = 1, scales = "free_y", axes = "all") +
  scale_fill_brewer(palette = "Set3") +
  xlab("Position on guide (1-based)") +
  ylab("# edits") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

p6 <- n_edits_sample_position_on_guide %>%
  ungroup() %>%
  filter(position_on_guide <= 30 & position_on_guide >= -19) %>%
  mutate(edit_color = ifelse(correct_edits, as.character(detected_edit), "off target")) %>%
  mutate(edit_color = factor(edit_color, levels = c("off target", "C -> T", "A -> G"), ordered= T)) %>%
  ggplot(aes(
    x = factor(position_on_guide_label, levels = axis_label),
    y = n_cells,
    fill = edit_color
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ editor, ncol = 1, scales = "free_y", axes = "all") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Position on guide (1-based)") +
  ylab("# edits") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.2, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.2, color = "black"),
    axis.text.x = element_text(size = rel(0.5), angle = 90)
  )

pdf(
  file.path(
    plot_dir, 
    "MF03_nanopore_edit_window_n_edits_3_cells_per_guide.pdf"
  ), 
  height = 9, 
  width = 6
)

p5
p6

dev.off()
