#!/usr/bin/env Rscript

# Color: A-G vs C-T
# intensity: percent edited

library(tidyverse)
library(readxl)
library(openxlsx)
plot_dir <- Sys.getenv("PLOT_DIR")

path <- "/dawson_genomics/Projects/PRC2_BE_screen/manuscript/data/be_optimization/DBE_data.xlsx"
excel_sheets(path)

sheets <- openxlsx::getSheetNames(path) 
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=path)

# Merge all experiments and all time points into one dataframe. 
names(data_frame) <- sheets
names(data_frame) <- gsub("_v14", "", names(data_frame))
be_opt <- bind_rows(data_frame, .id = "experiment_day")

# split column with experiment and day information into 2 columns.
be_opt <- be_opt %>%
  separate_wider_delim(cols = experiment_day, delim = "_", names = c("experiment", "day"), cols_remove = F) %>%
  mutate(version = gsub(".*V([0-9]+).*", "\\1", Batch))

# Format Ref allele (original sequence) and Nucleotide (what base is changed into) and position on sgRNA sequence. 
be_opt <- be_opt %>%
  pivot_longer(cols = -c(experiment, day, experiment_day, Batch, Nucleotide, version), names_to = "Ref_allele", values_to = "proportion") %>%
  mutate(Ref_allele = gsub("^(.).*", "\\1", Ref_allele)) %>%
  group_by(experiment_day, Batch, Nucleotide) %>%
  mutate(position = row_number()) %>%
  mutate(Batch = gsub("dual_BE_", "", Batch))

be_opt <- be_opt %>%
  mutate(allowed_edit = ifelse((Ref_allele == "A" & Nucleotide == "G") | (Ref_allele == "C" & Nucleotide == "T"), "on_target", "off_target"))

# order gRNA
be_opt <- be_opt %>%
  mutate(Batch = gsub("DBE_", "", Batch)) %>%
  mutate(sgrna = gsub("^V[0-9]*_?", "", Batch)) %>%
  mutate(sgrna = factor(sgrna, levels = rev(c("", "gRNA_con", "gRNA_F+E", "gRNA_2.1", "gRNA_cr772")), ordered=T))

# Only keep valid changes, ignoring non-changed and off-target effects. 
be_opt_allowed <- be_opt %>% filter(allowed_edit == "on_target")

# Order by BE version.
be_opt_allowed <- be_opt_allowed %>%
  mutate(version = as.numeric(version)) %>%
  arrange(desc(version), sgrna, position)


# Check proportions of facets so that tiles are same height
proportions_tiles <- be_opt_allowed %>%
  ungroup() %>%
  select(experiment, Batch) %>%
  unique() %>%
  group_by(experiment) %>%
  summarize(n = n()) %>% pull(n)

p <- be_opt_allowed %>%
  mutate(Batch = factor(Batch, levels = unique(be_opt_allowed$Batch), ordered = T)) %>%
  ggplot(aes(
    x = position,
    fill = Ref_allele,
    y = Batch,
    alpha = proportion,
  )) +
  geom_tile() +
  facet_grid(experiment ~ day, scales = "free_y") +
  # make squares same height
  ggh4x::force_panelsizes(rows = proportions_tiles) +
  # draw black box
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    colour = "black",
    fill = NA,
    inherit.aes = FALSE
  ) +
  scale_alpha_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.1),
    range = c(0, 1)
  ) +
  coord_cartesian(xlim = c(1, 23)) +
  scale_x_continuous(breaks = seq(1, 23, 1), expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#0010D9", "#B40011")) +
  theme_minimal() +
  theme(panel.grid = element_blank())

pdf(
  file.path(plot_dir, "be_optimization_raster_dbe.pdf"),
  width=15, height=6
)
p
dev.off()

# one per base editor/experiment with each individual scales
experiments <- c("AM5", "AM7", "AM8")

p_experiments <- lapply(experiments, function(selected_experiment) {

be_opt_allowed_experiment <- be_opt_allowed %>% filter(experiment == selected_experiment)
proportions_tiles <- be_opt_allowed_experiment %>%
  ungroup() %>%
  select(Batch) %>%
  unique() %>%
  nrow()

be_opt_allowed_experiment %>%
  mutate(Batch = factor(Batch, levels = unique(be_opt_allowed_experiment$Batch), ordered = T)) %>%
  ggplot(aes(
    x = position,
    fill = Ref_allele,
    y = Batch,
    alpha = proportion,
  )) +
  geom_tile() +
  facet_grid( ~ day) +
  # draw black box
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    colour = "black",
    fill = NA,
    inherit.aes = FALSE
  ) +
  scale_alpha_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.1),
    range = c(0, 1)
  ) +
  coord_cartesian(xlim = c(1, 23)) +
  scale_x_continuous(breaks = seq(1, 23, 1), expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#0010D9", "#B40011")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = proportions_tiles / 23)
})

be_opt_allowed_experiment <- be_opt_allowed %>% filter(experiment == "AM8", version != 6)
proportions_tiles <- be_opt_allowed_experiment %>%
  ungroup() %>%
  select(Batch) %>%
  unique() %>%
  nrow()

p2 <- be_opt_allowed_experiment %>%
  mutate(Batch = factor(Batch, levels = unique(be_opt_allowed_experiment$Batch), ordered = T)) %>%
  ggplot(aes(
    x = position,
    fill = Ref_allele,
    y = Batch,
    alpha = proportion,
  )) +
  geom_tile() +
  facet_grid( ~ day) +
  # draw black box
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    colour = "black",
    fill = NA,
    inherit.aes = FALSE
  ) +
  scale_alpha_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.1),
    range = c(0, 1)
  ) +
  coord_cartesian(xlim = c(1, 23)) +
  scale_x_continuous(breaks = seq(1, 23, 1), expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#0010D9", "#B40011")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        aspect.ratio = proportions_tiles / 23)


pdf(
  file.path(plot_dir, "be_optimization_raster_per_be_dbe.pdf"),
  width=15, height=4
)
p_experiments[[1]]
p_experiments[[2]]
p_experiments[[3]]
p2
dev.off()

