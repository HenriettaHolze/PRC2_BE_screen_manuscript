#!/usr/bin/env Rscript

library(tidyverse)
result_dir <- Sys.getenv("RESULT_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
if (plot_dir != "")
  dir.create(plot_dir, recursive = TRUE)


sample_alt_ref_df_guides_anno_filter_columns <- read_tsv(
  paste0(
    result_dir,
    "MF03/nanopore/MF03_sicelore_detected_edits_guides_vep_guide_clone_annotation_af_filters.tsv"
  ), 
  col_types = cols(alt_allele_frequency = col_double(), 
                   alt_umi_count = col_double(),
                   total_umi_count = col_double(),
                   position_on_guide_label = col_character())
)

filter_columns <- c(
  "edit_not_k562_heterozygous",
  "edit_matches_guide_and_be",
  "edit_3_cells_guide_clone",
  "edit_is_most_frequent",
  "guide_clone_single_guide",
  "guide_high_eed_ko",
  "guide_suz12_no_off_target_nmd",
  "guide_eed_no_off_target_nmd",
  "guide_ezh2_no_off_target_nmd"
)

variant_effect_colors = list(
  "synonymous_variant" = "#E69F00",
  "missense_variant" = "#56B4E9",
  "stop_gained" = "#009E73",
  "3_prime_UTR_variant" = "#F0E442",
  "5_prime_UTR_variant" = "#0072B2",
  "start_lost" = "#CC79A7",
  "splice_region_variant" = "#D55E00"
)

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

##############################
# Apply filtering and plot

# filter the bare minimum, i.e. edits that I'm guide confident were induced by guides and the base editors
sample_alt_ref_df_guides_anno_induced <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(
    edit_matches_guide_and_be,
    edit_3_cells_guide_clone,
    edit_not_k562_heterozygous,
    !is.na(guide_clone)
  )

# filter mutations that we are confident that they are disrupting PRC2
sample_alt_ref_df_guides_anno_functional <- sample_alt_ref_df_guides_anno_filter_columns %>%
  # true across all filtering columns
  filter(if_all(all_of(filter_columns), ~ . == TRUE))

plotSelectedEdits <- function(sample_alt_ref_df_guides_anno_induced,
                              sample_alt_ref_df_guides_anno_functional,
                              label_mutations = F,
                              color = "consequence") {
  data_dummy <-
    data.frame(list(
      aa_position = c(0, 441, 0, 746, 0, 739),
      gene = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12"),
      diff_eed_ko_score_normalized = c(0, 0, 0, 0, 0, 0)
    ))
  
  
  sample_alt_ref_df_guides_anno_plotting <- sample_alt_ref_df_guides_anno_induced %>%
    filter(aa_mutation_1letter != "-") %>%
    merge(
      sample_alt_ref_df_guides_anno_functional %>%
        dplyr::select(POS, sample, guide_clone) %>%
        mutate(is_functional = T),
      by = c("POS", "sample", "guide_clone"),
      all.x = T
    ) %>%
    mutate(is_functional = replace_na(data = is_functional, replace = FALSE)) %>%
    dplyr::select(
      sample,
      editor,
      guide,
      guide_clone,
      aa_position,
      aa_mutation_1letter,
      POS,
      Consequence_single,
      is_functional,
      gene,
      diff_eed_ko_score_normalized
    ) %>%
    unique()
  
  if (color == "consequence") {
    sample_alt_ref_df_guides_anno_plotting <- sample_alt_ref_df_guides_anno_plotting %>%
      mutate(label = ifelse(is_functional, paste0(aa_mutation_1letter, " ", editor), "")) %>%
      mutate(label = gsub("V14", "EC", label))
  } else if (color == "editor") {
    sample_alt_ref_df_guides_anno_plotting <- sample_alt_ref_df_guides_anno_plotting %>%
      mutate(label = ifelse(is_functional, paste0(aa_mutation_1letter), "")) %>%
      mutate(label = gsub("V14", "EC", label))
  }
  
  p1 <- sample_alt_ref_df_guides_anno_plotting %>%
    ggplot(aes(x = aa_position, y = diff_eed_ko_score_normalized)) +
    geom_blank(data = data_dummy) +
    geom_hline(
      yintercept = 0.06,
      linetype = 2,
      # linewidth = 0.5,
      alpha = 0.3
    )
  
  if (color == "consequence") {
    p1 <- p1 +
      geom_point(
        aes(
          alpha = is_functional,
          size = is_functional,
          color = Consequence_single
        ),
        shape = 19,
        stroke = 0
      ) +
      # scale_color_brewer(palette = "Dark2") +
      scale_color_manual(values = variant_effect_colors)
  } else if (color == "editor") {
    p1 <- p1 +
      geom_point(
        aes(
          alpha = is_functional,
          size = is_functional,
          color = sample
        ),
        shape = 19,
        stroke = 0
      ) +
      # scale_color_brewer(palette = "Dark2") +
      scale_color_manual(values = samples_color)
  }
  p1 <- p1 +
    # don't show alpha legend
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(alpha = "none", size = "none") +
    scale_size_manual(values = c(1.5, 2.5))
  
  if (label_mutations) {
    if (color == "consequence") {
      p1 <- p1 +
        ggrepel::geom_text_repel(
          aes(label = label, color = Consequence_single),
          size = 2,
          max.overlaps = Inf,
          min.segment.length = 0
        )
    } else if (color == "editor") {
      p1 <- p1 +
        ggrepel::geom_text_repel(
          aes(label = label, color = sample),
          size = 2,
          max.overlaps = Inf,
          min.segment.length = 0
        )
    }
  }
  
  p1 <- p1 +
    facet_wrap( ~ gene, ncol = 1, scales = "free_x") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "bottom", 
          strip.background = element_blank(),
          legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.5)) +
    ylab("Avevrage EED KO score")
  
  return(p1)
}

p1 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional,
  label_mutations = F
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, most frequent, single guide, high EED KO, no other target NMD")

p2 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional,
  label_mutations = T
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, most frequent, single guide, high EED KO, no other target NMD")



# same as before, but don't remove guide clones that have other target downregulated
sample_alt_ref_df_guides_anno_functional_no_nmd_filter <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(if_all(all_of(c("edit_not_k562_heterozygous",
                         "edit_matches_guide_and_be",
                         "edit_3_cells_guide_clone",
                         "edit_is_most_frequent",
                         "guide_clone_single_guide",
                         "guide_high_eed_ko")), ~ . == TRUE))


p3 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_no_nmd_filter,
  label_mutations = F
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, most frequent, single guide, high EED KO")

p4 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_no_nmd_filter,
  label_mutations = T
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, most frequent, single guide, high EED KO")

sample_alt_ref_df_guides_anno_functional_induced <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(if_all(all_of(c("edit_not_k562_heterozygous",
                         "edit_matches_guide_and_be",
                         "edit_3_cells_guide_clone",
                         "guide_high_eed_ko")), ~ . == TRUE))

# plot without and with labels, colored by consequence and editor
p5 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_induced,
  label_mutations = F
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO")

p6 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_induced,
  label_mutations = T
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO")

p7 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_induced,
  label_mutations = F,
  color = "editor"
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO")

p8 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced,
  sample_alt_ref_df_guides_anno_functional_induced,
  label_mutations = T,
  color = "editor"
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO")

sample_alt_ref_df_guides_anno_functional_induced_non_synonymous <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(if_all(all_of(c("edit_not_k562_heterozygous",
                         "edit_matches_guide_and_be",
                         "edit_3_cells_guide_clone",
                         "guide_high_eed_ko")), ~ . == TRUE)) %>%
  filter(Consequence_single != "synonymous_variant")

sample_alt_ref_df_guides_anno_induced_non_synonymous <- sample_alt_ref_df_guides_anno_filter_columns %>%
  filter(
    edit_matches_guide_and_be,
    edit_3_cells_guide_clone,
    edit_not_k562_heterozygous,
    !is.na(guide_clone)
  ) %>%
  filter(Consequence_single != "synonymous_variant")

p9 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced_non_synonymous,
  sample_alt_ref_df_guides_anno_functional_induced_non_synonymous,
  label_mutations = F
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO, non-synonymous")

p10 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced_non_synonymous,
  sample_alt_ref_df_guides_anno_functional_induced_non_synonymous,
  label_mutations = T
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO, non-synonymous")

p11 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced_non_synonymous,
  sample_alt_ref_df_guides_anno_functional_induced_non_synonymous,
  label_mutations = F,
  color = "editor"
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO, non-synonymous")

p12 <- plotSelectedEdits(
  sample_alt_ref_df_guides_anno_induced_non_synonymous,
  sample_alt_ref_df_guides_anno_functional_induced_non_synonymous,
  label_mutations = T,
  color = "editor"
) + ggtitle(label = element_blank(), subtitle = "K562 <10%, matches guide, 3 cells, high EED KO, non-synonymous")

pdf(
  file.path(
    plot_dir,
    "MF03_eed_ko_score_edits_protein_sequence_stringent_filters.pdf"
  ),
  height = 8.5,
  width = 8
)

p1
p2
p3
p4
p5
p6
p8
p9
p10
p11
p12

dev.off()
