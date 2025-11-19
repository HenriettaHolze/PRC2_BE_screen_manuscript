#!/usr/bin/env Rscript
# COSMIC zygosity v100
# Extract all missense variants reported in COSMIC in tartet proteins.
# Get all positions with either heterozygous or homozygous mutations reported in COSMIC.

library(tidyverse)
library(ggbreak)
library(patchwork)

# Load COSMIC v100 genomic screen and targeted screen.
# Information on target genes was extracted as follows:

# cd /dawson_genomics/Projects/PRC2_BE_screen/data/cosmic
# zgrep -E -e "^GENE_SYMBOL" -e "^EZH2" -e "^SUZ12" -e "^EED" /data/reference/dawson_labs/COSMIC/v100/Cosmic_CompleteTargetedScreensMutant_v100_GRCh38.tsv.gz |\
#   awk -F'\t' '$27 ~ "y" || $27 ~ "POSITIVE_SCREEN"' \
#   > cosmic_v100_CompleteTargetedScreensMutant_EED_EZH2_SUZ12.tsv
# zgrep -E -e "^GENE_SYMBOL" -e "^EZH2" -e "^SUZ12" -e "^EED" /data/reference/dawson_labs/COSMIC/v100/Cosmic_GenomeScreensMutant_v100_GRCh38.tsv.gz \
#   > cosmic_v100_GenomeScreensMutant_EED_EZH2_SUZ12.tsv

cosmic_db_genomic <- data.table::fread(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/cosmic/cosmic_v100_GenomeScreensMutant_EED_EZH2_SUZ12.tsv",
  sep = "\t"
)
cosmic_db_targeted <- data.table::fread(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/cosmic/cosmic_v100_CompleteTargetedScreensMutant_EED_EZH2_SUZ12.tsv",
  sep = "\t"
)
cosmic_db_samples <- data.table::fread(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/cosmic/COSMIC/v100/Cosmic_Sample_v100_GRCh38.tsv.gz",
  sep = "\t"
)
cosmic_db_classification <- data.table::fread(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/cosmic/COSMIC/v100/Cosmic_Classification_v100_GRCh38.tsv.gz",
  sep = "\t"
)

head(cosmic_db_classification)

# Remove column that indicates whether the mutation was detected or not, the data has already been filtered for positive hits.
cosmic_db_targeted <- cosmic_db_targeted %>% select(-POSITIVE_SCREEN)

# Make sure there is no overlap in sample IDs between the targeted and genome wide screens.
sum(cosmic_db_genomic$COSMIC_SAMPLE_ID %in% cosmic_db_targeted$COSMIC_SAMPLE_ID)

# Merge targeted and genomic screen information.
cosmic_db <- rbind(cosmic_db_genomic, cosmic_db_targeted)

# Merge sample information into the dataframe.
cosmic_db <-
  merge(
    cosmic_db,
    cosmic_db_samples,
    all.x = T,
    by = c("COSMIC_SAMPLE_ID", "SAMPLE_NAME", "COSMIC_PHENOTYPE_ID")
  )

# Merge phenotype information into the dataframe (site, histology)
cosmic_db <-
  merge(
    cosmic_db,
    cosmic_db_classification,
    all.x = T,
    by = c("COSMIC_PHENOTYPE_ID")
  )

# Select the correct isoform that matches the isoform in UniProt for consistent positions.
# EZH2 ENST00000320356.6 == NM_004456.5 == NP_004447.2 == 751aa
# EZH2 ENST00000350995.6 == NM_152998.3 == NP_694543.1 == XP_374578 == 707aa
# EZH2 ENST00000460911.5 == NM_001203247.2 == NP_001190176.1 == 746aa == Q15910
# -> use ENST00000460911.5

# Subset to transcripts that match UniProt.
# This might get rid of some variants that are only present in other isoforms.
cosmic_db <- cosmic_db %>%
  filter(
    TRANSCRIPT_ACCESSION %in% c(
      "ENST00000460911.5",
      "ENST00000322652.9",
      "ENST00000263360.10"
    )
  )

# Select mutations with a single aa missense mutation.
cosmic_db_mis <- cosmic_db %>%
  filter(grepl("missense_variant", MUTATION_DESCRIPTION)) %>%
  mutate(aa_mutation = gsub("^p.", "", MUTATION_AA)) %>%
  mutate(aa_position = gsub("^[A-Z]([0-9]*)[A-Z]$", "\\1", aa_mutation))

# 968 missense mutations reported, across 1,828 samples.
write_tsv(
  cosmic_db_mis,
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense.tsv"
)

cosmic_db_mis <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense.tsv"
)

# subset data to haematopoietic_neoplasm because they are known LoF
cosmic_db_mis_lof <- cosmic_db_mis %>% filter(PRIMARY_HISTOLOGY == "haematopoietic_neoplasm")
write_tsv(
  cosmic_db_mis_lof,
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense_haematopoietic_neoplasm.tsv"
)

cosmic_db_mis_lof <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/cosmic/cosmic_v100_GenomeScreensMutant_CompleteTargetedScreensMutant_EED_EZH2_SUZ12_uniprot_isoforms_missense_haematopoietic_neoplasm.tsv"
)


# Samples include surgery tissue and bone marrow but also cell lines
# table(cosmic_db_mis$SAMPLE_TYPE)

cosmic_db_mis_indiv_counts <- cosmic_db_mis %>%
  ungroup() %>%
  # consider an INDIVIDUAL patient ONLY ONCE
  select(GENE_SYMBOL, GENOMIC_MUTATION_ID, aa_position, INDIVIDUAL_ID, PRIMARY_HISTOLOGY, HISTOLOGY_SUBTYPE_1, HISTOLOGY_SUBTYPE_2) %>%
  unique()

# nrow(cosmic_db_mis_indiv_counts)

table(cosmic_db_mis_indiv_counts$PRIMARY_HISTOLOGY)
# main histologies are 
# haematopoietic_neoplasm (387/1949)
# lymphoid_neoplasm (534/1949)
# malignant_melanoma (139/1949)
# carcinoma (763/1949)

table(cosmic_db_mis_indiv_counts %>% filter(PRIMARY_HISTOLOGY == "carcinoma") %>% pull(HISTOLOGY_SUBTYPE_1))
# among carcinoma, most frequent are 
# adenocarcinoma 254
# squamous_cell_carcinoma 107
# endometrioid_carcinoma 80

cosmic_db_mis_indiv_counts %>%
  group_by(PRIMARY_HISTOLOGY) %>%
  mutate(n_indiv_histo = n()) %>%
  group_by(GENE_SYMBOL, PRIMARY_HISTOLOGY) %>%
  mutate(n_indiv_histo_gene = n(), prop_indiv_gene_histo = n_indiv_histo_gene / n_indiv_histo) %>%
  select(GENE_SYMBOL, PRIMARY_HISTOLOGY, n_indiv_histo, n_indiv_histo_gene, prop_indiv_gene_histo) %>%
  unique() %>%
  arrange(desc(n_indiv_histo), GENE_SYMBOL) %>% head(30) %>%
  as.data.frame()
# both lymphoid_neoplasm and haematopoietic_neoplasm 

cosmic_db_mis_indiv_counts %>%
  mutate(histo = paste0(PRIMARY_HISTOLOGY, "_", HISTOLOGY_SUBTYPE_1)) %>%
  group_by(histo) %>%
  mutate(n_indiv_histo = n()) %>%
  group_by(GENE_SYMBOL, histo) %>%
  mutate(n_indiv_histo_gene = n(), prop_indiv_gene_histo = n_indiv_histo_gene / n_indiv_histo) %>%
  select(GENE_SYMBOL, histo, n_indiv_histo, n_indiv_histo_gene, prop_indiv_gene_histo) %>%
  unique() %>%
  arrange(desc(n_indiv_histo), GENE_SYMBOL) %>% head(30) %>%
  as.data.frame()


# consider an INDIVIDUAL patient ONLY ONCE
cosmic_db_mis_counts <- cosmic_db_mis %>%
  ungroup() %>%
  select(GENE_SYMBOL, GENOMIC_MUTATION_ID, aa_position, INDIVIDUAL_ID) %>%
  unique() %>%
  mutate(aa_position = as.numeric(aa_position))

data_dummy <-
  data.frame(list(
    aa_position = c(0, 441, 0, 746, 0, 739),
    GENE_SYMBOL = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12")
  ))
# cosmic_db_mis_counts %>%
#   ggplot(aes(x = aa_position)) +
#   # geom_bar() +
#   geom_histogram(binwidth = 5) +
#   facet_wrap( ~ GENE_SYMBOL, scales = "free", ncol = 1) +
#   theme_minimal() +
#   theme(panel.grid = element_blank()) +
#   geom_blank(data = data_dummy)

# p_list <- lapply(c("EZH2", "EED", "SUZ12"), function(selected_gene) {
#   p <- cosmic_db_mis %>%
#     filter(GENE_SYMBOL == selected_gene) %>%
#     ungroup() %>%
#     # consider an INDIVIDUAL patient ONLY ONCE
#     select(GENE_SYMBOL, GENOMIC_MUTATION_ID, aa_position, INDIVIDUAL_ID) %>%
#     unique() %>%
#     mutate(aa_position = as.numeric(aa_position)) %>%
#     # group_by(GENE_SYMBOL) %>%
#     # summarise(max(aa_position))
#     ggplot(aes(x = aa_position)) +
#     # geom_bar() +
#     geom_histogram(binwidth = 5) +
#     # facet_wrap( ~ GENE_SYMBOL, scales = "free", ncol = 1) +
#     geom_blank(data = data_dummy %>% filter(GENE_SYMBOL == selected_gene)) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = "bottom",
#       strip.background = element_blank(),
#       panel.border = element_blank(),
#       axis.line = element_line(linewidth = 0.5)
#     )
#   
#   if (selected_gene == "EZH2") {
#     p <- p + scale_y_break(breaks = c(100, 410), scales = 1)
#   }
#   return(p)
# })

plotCosmicDistribution <- function(cosmic_db_mis_counts) {
  p1 <- cosmic_db_mis_counts %>%
    filter(GENE_SYMBOL == "EED") %>%
    ggplot(aes(x = aa_position)) +
    geom_histogram(binwidth = 5) +
    geom_blank(data = data_dummy %>% filter(GENE_SYMBOL == "EED")) +
    ggtitle("EED") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5)
    )
  
  p2 <- cosmic_db_mis_counts %>%
    filter(GENE_SYMBOL == "SUZ12") %>%
    ggplot(aes(x = aa_position)) +
    geom_histogram(binwidth = 5) +
    geom_blank(data = data_dummy %>% filter(GENE_SYMBOL == "SUZ12")) +
    ggtitle("SUZ12") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5)
    )
  
  p3 <- cosmic_db_mis_counts %>%
    filter(GENE_SYMBOL == "EZH2") %>%
    ggplot(aes(x = aa_position)) +
    geom_histogram(binwidth = 5) +
    geom_blank(data = data_dummy %>% filter(GENE_SYMBOL == "EZH2")) +
    ggtitle("EZH2") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5)
    )
  return(list(p1, p2, p3))
}

p_list <- plotCosmicDistribution(cosmic_db_mis_counts)
p1 <- p_list[[1]]
p2 <- p_list[[2]]
p3 <- p_list[[3]]

p3 <- p3 + scale_y_break(breaks = c(100, 375), ticklabels = c(375, 400, 425))


pdf(
  file.path(plot_dir, "cosmic_v100_missense_mutations.pdf"),
  height = 9,
  width = 9
)

wrap_plots(wrap_plots(p1, p2, ncol = 1), p3, ncol = 1, heights = c(2, 1))

dev.off()

# same for haematopoietic mallignancies (LoF)
cosmic_db_mis_lof_counts <- cosmic_db_mis_lof %>%
  ungroup() %>%
  # consider an INDIVIDUAL patient ONLY ONCE
  select(GENE_SYMBOL, GENOMIC_MUTATION_ID, aa_position, INDIVIDUAL_ID) %>%
  unique() %>%
  mutate(aa_position = as.numeric(aa_position))

p_list <- plotCosmicDistribution(cosmic_db_mis_lof_counts)
pp1 <- p_list[[1]]
pp2 <- p_list[[2]]
pp3 <- p_list[[3]]

pdf(
  file.path(plot_dir, "cosmic_v100_missense_mutations_haematopoietic_neoplasm.pdf"),
  height = 9,
  width = 9
)

wrap_plots(pp1, pp2, pp3, ncol = 1)

dev.off()


# same for lymphoid mallignancies (GoF)
cosmic_db_mis_gof <- cosmic_db_mis %>% filter(PRIMARY_HISTOLOGY == "lymphoid_neoplasm")

cosmic_db_mis_gof_counts <- cosmic_db_mis_gof %>%
  ungroup() %>%
  # consider an INDIVIDUAL patient ONLY ONCE
  select(GENE_SYMBOL, GENOMIC_MUTATION_ID, aa_position, INDIVIDUAL_ID) %>%
  unique() %>%
  mutate(aa_position = as.numeric(aa_position))

p_list <- plotCosmicDistribution(cosmic_db_mis_gof_counts)
ppp1 <- p_list[[1]]
ppp2 <- p_list[[2]]
ppp3 <- p_list[[3]]

wrap_plots(ppp1, ppp2, ppp3, ncol = 1)


# # Extract zygosity information on variants.
# # Example of a mutation with conflicting zygosity evidence.
# # Might be related to primary site/cancer subtype.
# cosmic_db_mis %>% filter(aa_mutation== "D317N")
#
# cosmic_db_zygosity <- cosmic_db_mis %>%
#   filter(MUTATION_ZYGOSITY!= "") %>%
#   select(GENE_SYMBOL, aa_position, MUTATION_ZYGOSITY) %>%
#   # # could filter for positions that have at least 2 samples with the same zygosity
#   # # would remove majority of mutations
#   # group_by(`Gene name`, aa_position, `Mutation zygosity`) %>%
#   # mutate(n_evidence = n()) %>%
#   # filter(n_evidence > 1) %>%
#   unique() %>%
#   # remove aa positions with both heterozygous and homozygous examples
#   group_by(GENE_SYMBOL, aa_position) %>%
#   # 1 if only hom or only het, 2 if both are reported
#   mutate(zygosity_detected = n()) %>%
#   ungroup()
#
# # Check which positions are ambiguous. 12 positions.
# cosmic_db_zygosity %>%
#   filter(zygosity_detected == 2) %>%
#   select(-zygosity_detected, -MUTATION_ZYGOSITY) %>%
#   unique() %>%
#   arrange(GENE_SYMBOL, as.numeric(aa_position))
#
# # Only select non-ambiguous positions.
# cosmic_db_zygosity <- cosmic_db_zygosity %>%
#   filter(zygosity_detected == 1) %>%
#   select(-zygosity_detected) %>%
#   arrange(GENE_SYMBOL, as.numeric(aa_position))
#
# # Print positions for each gene to visualize in PyMol
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "EZH2") %>%
#   filter(MUTATION_ZYGOSITY == "hom") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
#
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "EZH2") %>%
#   filter(MUTATION_ZYGOSITY== "het") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
#
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "SUZ12") %>%
#   filter(MUTATION_ZYGOSITY == "hom") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
#
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "SUZ12") %>%
#   filter(MUTATION_ZYGOSITY == "het") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
#
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "EED") %>%
#   filter(MUTATION_ZYGOSITY == "hom") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
#
# cosmic_db_zygosity %>%
#   filter(GENE_SYMBOL == "EED") %>%
#   filter(MUTATION_ZYGOSITY == "het") %>%
#   pull(aa_position) %>%
#   paste0(collapse = "+")
