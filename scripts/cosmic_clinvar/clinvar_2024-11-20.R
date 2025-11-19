#!/usr/bin/env Rscript
# DESCRIPTION
# Check mutations in target genes in ClinVar

# Downloaded from `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/` 25/11/24 Henrietta Holze
# Release from 2024-11-20.
#
# `zcat variant_summary.txt.gz | awk -F  $'\t' '$5=="EED"' | grep GRCh38 > variant_summary_eed_hg38.txt`

library(tidyverse)
library(ggpubr)
library(seqinr)
library(patchwork)

con <- gzfile(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/ClinVar/variant_summary.txt.gz",
  "rt"
)
header_clinvar <- readLines(con, n = 1)
close(con)

clinvar_ezh2 <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/ClinVar/variant_summary_ezh2_hg38.txt",
  col_names = str_split(header_clinvar, "\t")[[1]]
)
clinvar_eed <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/ClinVar/variant_summary_eed_hg38.txt",
  col_names = str_split(header_clinvar, "\t")[[1]]
)
clinvar_suz12 <- read_tsv(
  "/dawson_genomics/Projects/PRC2_BE_screen/data/ClinVar/variant_summary_suz12_hg38.txt",
  col_names = str_split(header_clinvar, "\t")[[1]]
)

clinvar <- rbind(clinvar_ezh2, clinvar_suz12, clinvar_eed)

# ClinVar only has aa mutation annotation on the 5aa longer isoform of EZH2.
# Correct the aa position.
# Select missense variants because that's what we're interested in. 
# Stop codons are not informative what aa are important in the complex. 

clinvar_missense_corrected_aa <- clinvar %>%
  # filter for missense variant pattern in mutation name
  filter(grepl(".*p\\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}\\)", Name)) %>%
  # extract the mutation annotation
  mutate(aa_mutation = gsub(".*p\\.(.*)\\)", "\\1", Name)) %>%
  # remove stop codons, only keep missense variants
  filter(!grepl("Ter", aa_mutation)) %>%
  mutate(aa_position = as.numeric(
    gsub("[A-Z][a-z]{2}([0-9]+)[A-Z][a-z]{2}", "\\1", aa_mutation)
  )) %>%
  # create the 1 letter aa mutation code
  mutate(aa1 = a(gsub("[0-9]+.*", "", aa_mutation)), aa2 = a(gsub(".*[0-9]+", "", aa_mutation))) %>%
  mutate(aa2 = ifelse(grepl("Ter", aa_mutation), "*", aa2)) %>%
  mutate(aa2 = replace_na(aa2, "")) %>%
  mutate(aa_mutation_1letter = paste0(aa1, aa_position, aa2)) %>%
  select(-aa2, -aa1) %>%
  # change aa position to correct isoform for EZH2
  mutate(aa_mutation_1letter = ifelse(
    GeneSymbol == "EZH2",
    str_replace(aa_mutation_1letter, "(?<=^[A-Z])(\\d+)(?=[A-Z]$)", function(x) {
      pos <- as.numeric(x)
      ifelse(pos > 297, pos - 5, pos)
    }),
    aa_mutation_1letter
  )) %>%
  mutate(aa_mutation = ifelse(
    GeneSymbol == "EZH2",
    str_replace(aa_mutation, "(?<=^[A-Z][a-z]{2})(\\d+)(?=[A-Z][a-z]{2}$)", function(x) {
      pos <- as.numeric(x)
      ifelse(pos > 297, pos - 5, pos)
    })    ,
    aa_mutation
  )) %>%
  mutate(aa_position = ifelse(
    GeneSymbol == "EZH2" &
      aa_position > 297,
    aa_position - 5,
    aa_position
  )) %>%
  select(
    `#AlleleID`,
    Type,
    GeneSymbol,
    ClinicalSignificance,
    ClinSigSimple,
    aa_position,
    aa_mutation,
    aa_mutation_1letter,
    NumberSubmitters,
    PhenotypeList
  )

# Filter variants for at least 2 submitters and for pathogenic variants. 
# Only 4 variants.

clinvar_missense_corrected_aa %>% as.data.frame() %>%
  filter(ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic"),
         NumberSubmitters > 1) %>%
  select(
    GeneSymbol,
    aa_mutation_1letter,
    NumberSubmitters,
    ClinicalSignificance,
    NumberSubmitters,
    PhenotypeList
  )

# GeneSymbol aa_mutation_1letter NumberSubmitters ClinicalSignificance                                                     PhenotypeList
# 1       EZH2               E740K                3           Pathogenic                             Weaver syndrome|EZH2-related disorder
# 2       EZH2               R679C                8           Pathogenic                Weaver syndrome|not provided|EZH2-related disorder
# 3       EZH2               V621M                4           Pathogenic                             Weaver syndrome|EZH2-related disorder
# 4       EZH2               Y641N                2    Likely pathogenic Malignant melanoma of skin|Lymphoma|Non-Hodgkin lymphoma|Neoplasm

clinvar_missense_corrected_aa %>% as.data.frame() %>%
  filter(ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic")) %>%
  group_by(GeneSymbol, aa_position) %>%
  summarize(NumberSubmitters_aa_position = sum(NumberSubmitters), PhenotypeList = paste0(PhenotypeList, collapse = ";")) %>%
  filter(NumberSubmitters_aa_position > 1)

# there are 9 amino acid positions with multiple submissions.

# GeneSymbol aa_position NumberSubmitters_aa_position PhenotypeList                                                                                                               
# 1 EED                302                            2 Cohen-Gibson syndrome;Cohen-Gibson syndrome                                                                                 
# 2 EZH2               132                            2 Weaver syndrome;not provided                                                                                                
# 3 EZH2               621                            4 Weaver syndrome|EZH2-related disorder                                                                                       
# 4 EZH2               636                            2 Malignant melanoma of skin|Non-Hodgkin lymphoma;Non-Hodgkin lymphoma|Malignant melanoma of skin                             
# 5 EZH2               641                            5 Lymphoma|Non-Hodgkin lymphoma|Malignant melanoma of skin;Non-Hodgkin lymphoma|Malignant melanoma of skin|Lymphoma;Malignant…
# 6 EZH2               679                            8 Weaver syndrome|not provided|EZH2-related disorder                                                                          
# 7 EZH2               724                            2 Weaver syndrome;Weaver syndrome                                                                                             
# 8 EZH2               740                            4 Weaver syndrome|EZH2-related disorder;Weaver syndrome                                                                       
# 9 SUZ12              610                            2 Imagawa-Matsumoto syndrome;Imagawa-Matsumoto syndrome


# Mutations are a bit more reliable if associated to a known disease that is associated with LoF.
# Imagawa-Matsumoto SUZ12
# Weaver EZH2
# Cohen-Gibson EED
clinvar_missense_corrected_aa %>% as.data.frame() %>%
  filter(ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic")) %>%
  filter(grepl("Weaver|Cohen-Gibson|Imagawa-Matsumoto", PhenotypeList)) %>%
  arrange(GeneSymbol, aa_position) %>%
  select(
    GeneSymbol,
    aa_mutation_1letter,
    NumberSubmitters,
    ClinicalSignificance,
    NumberSubmitters,
    PhenotypeList
  )

write_tsv(
  clinvar_missense_corrected_aa %>% as.data.frame() %>%
    filter(ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic")) %>%
    filter(grepl("Weaver|Cohen-Gibson|Imagawa-Matsumoto", PhenotypeList)),
  "/dawson_genomics/Projects/PRC2_BE_screen/results/external_datasets/clinvar/variant_summary_hg38_cosmic_EED_EZH2_SUZ12_likely_path_path_lof_diseases.tsv"
)

# GeneSymbol aa_mutation_1letter NumberSubmitters ClinicalSignificance                                      PhenotypeList
# 1         EED               R236T                1           Pathogenic                              Cohen-Gibson syndrome
# 2         EED               H258Y                1           Pathogenic                              Cohen-Gibson syndrome
# 3         EED               R302S                1           Pathogenic                              Cohen-Gibson syndrome
# 4         EED               R302G                1           Pathogenic                              Cohen-Gibson syndrome
# 5         EED               Y308C                1           Pathogenic                              Cohen-Gibson syndrome
# 6         EED               M366T                1    Likely pathogenic                              Cohen-Gibson syndrome
# 7         EED               A378V                1    Likely pathogenic                              Cohen-Gibson syndrome
# 8        EZH2               H129R                1    Likely pathogenic                                    Weaver syndrome
# 9        EZH2               P132S                1           Pathogenic                                    Weaver syndrome
# 10       EZH2               H158Y                1    Likely pathogenic                                    Weaver syndrome
# 11       EZH2               V621M                4           Pathogenic              Weaver syndrome|EZH2-related disorder
# 12       EZH2               G623S                1    Likely pathogenic                                    Weaver syndrome
# 13       EZH2               D659Y                1    Likely pathogenic                                    Weaver syndrome
# 14       EZH2               M662T                1    Likely pathogenic                                    Weaver syndrome
# 15       EZH2               F667C                1    Likely pathogenic                                    Weaver syndrome
# 16       EZH2               V674L                1    Likely pathogenic                                    Weaver syndrome
# 17       EZH2               T678N                1    Likely pathogenic                                    Weaver syndrome
# 18       EZH2               R679C                8           Pathogenic Weaver syndrome|not provided|EZH2-related disorder
# 19       EZH2               H689Y                1           Pathogenic                                    Weaver syndrome
# 20       EZH2               H706L                1    Likely pathogenic                                    Weaver syndrome
# 21       EZH2               F724L                1    Likely pathogenic                                    Weaver syndrome
# 22       EZH2               F724L                1           Pathogenic                                    Weaver syndrome
# 23       EZH2               A733D                1    Likely pathogenic                                    Weaver syndrome
# 24       EZH2               E740K                3           Pathogenic              Weaver syndrome|EZH2-related disorder
# 25       EZH2               E740D                1           Pathogenic                                    Weaver syndrome
# 26       EZH2               R741G                1    Likely pathogenic                                    Weaver syndrome
# 27      SUZ12               S350F                1    Likely pathogenic                         Imagawa-Matsumoto syndrome
# 28      SUZ12               F603L                1           Pathogenic                         Imagawa-Matsumoto syndrome
# 29      SUZ12               E610V                1           Pathogenic                         Imagawa-Matsumoto syndrome
# 30      SUZ12               E610Q                1           Pathogenic                         Imagawa-Matsumoto syndrome
# 31      SUZ12               L660P                1    Likely pathogenic                         Imagawa-Matsumoto syndrome



clinvar_missense_corrected_aa_plot <- clinvar_missense_corrected_aa %>% as.data.frame() %>%
  filter(ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic")) %>%
  filter(grepl("Weaver|Cohen-Gibson|Imagawa-Matsumoto", PhenotypeList)) %>%
  arrange(GeneSymbol, aa_position) %>%
  select(
    GeneSymbol,
    aa_position,
    NumberSubmitters
  )

data_dummy <-
  data.frame(list(
    aa_position = c(0, 441, 0, 746, 0, 739),
    GeneSymbol = c("EED", "EED", "EZH2", "EZH2", "SUZ12", "SUZ12"),
    NumberSubmitters = c(0, 0, 0, 0, 0, 0)
  ))



p1 <- clinvar_missense_corrected_aa_plot %>%
  filter(GeneSymbol == "EED") %>%
  ggplot(aes(x = aa_position)) +
  geom_histogram(binwidth = 5, aes(weight = NumberSubmitters)) +
  geom_blank(data = data_dummy %>% filter(GeneSymbol == "EED")) +
  ggtitle("EED") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

p2 <- clinvar_missense_corrected_aa_plot %>%
  filter(GeneSymbol == "SUZ12") %>%
  ggplot(aes(x = aa_position)) +
  geom_histogram(binwidth = 5, aes(weight = NumberSubmitters)) +
  geom_blank(data = data_dummy %>% filter(GeneSymbol == "SUZ12")) +
  ggtitle("SUZ12") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

p3 <- clinvar_missense_corrected_aa_plot %>%
  filter(GeneSymbol == "EZH2") %>%
  ggplot(aes(x = aa_position)) +
  geom_histogram(binwidth = 5, aes(weight = NumberSubmitters)) +
  geom_blank(data = data_dummy %>% filter(GeneSymbol == "EZH2")) +
  ggtitle("EZH2") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5)
  )

wrap_plots(p1, p2, p3, ncol = 1)
