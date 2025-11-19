
# DESCRIPTION

# INPUT
# - Archived EED KO data
#     `/dawson_genomics/Mediaflux/archived_staff_data/Marian_MHC/` and `/dawson_genomics/Mediaflux/archived_projects/Menin/`
# - List of genes in EED KO gene module

library(tidyverse)
plot_dir <- Sys.getenv("PLOT_DIR")
result_dir <- Sys.getenv("RESULT_DIR")

# Christina's EED KO
eed_ko_deg <-
  read.table(
    "/dawson_genomics/Mediaflux/archived_projects/Menin/RNAseq_200409_NB501056_0484_AHJMKYBGXF/analysis/RNAseq_analysis_200409/results/DE/DE_DMSO_EEDKO_vs_DMSO_wt.txt",
    row.names = 1,
    sep = "\t",
    header = T
  )

eed_ko_deg <- rownames_to_column(eed_ko_deg, "ensembl")
eed_ko_deg <- eed_ko_deg %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl, hgnc_symbol))


# Marian's EED KO https://doi.org/10.1016/j.ccell.2019.08.008 Figure 7a.
# Samples: Ctrl sgRNA, Ctrl sgRNA + IFN, EED sgRNA, EED sgRNA + IFN.

# IFN effect in normal and EED KO cells. 
eed_ifn_vs_eed <-
  read.table(
    "/dawson_genomics/Mediaflux/archived_staff_data/Marian_MHC/RNAseq/K562/DE/DE_EED_IFN_vs_EED_noIFN.txt",
    row.names = 1,
    sep = "\t",
    header = T
  )

eed_ifn_vs_eed <- rownames_to_column(eed_ifn_vs_eed, "ensembl")
eed_ifn_vs_eed <- eed_ifn_vs_eed %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl, hgnc_symbol))

ctrl_ifn_vs_ctrl <-
  read.table(
    "/dawson_genomics/Mediaflux/archived_staff_data/Marian_MHC/RNAseq/K562/DE/DE_Safe_IFN_vs_Safe_noIFN.txt",
    row.names = 1,
    sep = "\t",
    header = T
  )

ctrl_ifn_vs_ctrl <- rownames_to_column(ctrl_ifn_vs_ctrl, "ensembl")
ctrl_ifn_vs_ctrl <- ctrl_ifn_vs_ctrl %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl, hgnc_symbol))

# EED KO effect with and without IFN stimulation. 
eed_vs_ctrl <-
  read.table(
    "/dawson_genomics/Mediaflux/archived_staff_data/Marian_MHC/RNAseq/K562/DE/DE_EED_noIFN_vs_Safe_noIFN.txt",
    row.names = 1,
    sep = "\t",
    header = T
  )

eed_vs_ctrl <- rownames_to_column(eed_vs_ctrl, "ensembl")
eed_vs_ctrl <- eed_vs_ctrl %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl, hgnc_symbol))

eed_ifn_vs_ctrl_ifn <-
  read.table(
    "/dawson_genomics/Mediaflux/archived_staff_data/Marian_MHC/RNAseq/K562/DE/DE_EED_IFN_vs_Safe_IFN.txt",
    row.names = 1,
    sep = "\t",
    header = T
  )

eed_ifn_vs_ctrl_ifn <- rownames_to_column(eed_ifn_vs_ctrl_ifn, "ensembl")
eed_ifn_vs_ctrl_ifn <- eed_ifn_vs_ctrl_ifn %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl, hgnc_symbol))

colnames(eed_ifn_vs_ctrl_ifn) <-
  c("ensembl", paste0(colnames(eed_ifn_vs_ctrl_ifn)[-1], "_marian_ifn"))


eed_ko_gene_module_robust <- read_lines(
  paste0(result_dir, "christina_marian_eed_ko_gene_module.txt")
)



p1 <- eed_ko_deg %>%
  mutate(color = ifelse(
    hgnc_symbol %in% eed_ko_gene_module_robust,
    "red",
    "gray"
  )) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 2,
             linetype = 2,
             alpha = 0.5) +
  geom_vline(xintercept = -2,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.5, color = "gray", size = 1) +
  geom_point(data = . %>% filter(color == "red"),
             color = "red", size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Christina's EED KO vs DMSO")

p2 <- eed_vs_ctrl %>%
  mutate(color = ifelse(
    hgnc_symbol %in% eed_ko_gene_module_robust,
    "red",
    "gray"
  )) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 2,
             linetype = 2,
             alpha = 0.5) +
  geom_vline(xintercept = -2,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.5, color = "gray", size = 1) +
  geom_point(data = . %>% filter(color == "red"),
             color = "red", size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Marian's EED KO vs Ctrl")


p3 <- eed_ifn_vs_ctrl_ifn %>%
  mutate(color = ifelse(
    hgnc_symbol_marian_ifn %in% eed_ko_gene_module_robust,
    "red",
    "gray"
  )) %>%
  ggplot(aes(x = log2FoldChange_marian_ifn, y = -log10(padj_marian_ifn))) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 2,
             linetype = 2,
             alpha = 0.5) +
  geom_vline(xintercept = -2,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.5, color = "gray", size = 1) +
  geom_point(data = . %>% filter(color == "red"),
             color = "red", size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Marian's EED KO IFN vs Ctrl IFN")

p4 <- eed_ifn_vs_eed %>%
  mutate(color = ifelse(
    hgnc_symbol %in% eed_ko_gene_module_robust,
    "red",
    "gray"
  )) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 2,
             linetype = 2,
             alpha = 0.5) +
  geom_vline(xintercept = -2,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.5, color = "gray", size = 1) +
  geom_point(data = . %>% filter(color == "red"),
             color = "red", size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Marian's EED KO IFN vs EED KO")

p5 <- ctrl_ifn_vs_ctrl %>%
  mutate(color = ifelse(
    hgnc_symbol %in% eed_ko_gene_module_robust,
    "red",
    "gray"
  )) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 2,
             linetype = 2,
             alpha = 0.5) +
  geom_vline(xintercept = -2,
             linetype = 2,
             alpha = 0.5) +
  geom_point(alpha = 0.5, color = "gray", size = 1) +
  geom_point(data = . %>% filter(color == "red"),
             color = "red", size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Marian's Ctrl IFN vs Ctrl")

# Plot the gene signature in the DEG results of the published data
pdf(file.path(plot_dir, "eed_ko_gene_module_volcano.pdf"), width = 12, height = 8)

gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 3)

dev.off()

# save as png in high resolution because pdf is huge
png(
  file.path(plot_dir, "eed_ko_gene_module_volcano.png"),
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 3)

dev.off()



