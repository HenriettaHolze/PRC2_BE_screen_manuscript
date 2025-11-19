#!/usr/bin/Rscript

# INPUT
# excel file with MHC1 guide information
# DESCRIPTION
# Read sheet1 of excel file, add column with reverse complement guide, write data to csv.
# Write fasta file with reverse complement guides for BARtab analysis. 20bp and 18bp version
# OUTPUT
# csv and fasta file

library(xlsx)
library(Biostrings)

# read excel file
guide_excel <- "/dawson_genomics/Projects/PRC2_BE_screen/data/mhc1_guides/MHC1_base_editing_library.xlsx"
guides <- read.xlsx(guide_excel, sheetIndex = 1)
# remove lines with NA (excel artifact)
guides <- na.omit(guides)

# read data with additional information per guide
guide_annotation_excel <- "/dawson_genomics/Projects/PRC2_BE_screen/data/mhc1_guides/MHC1-sgrna-designs_custom.xlsx"
guide_annotation <- read.xlsx(guide_annotation_excel, sheetIndex = 1)

# merge both metadata tables
guides <- merge(guides, guide_annotation, by="sgRNA.Sequence")

# read in library count
library_count <- read.csv("/dawson_genomics/Projects/PRC2_BE_screen/data/mhc1_guides/library_count.csv")
colnames(library_count) <- c("sgRNA.Sequence", "library_count")

# merge library count into table
guides <- merge(guides, library_count, by="sgRNA.Sequence")

# get reverse complement guide sequences
rev_comp_guide <- reverseComplement(DNAStringSet(guides$sgRNA.Sequence))
# add as new column to table
guides$sgRNA.Sequence.reverse_complement <- as.character(rev_comp_guide)
readr::write_tsv(guides, "results/MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12.tsv")

# write fasta 20bp reverse compl
seqinr::write.fasta(sequences = as.list(guides$sgRNA.Sequence.reverse_complement), 
                    names = guides$sgRNA_ID, 
                    "results/MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12_20bp.fasta",
                    open = "w", 
                    as.string = T)

# write fasta 18bp reverse compl for scRNA-seq data
seqinr::write.fasta(sequences = as.list(gsub("..$", "", guides$sgRNA.Sequence.reverse_complement)), 
                    names = guides$sgRNA_ID, 
                    "results/MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12_18bp.fasta",
                    open = "w", 
                    as.string = T)

# write fasta for bulk CRISPR screen data, not reverse compl
seqinr::write.fasta(sequences = as.list(guides$sgRNA.Sequence), 
                    names = guides$sgRNA_ID, 
                    "results/MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12_20bp_not_rc.fasta",
                    open = "w", 
                    as.string = T)

# write fasta 18bp for bulk CRISPR screen data, not reverse compl
seqinr::write.fasta(sequences = as.list(gsub("^..", "", guides$sgRNA.Sequence)),
                    names = guides$sgRNA_ID,
                    "results/MHC1_guide_library/MHC1_base_editing_library_EZH2_EED_SUZ12_18bp_not_rc.fasta",
                    open = "w",
                    as.string = T)

