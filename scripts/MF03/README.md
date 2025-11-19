# MF03 analysis

## Preprocessing

Directory `preprocessing/` contains bash scripts used to run Cell Ranger, Cite-seq-Count, Sicelore, cellsnp-lite and BARtab to preprocess the 10X scRNA-seq, ONT long-read and sgRNA data. 

## R scripts

Set environment variable for all scripts with  
`export RESULT_DIR="/scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/manuscript/results/"`  
or in R
`Sys.setenv(RESULT_DIR="/scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/manuscript/results/")`

Directory for figures
`Sys.setenv(PLOT_DIR="/dawson_genomics/Projects/PRC2_BE_screen/manuscript/figures/")`

Required input to run end-to-end (in this order)

- Filtered count matrices and cell barcode lists of Cell Ranger results (see preprocessing)
- Cite-seq-count LMO count matrix (see preprocessing)
- List of samples and LMOs
- Guide count tables from BARtab results (see preprocessing)
- List of EED KO gene module genes, that are upregulated upon EED KO in K562
    but robust to Interferon stimulation
- sicelore results (gene count matrix) of Nanopore data
- Nanopore variant calling data from cellsnp-lite (see preprocessing)
- VEP results of edits detected by cellsnp-lite

See session info for dependencies. Mainly, `Matrix_1.6-1.1` `SeuratObject_4.1.4` `Seurat_4.3.0`.

