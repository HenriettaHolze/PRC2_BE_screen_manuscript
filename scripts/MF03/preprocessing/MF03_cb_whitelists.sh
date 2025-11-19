zcat $RESULT_DIR/MF03/scRNA/cellranger_results/MF03_P1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/-1$//g' >  $RESULT_DIR/MF03/scRNA/bartab_guides/whitelists/MF03-1-CRISPR_whitelist.tsv
zcat $RESULT_DIR/MF03/scRNA/cellranger_results/MF03_P2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/-1$//g' >  $RESULT_DIR/MF03/scRNA/bartab_guides/whitelists/MF03-2-CRISPR_whitelist.tsv
