#!/bin/bash

# Commands to create sicelore pipeline inputs

# pipeline prefers refFlat file instead of gtf annotation file, see for instructions
# https://github.com/ucagenomix/sicelore-2.1#parameters-2

# install package to convert gtf to refFlat
conda create -c bioconda -n gtftogenepred ucsc-gtftogenepred
conda activate gtftogenepred

# unzip gtf file (same as for Cell Ranger and for epi2me)
cd /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore
zcat /data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/gencode.v32.primary_assembly.annotation.gtf.gz > gencode.v32.primary_assembly.annotation.gtf

# convert gtf to refFlat
gtfToGenePred -genePredExt -geneNameAsName2 gencode.v32.primary_assembly.annotation.gtf gencode.v32.primary_assembly.annotation.refflat.txt
paste <(cut -f 12 gencode.v32.primary_assembly.annotation.refflat.txt) <(cut -f 1-10 gencode.v32.primary_assembly.annotation.refflat.txt) > gencode.v32.refFlat
rm gencode.v32.primary_assembly.annotation.refflat.txt gencode.v32.primary_assembly.annotation.gtf

# pipeline needs BED file of junctions to run minimap2
# from minimap2 documentation:
# --junc-bed FILE
#  	Gene annotations in the BED12 format (aka 12-column BED), or intron positions in 5-column BED. With this option, minimap2 prefers splicing in annotations. 
#     BED12 file can be converted from GTF/GFF3 with ‘paftools.js gff2bed anno.gtf’
# paftools are available in minimap2 github repo
# following installation instructions from https://github.com/lh3/minimap2/blob/master/misc/README.md#introduction
# install paf tools
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8
git clone https://github.com/lh3/minimap2.git
export PATH="$PATH:/home/hholze/bin/"

# convert gtf to BED
cd /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore
zcat /data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/gencode.v32.primary_assembly.annotation.gtf.gz \
  > gencode.v32.primary_assembly.annotation.gtf

/dawson_genomics/Projects/PRC2_BE_screen/scripts/repositories/minimap2/misc/paftools.js gff2bed \
  gencode.v32.primary_assembly.annotation.gtf \
  > gencode.v32.primary_assembly.annotation.bed

rm gencode.v32.primary_assembly.annotation.gtf

# prepare list of cellranger 10X barcodes from short-read data
cd /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore
zcat /dawson_genomics/Projects/PRC2_BE_screen/results/MF03_scRNA/cellranger_results/MF03_P1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz /dawson_genomics/Projects/PRC2_BE_screen/results/MF03_scRNA/cellranger_results/MF03_P2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > MF03_scRNA_cellranger_P1_P2_barcodes.tsv

# Lastly, modify /dawson_genomics/Projects/PRC2_BE_screen/scripts/repositories/sicelore-2.1/Jar/config.xml for 5prime sequencing
# set umi_length to 10 and fileWithAllPossibleTenXbarcodes to 737k-august-2016.txt

