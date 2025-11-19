#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --job-name=sicelore_cellsnp_b2
#SBATCH --time=00-20:00:00
#SBATCH --mem=16GB
#SBATCH --mail-user='henrietta.holze@petermac.org'
#SBATCH --mail-type=ALL
#SBATCH --output='logs/%x.%j.out'
#SBATCH --error='logs/%x.%j.err'
#SBATCH --partition=prod

## run cellSNP in mode 2b
# first step of tabulating variants in nanopore scRNA-seq data
# use BAM file with consensus sequences from sicelore
# use human reference hg38 as reference allele

source activate vireo
ulimit -n 10240
ulimit -v 33170449147

## Create BAM file with SUZ12_EED_EZH2 reads from nanopore scRNA-seq data (sicelore)

# Genomic regions
# SUZ12 chr17:31937007-32001038
# EZH2 chr7:148807383-148884291
# EED chr11:86244753-86278810

indir_sicelore="results/MF03/nanopore/sicelore/output"
indir_cellsnp="results/MF03/nanopore"

# get all reads from each capture that map to the correct region and are aligned to SUZ12_EED_EZH2
samtools view -@ 4 -h $indir_sicelore/04b.matrices/molecules.tags.GE.bam \
  chr17:31937007-32001038 chr11:86244753-86278810 chr7:148807383-148884291 |\
  grep -e SUZ12 -e EZH2 -e EED -e '^@' |\
  samtools view -b > /scratch/teams$indir_cellsnp/sicelore_molecules.tags.GE_SUZ12_EED_EZH2.bam

samtools index /scratch/teams$indir_cellsnp/sicelore_molecules.tags.GE_SUZ12_EED_EZH2.bam

# run cellsnp to get a table of positions, ref allele and most frequent alternative allele
# positions with >= 10 UMI and at least 1 in 10k UMI the alternative allele. 
# That means for SUZ12 (154,311 UMI in sicelore output), at least 16 UMI with alternative allele
# important to set minimum quality 

cellsnp-lite \
  -s $indir_cellsnp/sicelore_molecules.tags.GE_SUZ12_EED_EZH2.bam \
  -I MF03_nanopore_sicelore_SUZ12_EED_EZH2 \
  -O $indir_cellsnp/sicelore_cellsnp_2b_SUZ12_EED_EZH2/ \
  --refseq /data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  -p 22 \
  --minMAF 0.0001 \
  --minCOUNT 10 \
  --cellTAG None \
  --UMItag None \
  --gzip \
  --minMAPQ 0

