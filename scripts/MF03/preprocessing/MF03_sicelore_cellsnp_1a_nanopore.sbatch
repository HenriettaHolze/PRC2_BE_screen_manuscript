#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --job-name=sicelore_cellsnp_1a
#SBATCH --time=01-00:00:00
#SBATCH --mem=8GB
#SBATCH --mail-user='henrietta.holze@petermac.org'
#SBATCH --mail-type=ALL
#SBATCH --output='logs/R-%x.%j.out'
#SBATCH --error='logs/R-%x.%j.err'
#SBATCH --partition=prod

## run cellSNP in mode 1a

source activate vireo
ulimit -n 10240
ulimit -v 33170449147

indir_cellsnp="results/MF03/nanopore"
indir_sicelore="results/MF03/nanopore/sicelore/output"

cellsnp-lite \
  -s /scratch/teams$indir_cellsnp/sicelore_molecules.tags.GE_SUZ12_EED_EZH2.bam \
  -b /scratch/teams$indir_cellsnp/barcodes.tsv.gz \
  -O /scratch/teams$indir_cellsnp/sicelore_cellsnp_1a_SUZ12_EED_EZH2/ \
  -R /scratch/teams$indir_cellsnp/sicelore_cellsnp_2b_SUZ12_EED_EZH2/cellSNP.base.vcf.gz \
  -p 30 \
  --minMAF 0.0001 \
  --minCOUNT 10 \
  --gzip \
  --UMItag U8 \
  --cellTAG BC \
  --refseq /data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa

