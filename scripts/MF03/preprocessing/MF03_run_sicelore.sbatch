#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --job-name=sicelore
#SBATCH --time=01-00:00:00
#SBATCH --mem=128GB
#SBATCH --mail-user='henrietta.holze@petermac.org'
#SBATCH --mail-type=ALL
#SBATCH --output='logs/%x.%j.out'
#SBATCH --error='logs/%x.%j.err'
#SBATCH --partition=prod_med
#SBATCH --nodelist=stpr-res-compute02

# run sicelore pipeline to create consensus sequences for gene, cell barcode, UMI combinations
# PCR amplification yields >=3 reads per molecule for about 2/3rd of molecules
# consensus sequences will remove sequencing errors for SNV calling

# pipeline needs at least 16GB memory to align to human genome, even for test dataset with few reads.
# pipeline needs >64GB for 10M reads if PCR amplified genes (i.e. many reads at same position)

# pipeline needs to be run from within pipeline directory to find jar files

module load nextflow/23.04.1
echo $PATH

cd /dawson_genomics/Projects/PRC2_BE_screen/scripts/repositories/sicelore-2.1/

tmp_dir="results/MF03/nanopore/sicelore/output/tmp/"
mkdir $tmp_dir

nextflow run sicelore-nf/main.nf \
  -profile conda \
  --fastqdir /pipeline/Runs/Nanopore/20240905_1547_MN22007_FAY33333/no_sample_id/20240905_1547_MN22007_FAY33333_71d723fd/fastq_pass/ \
  --outdir results/MF03/nanopore/sicelore/output \
  --minimapfasta /data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  --refflat /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore/gencode.v32.refFlat \
  --juncbed /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore/gencode.v32.primary_assembly.annotation.bed \
  --cellRangerBCs /scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MF03_nanopore/sicelore/MF03_scRNA_cellranger_P1_P2_barcodes.tsv \
  --fivePbc true \
  --noPolyARequired true \
  --tmpdir $tmp_dir \
  --max_memory '128.GB' --max_cpus 18 \
  -w results/MF03/nanopore/sicelore/work \
  -resume

