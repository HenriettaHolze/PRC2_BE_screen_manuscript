#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name=bartab_guides
#SBATCH --time=00-02:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user='HenriettaHolze@petermac.org'
#SBATCH --mail-type=ALL
#SBATCH --output='logs/%x.%j.out'
#SBATCH --error='logs/%x.%j.err'
#SBATCH --partition=prod_short

#### Load modules ####
# Clean environment
module purge
# load modules
module load nextflow/23.04.1
module load singularity/3.7.3

# running bartab v1.5

export NXF_SINGULARITY_LIBRARYDIR=/dawson_genomics/Projects/bartools_bartab_paper/singularity/

nextflow run DaneVass/bartab -r fbed213c6dabe4357575365a872da88c1692740d \
  -profile singularity \
  -params-file MF03_bartab_guide_params.yaml \
  -w "/scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/MF03_scRNA/bartab_guides/work/" \
  -resume
