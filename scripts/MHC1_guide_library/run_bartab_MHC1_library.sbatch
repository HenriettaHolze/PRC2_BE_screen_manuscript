#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name=bartab
#SBATCH --time=00-2:00:00
#SBATCH --mem=32GB
#SBATCH --mail-user='henrietta.holze@petermac.org'
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

nextflow run /researchers/henrietta.holze/splintr_tools/BARtab/BARtab.nf \
  -profile singularity \
  -params-file ./bartab_params_MHC1_18bp.yaml \
  -w "/scratch/teams/dawson_genomics/Projects/PRC2_BE_screen/results/MHC1_guide_library/bartab/work/" \
  -resume

