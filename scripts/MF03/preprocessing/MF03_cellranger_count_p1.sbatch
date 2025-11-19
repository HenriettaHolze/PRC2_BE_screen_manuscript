#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --job-name="cellranger_p1"
#SBATCH --time=2-0:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henrietta.holze@petermac.org
#SBATCH --mem=180G
#SBATCH --output="./logs/%x-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/%x-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

module purge
module load cellranger/7.1.0

fastqdir="/pipeline/Runs/NextSeq/240821_VH01624_155_2223HTKNX/ProjectFolders/Project_Maria-Faleeva/Sample_MF03-1-GEX/"
id="MF03_P1"

transcriptome='/data/reference/dawson_labs/genomes/cellranger_reference_GRCh38-2020-A/refdata-gex-GRCh38-2020-A'
localcores=20
localmem=120

cd $RESULT_DIR/MF03/scRNA/cellranger_results

cmd="cellranger count --id ${id} \
    --fastqs ${fastqdir} \
    --transcriptome=${transcriptome} \
    --localcores=${localcores} \
    --localmem=${localmem}"

echo $cmd
$cmd

