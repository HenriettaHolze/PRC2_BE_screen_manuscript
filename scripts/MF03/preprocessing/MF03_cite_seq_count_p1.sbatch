#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --ntasks=12
#SBATCH --job-name=cite-seq-count-p1
#SBATCH --time=00-4:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=henrietta.holze@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --output='./logs/CITE-seq_count-pool1.%j.out'
#SBATCH --error='./logs/CITE-seq_count-pool1.%j.err'
#SBATCH --partition=prod_med

#### Load default modules ####
module purge

# load CITE-seq-Count-1.4.5 
source activate cite-seq-count

#### RUN Parameters: ####

data_dir="/pipeline/Runs/NextSeq/240821_VH01624_155_2223HTKNX/ProjectFolders/Project_Maria-Faleeva/Sample_MF03-1-MULTIseq"
READ1="$data_dir/MF03-1-MULTIseq_S5_R1_001.fastq.gz"
READ2="$data_dir/MF03-1-MULTIseq_S5_R2_001.fastq.gz"
SAMPLE="MF03_P1"
WL="$RESULT_DIR/MF03/scRNA/bartab_guides/whitelists/MF03-1-CRISPR_whitelist.tsv"

TAGS="$RESULT_DIR/MF03/scRNA/lmo_tags_MF03.csv"
OUT="$RESULT_DIR/MF03/scRNA/lmo_results/${SAMPLE}"
UNKNOWNS="unknowns.txt"
CELLS=23000

CITE-seq-Count -T 12 -R1 $READ1 -R2 $READ2 -t $TAGS -o $OUT \
    -u $UNKNOWNS --dense -cells $CELLS \
    -cbf 1 -cbl 16 -umif 17 -umil 26 \
    --max-errors 2 -wl $WL
