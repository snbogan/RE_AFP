#!/bin/bash
#SBATCH --job-name=bsig_rm
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=bsig_rm.out
#SBATCH --error=bsig_rm.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --array=1-2

# Load RepeatModeler
module load repeatmodeler

# Set working directory

# List of input FASTA files
# Bsig Lconc Mpamm Apurp

FASTA_FILES=(
    "/hb/home/snbogan/PolarFish/Genome_Assem/Raw_Reads2/b_signatus/"
    "/hb/home/snbogan/PolarFish/Genome_Assem/Raw_Reads2/b_signatus/"
)

# Get the FASTA file for the current task
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Create output directory
OUTPUT_BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$(basename "$FASTA_FILE" .fasta)_RepeatModeler"
mkdir -p "$OUTPUT_DIR"

# Run RepeatModeler
BuildDatabase -name "$OUTPUT_DIR/$(basename "$FASTA_FILE" .fasta)" "$FASTA_FILE"
RepeatModeler -threads 8 -database "$OUTPUT_DIR/$(basename "$FASTA_FILE" .fasta)"



