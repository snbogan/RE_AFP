#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=04_re_mafft
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=04_re_mafft_out/04_re_mafft_%A_%a.out
#SBATCH --error=04_re_mafft_err/04_re_mafft_%A_%a.err
#SBATCH --array=1-104
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load modules
module load miniconda3
conda activate mafft_beast

# Directories
FASTA_DIR="/hb/home/snbogan/PolarFish/RE_AFP/re_fasta"
OUT_DIR="/hb/home/snbogan/PolarFish/RE_AFP/re_mafft"

# Create output directory if needed
mkdir -p "$OUT_DIR"

# Build an indexed list of fasta files
FASTA_LIST=$(ls "$FASTA_DIR"/*.fa | sort)

# Select the fasta for this array task
FASTA_FILE=$(echo "$FASTA_LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Safety check
if [[ -z "$FASTA_FILE" ]]; then
    echo "ERROR: No fasta file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Output filename
BASENAME=$(basename "$FASTA_FILE" .fa)
OUT_FILE="${OUT_DIR}/${BASENAME}.mafft.fa"

echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "Input:  $FASTA_FILE"
echo "Output: $OUT_FILE"

# Run MAFFT
mafft --auto --thread ${SLURM_CPUS_PER_TASK} "$FASTA_FILE" > "$OUT_FILE"

