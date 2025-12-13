#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri-hmem
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=04_big_re_mafft
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=04_big_re_mafft_out/04_big_re_mafft_%A_%a.out
#SBATCH --error=04_big_re_mafft_err/04_big_re_mafft_%A_%a.err
#SBATCH --array=1-24
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G

# Load modules
module load miniconda3
conda activate mafft_beast

# Paths
LIST_FILE="/hb/home/snbogan/PolarFish/RE_AFP/big_re_mafft_list_fix.txt"
OUT_DIR="/hb/home/snbogan/PolarFish/RE_AFP/big_re_mafft"

# Create output directory if needed
mkdir -p "$OUT_DIR"

# Get the input fasta for this array task
FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE")

# Safety check
if [[ -z "$FASTA_FILE" ]]; then
    echo "ERROR: No file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Build output filename
BASENAME=$(basename "$FASTA_FILE")
BASENAME=${BASENAME%.fa}
OUT_FILE="${OUT_DIR}/${BASENAME}.mafft.fa"

echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Input: $FASTA_FILE"
echo "Output: $OUT_FILE"

# Run MAFFT
mafft --auto --thread ${SLURM_CPUS_PER_TASK} "$FASTA_FILE" > "$OUT_FILE"
