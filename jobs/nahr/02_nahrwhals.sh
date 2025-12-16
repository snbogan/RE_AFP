#!/bin/bash
#SBATCH --job-name=nahrwals
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=nahrwhals_out/nahrwhals_%A_%a.out
#SBATCH --error=nahrwhals_err/nahrwhals_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70GB
#SBATCH --array=1-15

# Load conda and activate environment
module load bedtools
module load minimap2
module load r
module load miniconda3
conda activate nahrwhals

## To create .bed for whole ref haplotype, this command was previously run
# cut -f1,2 reordered.ref.fasta.fai | awk '{print $1"\t0\t"$2}' > all_regions.bed

# Base directory containing species subfolders
BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes"

# Get list of species directories
SPECIES_DIR=$(ls -d ${BASE_DIR}/*/ | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Identify hap1 and hap2 fasta files (handles .fa, .fna, .fasta)
HAP1=$(find "${SPECIES_DIR}/hap1" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)
HAP2=$(find "${SPECIES_DIR}/hap2" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)

BED=$(find "${SPECIES_DIR}/hap1" -type f \( -name "*.bed" \) | head -n 1)

# Check that both files exist and are non-empty
if [ ! -s "$HAP1" ] || [ ! -s "$HAP2" ]; then
    echo "Missing or empty haplotype FASTA in ${SPECIES_DIR}, skipping."
    exit 1
fi

# Index the FASTA files before running nahrwals
samtools faidx "$HAP1"
samtools faidx "$HAP2"

# Make output directory
OUTDIR="${SPECIES_DIR}/nahrwals_res"

# Run nahrwals in R
Rscript --vanilla - <<EOF
library(nahrwhals)
nahrwhals(
  ref_fa = "$HAP1",
  asm_fa = "$HAP2",
  regionfile = "$BED",
  outdir = "$OUTDIR",
  threads = 8
)
EOF
