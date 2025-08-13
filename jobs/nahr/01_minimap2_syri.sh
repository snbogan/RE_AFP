#!/bin/bash
#SBATCH --job-name=syri_sv_call
#SBATCH --time=2-00:00:00           # 2 days should be enough for 800 Mb genomes
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=syri_out/syri_%A_%a.out
#SBATCH --error=syri_err/syri_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          # use 16 threads for faster alignment
#SBATCH --mem=100GB                 # 100 GB RAM should be sufficient for 800 Mb genome
#SBATCH --array=1-15

# Load necessary modules
module load minimap2
module load samtools
module load miniconda3
conda activate syri          # your conda env with syri, chroder, and dependencies installed

# Base directory containing species subfolders
BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes"

# Get species directory for this task
SPECIES_DIR=$(ls -d ${BASE_DIR}/*/ | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Identify hap1 and hap2 fasta files (handles .fa, .fna, .fasta)
HAP1=$(find "${SPECIES_DIR}/hap1" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)
HAP2=$(find "${SPECIES_DIR}/hap2" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)

# Check that FASTAs exist and are not empty
if [ ! -s "$HAP1" ] || [ ! -s "$HAP2" ]; then
    echo "Missing or empty haplotype FASTA in ${SPECIES_DIR}, skipping."
    exit 1
fi

# Output directory
OUTDIR="${SPECIES_DIR}/syri_results"
mkdir -p "$OUTDIR"

########################################
# Step 1: Initial alignment for chroder
########################################
echo "[$(date)] Running minimap2 to get coords for chroder..."
minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} "$HAP1" "$HAP2" \
    | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${OUTDIR}/hap2_vs_hap1.bam"
samtools index "${OUTDIR}/hap2_vs_hap1.bam"

########################################
# Step 2: Reorder query with chroder
########################################
echo "[$(date)] Running chroder to create pseudo-chromosome hap2..."
chroder -F B -o "${OUTDIR}/reordered" -noref \
    "${OUTDIR}/hap2_vs_hap1.bam" "$HAP1" "$HAP2"

# chroder will output hap2_reordered.fa
REORDERED_QRY="${OUTDIR}/reordered.qry.fasta"
REORDERED_REF="${OUTDIR}/reordered.ref.fasta"

########################################
# Step 3: Final alignment with reordered hap2
########################################
echo "[$(date)] Running minimap2 with reordered hap2..."
minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} "$REORDERED_REF" "$REORDERED_QRY" \
    > "${OUTDIR}/reordered.sam"

########################################
# Step 4: Run SyRI
########################################
echo "[$(date)] Running SyRI..."
syri -c "${OUTDIR}/reordered.sam" \
     -r "$REORDERED_REF" \
     -q "$REORDERED_QRY" \
     -F S \
     -k \
     --dir "$OUTDIR" \
     --prefix syri \
     --nc ${SLURM_CPUS_PER_TASK}

echo "[$(date)] SyRI SV calling complete for ${SPECIES_DIR}"
