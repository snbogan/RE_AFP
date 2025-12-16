#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=redo_syri
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=redo_syri_out/redo_syri_%A_%a.out
#SBATCH --error=redo_syri_err/redo_syri_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --array=1-15

set -euo pipefail

module load minimap2
module load samtools
module load miniconda3
conda activate syri

# -----------------------------
# Paths
# -----------------------------
BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes_redo"
SPECIES_DIR=$(ls -d ${BASE_DIR}/*/ | sed -n "${SLURM_ARRAY_TASK_ID}p")
OUTDIR="${SPECIES_DIR}/syri_results"

cd "$OUTDIR"

REF_IN="reordered.ref.fasta"
QRY_IN="reordered.qry.fasta"

if [[ ! -s "$REF_IN" || ! -s "$QRY_IN" ]]; then
    echo "ERROR: Missing reordered FASTAs"
    exit 1
fi

# -----------------------------
# Step 1: Enforce identical chromosome order
# -----------------------------
echo "[$(date)] Enforcing identical chromosome order..."

samtools faidx "$REF_IN"
samtools faidx "$QRY_IN"

cut -f1 "${REF_IN}.fai" > ref.order

REF_FINAL="ref.final.fa"
QRY_FINAL="qry.final.fa"

samtools faidx "$REF_IN" $(cat ref.order) > "$REF_FINAL"
samtools faidx "$QRY_IN" $(cat ref.order) > "$QRY_FINAL"

# sanity check
if [[ "$(grep -c '^>' "$REF_FINAL")" -ne "$(grep -c '^>' "$QRY_FINAL")" ]]; then
    echo "ERROR: FASTAs differ after reordering"
    exit 1
fi

# -----------------------------
# Step 2: Alignment
# -----------------------------
echo "[$(date)] Running minimap2..."

minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} \
    "$REF_FINAL" "$QRY_FINAL" \
    > reordered.sam

# -----------------------------
# Step 3: SyRI (disable chr matching)
# -----------------------------
echo "[$(date)] Running SyRI with --no-chrmatch..."

syri \
  -c reordered.sam \
  -r "$REF_FINAL" \
  -q "$QRY_FINAL" \
  --no-chrmatch \
  -F S \
  -k \
  --nc ${SLURM_CPUS_PER_TASK} \
  --dir "$OUTDIR" \
  --prefix syri

echo "[$(date)] SyRI completed successfully for ${SPECIES_DIR}"
