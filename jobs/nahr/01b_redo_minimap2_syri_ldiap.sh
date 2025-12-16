#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=syri_l_diap
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=syri_ldiap.out
#SBATCH --error=syri_ldiap.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB

set -euo pipefail

# -----------------------------
# Environment
# -----------------------------
module load minimap2
module load samtools
module load miniconda3
conda activate syri

# -----------------------------
# Paths (single pair)
# -----------------------------
BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes_redo/l_diapterus"

HAP1_DIR="${BASE_DIR}/hap1"
HAP2_DIR="${BASE_DIR}/hap2"

OUTDIR="${BASE_DIR}/syri_results"
mkdir -p "$OUTDIR"

echo "[$(date)] Running SyRI for l_diapterus"

# -----------------------------
# Input haplotypes
# -----------------------------
HAP1=$(find "$HAP1_DIR" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)
HAP2=$(find "$HAP2_DIR" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)

if [[ ! -s "$HAP1" || ! -s "$HAP2" ]]; then
    echo "ERROR: Missing haplotype FASTA(s)"
    exit 1
fi

echo "HAP1: $HAP1"
echo "HAP2: $HAP2"

# -----------------------------
# Step 1: Initial alignment
# -----------------------------
echo "[$(date)] Initial minimap2 alignment (hap2 vs hap1)..."

minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} \
    "$HAP1" "$HAP2" \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} \
    -o "${OUTDIR}/hap2_vs_hap1.bam"

samtools index "${OUTDIR}/hap2_vs_hap1.bam"

# -----------------------------
# Step 2: chroder
# -----------------------------
echo "[$(date)] Running chroder..."

chroder -F B -o "${OUTDIR}/reordered" -noref \
    "${OUTDIR}/hap2_vs_hap1.bam" "$HAP1" "$HAP2"

REF_IN="${OUTDIR}/reordered.ref.fasta"
QRY_IN="${OUTDIR}/reordered.qry.fasta"

if [[ ! -s "$REF_IN" || ! -s "$QRY_IN" ]]; then
    echo "ERROR: chroder failed"
    exit 1
fi

cd "$OUTDIR"

# -----------------------------
# Step 3: Enforce identical chromosome order
# -----------------------------
echo "[$(date)] Enforcing identical chromosome order (samtools)..."

samtools faidx "$REF_IN"
samtools faidx "$QRY_IN"

cut -f1 "${REF_IN}.fai" > ref.order

REF_FINAL="ref.final.fa"
QRY_FINAL="qry.final.fa"

samtools faidx "$REF_IN" $(cat ref.order) > "$REF_FINAL"
samtools faidx "$QRY_IN" $(cat ref.order) > "$QRY_FINAL"

if [[ "$(grep -c '^>' "$REF_FINAL")" -ne "$(grep -c '^>' "$QRY_FINAL")" ]]; then
    echo "ERROR: FASTAs differ after reordering"
    exit 1
fi

# -----------------------------
# Step 4: Final alignment
# -----------------------------
echo "[$(date)] Final minimap2 alignment..."

minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} \
    "$REF_FINAL" "$QRY_FINAL" \
    > reordered.sam

# -----------------------------
# Step 5: SyRI
# -----------------------------
echo "[$(date)] Running SyRI (--no-chrmatch)..."

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

echo "[$(date)] SyRI completed successfully for l_diapterus"


