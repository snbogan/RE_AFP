#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=syri_sv_call_redo
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=syri_redo_out/syri_redo_%A_%a.out
#SBATCH --error=syri_redo_err/syri_redo_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --array=1-15

set -euo pipefail

module load minimap2
module load samtools
module load miniconda3
conda activate syri

BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes_redo"
SPECIES_DIR=$(ls -d ${BASE_DIR}/*/ | sed -n "${SLURM_ARRAY_TASK_ID}p")

HAP1=$(find "${SPECIES_DIR}/hap1" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)
HAP2=$(find "${SPECIES_DIR}/hap2" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n 1)

if [ ! -s "$HAP1" ] || [ ! -s "$HAP2" ]; then
    echo "Missing or empty haplotype FASTA in ${SPECIES_DIR}, skipping."
    exit 1
fi

OUTDIR="${SPECIES_DIR}/syri_results"
mkdir -p "$OUTDIR"

# Step 1: Initial alignment
echo "[$(date)] Running minimap2 for chroder..."
minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} "$HAP1" "$HAP2" \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${OUTDIR}/hap2_vs_hap1.bam"
samtools index "${OUTDIR}/hap2_vs_hap1.bam"

# Step 2: chroder
echo "[$(date)] Running chroder..."
chroder -F B -o "${OUTDIR}/reordered" -noref \
    "${OUTDIR}/hap2_vs_hap1.bam" "$HAP1" "$HAP2"

REORDERED_REF="${OUTDIR}/reordered.ref.fasta"
REORDERED_QRY="${OUTDIR}/reordered.qry.fasta"

# Step 2.5: Harmonize chromosomes
echo "[$(date)] Harmonizing chromosomes between ref and qry..."

REF_LIST="${OUTDIR}/ref.contigs"
QRY_LIST="${OUTDIR}/qry.contigs"
COMMON_LIST="${OUTDIR}/common.contigs"

grep "^>" "$REORDERED_REF" | sed 's/^>//' > "$REF_LIST"
grep "^>" "$REORDERED_QRY" | sed 's/^>//' > "$QRY_LIST"

comm -12 <(sort "$REF_LIST") <(sort "$QRY_LIST") > "$COMMON_LIST"

N_COMMON=$(wc -l < "$COMMON_LIST")
if [ "$N_COMMON" -eq 0 ]; then
    echo "ERROR: No shared chromosomes after chroder"
    exit 1
fi

echo "[$(date)] Keeping $N_COMMON shared chromosomes"

# Subset FASTAs
seqkit grep -f "$COMMON_LIST" "$REORDERED_REF" > "${OUTDIR}/ref.sub.fa"
seqkit grep -f "$COMMON_LIST" "$REORDERED_QRY" > "${OUTDIR}/qry.sub.fa"

# Force identical names + order
awk '{print $0 "\tchr" NR}' "$COMMON_LIST" > "${OUTDIR}/rename.map"

seqkit replace -k "${OUTDIR}/rename.map" "${OUTDIR}/ref.sub.fa" \
  | seqkit sort -n > "${OUTDIR}/ref.final.fa"

seqkit replace -k "${OUTDIR}/rename.map" "${OUTDIR}/qry.sub.fa" \
  | seqkit sort -n > "${OUTDIR}/qry.final.fa"

FINAL_REF="${OUTDIR}/ref.final.fa"
FINAL_QRY="${OUTDIR}/qry.final.fa"

# Step 3: Final alignment
echo "[$(date)] Running minimap2 on harmonized assemblies..."
minimap2 -ax asm5 --eqx -t ${SLURM_CPUS_PER_TASK} "$FINAL_REF" "$FINAL_QRY" \
    > "${OUTDIR}/reordered.sam"

# Step 4: SyRI
echo "[$(date)] Running SyRI..."
syri -c "${OUTDIR}/reordered.sam" \
     -r "$FINAL_REF" \
     -q "$FINAL_QRY" \
     -F S \
     -k \
     --dir "$OUTDIR" \
     --prefix syri \
     --nc ${SLURM_CPUS_PER_TASK}

echo "[$(date)] SyRI SV calling complete for ${SPECIES_DIR}"
