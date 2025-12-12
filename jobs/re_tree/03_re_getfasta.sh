#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=03_re_getfasta
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=03_re_getfasta.out
#SBATCH --error=03_re_getfasta.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# -----------------------------
# Load modules and environment
# -----------------------------
module load miniconda3
conda activate mafft_beast
module load seqkit
module load bedtools
module load iq-tree

# -----------------------------
# Directories
# -----------------------------
RM_DIR="/hb/home/snbogan/PolarFish/RE_AFP/repeatmasker_out"
GENOME_DIR="/hb/home/snbogan/PolarFish/RE_AFP/genomes"
WORKDIR="/hb/home/snbogan/PolarFish/RE_AFP/iqtree_runs"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# -----------------------------
# 1. Build list of unique RE families
# -----------------------------
FAMILY_LIST="$WORKDIR/all_repeat_families.txt"
rm -f "$FAMILY_LIST"

echo "Collecting RE families..."
for SPECIES_DIR in "$RM_DIR"/*; do
    SPECIES=$(basename "$SPECIES_DIR")
    for OUTFILE in "$SPECIES_DIR"/*.out; do
        awk 'NF>0 && $11 !~ /^#/ && $11 != "repeat" { print $11 }' "$OUTFILE" >> "$FAMILY_LIST"
    done
done

sort -u "$FAMILY_LIST" > tmp && mv tmp "$FAMILY_LIST"
echo "Total RE families found: $(wc -l < "$FAMILY_LIST")"

# -----------------------------
# 2. Loop through RE families
# -----------------------------
while read -r FAMILY; do
    echo "---------------------------------------------"
    echo "Processing RE family: $FAMILY"

    # Make safe name for directories/files
    SAFE_FAMILY="${FAMILY//\//_}"          # replace / with _
    FAMILY_BASE=$(basename "$SAFE_FAMILY")  # subfamily name
    CATEGORY=$(echo "$SAFE_FAMILY" | cut -d'_' -f1)

    FAMILY_DIR="$WORKDIR/$CATEGORY/$FAMILY_BASE"
    mkdir -p "$FAMILY_DIR"
    SEQ_FASTA="$FAMILY_DIR/${FAMILY_BASE}.fa"
    rm -f "$SEQ_FASTA"

    # -----------------------------
    # 2a. Extract sequences for each species
    # -----------------------------
    for SPECIES_DIR in "$RM_DIR"/*; do
        SPECIES=$(basename "$SPECIES_DIR")
        RM_FILE=$(find "$SPECIES_DIR" -maxdepth 1 -name "*.out" | head -n1)
        if [ ! -f "$RM_FILE" ]; then
            echo " ⚠️ RepeatMasker .out file not found for $SPECIES; skipping"
            continue
        fi

        # Genome FASTA
        GENOME_FILE=$(find "$GENOME_DIR/$SPECIES" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | head -n1)
        if [ ! -f "$GENOME_FILE" ]; then
            echo " ⚠️ Genome FASTA not found for $SPECIES; skipping"
            continue
        fi

        BED_FILE="$FAMILY_DIR/${SPECIES}.bed"

        # Extract BED regions for this family (use original FAMILY with slash)
        awk -v fam="$FAMILY" '
            $11 == fam {
                start=$6; end=$7;
                if(start<end) { print $5"\t"start-1"\t"end"\t"fam }
                else { print $5"\t"end-1"\t"start"\t"fam }
            }' "$RM_FILE" > "$BED_FILE"

        # Extract sequences if BED is not empty
        if [[ -s "$BED_FILE" ]]; then
            bedtools getfasta -fi "$GENOME_FILE" -bed "$BED_FILE" -name | sed "s/>/>${SPECIES}_/" >> "$SEQ_FASTA"
        fi
    done

    # Skip if no sequences found
    if [[ ! -s "$SEQ_FASTA" ]]; then
        echo " ⚠️ No sequences found for $FAMILY; skipping"
        continue
    fi

    echo "Sequences extracted for $FAMILY. Stopping before alignment step."
    
done < "$FAMILY_LIST"

echo "All RE families processed."


