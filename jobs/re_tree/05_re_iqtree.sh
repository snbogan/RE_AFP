#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=05_re_iqtree
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=05_re_iqtree_out/05_re_iqtree_%A_%a.out
#SBATCH --error=05_re_iqtree_err/05_re_iqtree_%A_%a.err
#SBATCH --array=1-81
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

# Load modules
module load iq-tree/2.2.2.6

# Directories
ALIGN_DIR="/hb/home/snbogan/PolarFish/RE_AFP/re_mafft"
OUT_DIR="/hb/home/snbogan/PolarFish/RE_AFP/re_iqtree"

mkdir -p "$OUT_DIR"

# Pick alignment
ALIGNMENT=$(ls ${ALIGN_DIR}/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
GENE=$(basename "$ALIGNMENT" .fa)
OUT_PREFIX="${OUT_DIR}/${GENE}"

echo "=========================================="
echo "Processing TE family: $GENE"
echo "Input alignment: $ALIGNMENT"
echo "Output prefix: $OUT_PREFIX"
echo "=========================================="

# Step 1: Check alignment properties
NTAX=$(grep -c ">" "$ALIGNMENT")
NCHAR=$(grep -v ">" "$ALIGNMENT" | head -n1 | awk '{print length($0)}')
echo "Alignment stats: $NTAX sequences, $NCHAR sites"

# Minimum requirements check
if [ $NTAX -lt 4 ]; then
    echo "ERROR: Too few sequences ($NTAX). Need at least 4 for meaningful tree."
    exit 1
fi

if [ $NCHAR -lt 100 ]; then
    echo "WARNING: Alignment is very short ($NCHAR sites). Results may be unreliable."
fi

# Step 2: Run IQ-TREE to build ML tree
# Using mutation rate: 1.57e-9 substitutions/site/year (from your MrBayes script)
MUTATION_RATE="1.57e-9"
THREADS=$SLURM_CPUS_PER_TASK

echo "Step 1: Building ML tree with IQ-TREE..."
echo "Using mutation rate: $MUTATION_RATE"
echo "Start time: $(date)"

# Run IQ-TREE with model selection and ultrafast bootstrap
iqtree2 -s "$ALIGNMENT" \
    -m GTR+G \
    -bb 1000 \
    -nt $THREADS \
    -alrt 1000 \
    --prefix "$OUT_PREFIX"
