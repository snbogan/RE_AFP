#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=02_repeatmasker
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=02_repeatmasker_out/02_repeatmasker_%A_%a.out
#SBATCH --error=02_repeatmasker_err/02_repeatmasker_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --array=1-20

module load miniconda3
conda activate repeatmasker
module load perl/5.38.0

cd /hb/home/snbogan/PolarFish/RE_AFP/

# Get list of genome fasta files (*.genomic.fna or *.fa), excluding p_angeloi and z_viviparus
GENOMES=($(ls genomes/*/*genomic.fna genomes/*/*.fa 2>/dev/null \
  | grep -v "/p_angeloi/" \
  | grep -v "/z_viviparus/"))

# Pick the genome corresponding to this array task
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID-1]}

# Extract species directory name
SPECIES=$(basename "$(dirname "$GENOME")")

# Find the consensi.fa.classified file in the unique subdirectory
LIBRARY=$(find "repeatmodeler_out/${SPECIES}" -type f -name "consensi.fa.classified" | head -n 1)

# Safety check: make sure library exists
if [[ ! -f "$LIBRARY" ]]; then
  echo "ERROR: Library file not found for ${SPECIES}" >&2
  exit 1
fi

# Output directory for RepeatMasker
OUTDIR=repeatmasker_out/${SPECIES}
mkdir -p "$OUTDIR"

# Run RepeatMasker with species-specific library
RepeatMasker \
  -pa $SLURM_CPUS_PER_TASK \
  -lib "$LIBRARY" \
  -dir "$OUTDIR" \
  "$GENOME"
