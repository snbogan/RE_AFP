#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=01_repeatmodeler
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=01_repeatmodeler_out/01_repeatmodeler_%A_%a.out
#SBATCH --error=01_repeatmodeler_err/01_repeatmodeler_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --array=1-22

module load repeatmodeler

cd /hb/home/snbogan/PolarFish/RE_AFP/

# Get list of genome fasta files (*.genomic.fna or *.fa)
GENOMES=($(ls /hb/home/snbogan/PolarFish/RE_AFP/genomes/*/*genomic.fna /hb/home/snbogan/PolarFish/RE_AFP/genomes/*/*.fa 2>/dev/null))

# Pick the genome corresponding to this array task
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID-1]}

# Extract species directory name for output naming
SPECIES=$(basename $(dirname "$GENOME"))

# Make output directory for this species
OUTDIR=repeatmodeler_out/${SPECIES}
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Build RepeatModeler database
BuildDatabase -name ${SPECIES}_db -engine ncbi "$GENOME"

# Run RepeatModeler
RepeatModeler -threads $SLURM_CPUS_PER_TASK -engine ncbi -database ${SPECIES}_db
