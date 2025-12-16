#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=02b_repeatmasker_ldsorted
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=02b_repeatmasker_ldsorted_out/02b_repeatmasker_ldsorted_%A_%a.out
#SBATCH --error=02b_repeatmasker_ldsorted_err/02b_repeatmasker_ldsorted_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB

module load miniconda3
conda activate repeatmasker
module load perl/5.38.0

cd /hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_dearborni/syri_results/

# Get list of genome fasta files (*.genomic.fna or *.fa), excluding p_angeloi and z_viviparus
GENOME="/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_dearborni/syri_results/reordered.ref.fasta"

# Find the consensi.fa.classified file in the unique subdirectory
LIBRARY="/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_dearborni/RM_3289049.MonAug111436332025/consensi.fa.classified"

# Output directory for RepeatMasker
OUTDIR="/hb/home/snbogan/PolarFish/RE_AFP/repeatmasker_out/l_dearborni_hap1_reordered_align/"
mkdir -p "$OUTDIR"

# Run RepeatMasker with species-specific library
RepeatMasker \
  -pa $SLURM_CPUS_PER_TASK \
  -lib "$LIBRARY" \
  -dir "$OUTDIR" \
  -a \
  "$GENOME"
