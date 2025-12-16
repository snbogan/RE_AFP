#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=02c_repeatmasker_sorted
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=02c_repeatmasker_sorted_out/02c_repeatmasker_sorted_%A_%a.out
#SBATCH --error=02b_repeatmasker_sorted_err/02c_repeatmasker_sorted_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --array=1-14

module load miniconda3
conda activate repeatmasker
module load perl/5.38.0

# Define arrays of genomes and corresponding RepeatModeler libraries
GENOMES=(
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/a_lupus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/a_minor/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/a_flavidus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/c_maculatus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/c_violaceus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_concolor/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_diapterus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_maculatus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_pacificus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/l_platyrhina/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/m_gelatinosum/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/m_pammelas/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/p_gunnellus/syri_results/reordered.ref.fasta"
"/hb/home/snbogan/PolarFish/RE_AFP/nahr/phased_genomes/z_americanus/syri_results/reordered.ref.fasta"
)

LIBRARIES=(
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/a_lupus/RM_3289036.MonAug111436322025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/a_minor/RM_3289034.MonAug111436322025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/a_flavidus/RM_3289017.MonAug111436302025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/c_maculatus/RM_3289061.MonAug111436342025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/c_violaceus/RM_3289027.MonAug111436312025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_concolor/RM_3289071.MonAug111436342025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_diapterus/RM_1891132.MonAug111436332025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_maculatus/RM_1891140.MonAug111436342025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_pacificus/RM_1891136.MonAug111436342025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/l_platyrhina/RM_1891133.MonAug111436332025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/m_gelatinosum/RM_1891134.MonAug111436332025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/m_pammelas/RM_1891127.MonAug111436332025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/p_gunnellus/RM_1891124.MonAug111436322025/consensi.fa.classified"
"/hb/home/snbogan/PolarFish/RE_AFP/repeatmodeler_out/z_americanus/RM_1891172.MonAug111436352025/consensi.fa.classified"
)

# Select genome/library for this task
IDX=$((SLURM_ARRAY_TASK_ID-1))
GENOME=${GENOMES[$IDX]}
LIBRARY=${LIBRARIES[$IDX]}

# Create output directory
SPECIES=$(basename $(dirname $(dirname "$GENOME")))
OUTDIR="/hb/home/snbogan/PolarFish/RE_AFP/repeatmasker_out/${SPECIES}_hap1_reordered_align/"
mkdir -p "$OUTDIR"

# Run RepeatMasker
RepeatMasker \
  -pa $SLURM_CPUS_PER_TASK \
  -lib "$LIBRARY" \
  -dir "$OUTDIR" \
  -a \
  "$GENOME"

