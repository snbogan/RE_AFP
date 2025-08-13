#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=rm2
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=rm2%A_%a.out
#SBATCH --error=rm2%A_%a.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48GB
#SBATCH --array=1-25

# Load RepeatModeler
module load repeatmodeler

# List of input FASTA files
FASTA_FILES=(
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/a_lupus/a_lupus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/a_lupus/a_lupus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/a_minor/a_minor_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/a_minor/a_minor_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_maculatus/l_maculatus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_maculatus/l_maculatus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_platyrhina/l_platyrhina_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_platyrhina/l_platyrhina_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/p_gunnellus/p_gunnellus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/p_gunnellus/p_gunnellus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/z_americanus/z_americanus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/z_americanus/z_americanus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/z_viviparus/Zviv_GCA_040110945.1_Zoavi_1_2_genomic.fa"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/c_maculatus/c_maculatus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/c_maculatus/c_maculatus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/c_violaceous/c_violaceous_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/c_violaceous/c_violaceous_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_diapterus/l_diapterus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_diapterus/l_diapterus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_pacificus/l_pacificus_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_pacificus/l_pacificus_hap2_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/m_gelatinosum/m_gelatinosum_hap1_ctg.fasta"
    "/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/m_gelatinosum/m_gelatinosum_hap2_ctg.fasta"
)

# Get the FASTA file for the current task
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Create output directory
mkdir /hb/home/snbogan/PolarFish/RE_AFP/re_rm_redo

OUTPUT_BASE_DIR="/hb/home/snbogan/PolarFish/RE_AFP/re_rm_redo"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$(basename "$FASTA_FILE" .fasta)_RepeatModeler"
mkdir -p "$OUTPUT_DIR"

# Run RepeatModeler
BuildDatabase -name "$OUTPUT_DIR/$(basename "$FASTA_FILE" .fasta)" "$FASTA_FILE"
RepeatModeler -threads 16 -database "$OUTPUT_DIR/$(basename "$FASTA_FILE" .fasta)"



