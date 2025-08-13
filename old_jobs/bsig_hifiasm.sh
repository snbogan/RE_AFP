#!/bin/bash
#SBATCH --job-name=bsig_hifiasm
#SBATCH --time=7-00:00:00
#SBATCH --partition=512x64-ib
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=bsig_hifiasm.out
#SBATCH --error=bsig_hifiasm.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120GB

# Move to the working directory
cd /hb/home/snbogan/PolarFish/Genome_Assem/Raw_Reads2/b_signatus/

# Load dependencies
module load miniconda3
conda activate hifiasm

# Run the task
hifiasm -o /hb/home/snbogan/PolarFish/Genome_Assem/Raw_Reads2/b_signatus/bsig_hifi -t20 \
/hb/groups/kelley_lab/transfer/eelpout/Original_fastq/Bathymaster_signatus/6085/reads/m64047_220903_024446.reads.fastq.gz \
2> b_signatus.log



