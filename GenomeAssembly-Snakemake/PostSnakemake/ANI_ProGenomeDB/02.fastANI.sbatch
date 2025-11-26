#!/usr/bin/env bash
#SBATCH --job-name=fastANI
#SBATCH --time 12:00:00                     
#SBATCH --partition mit_normal  # edit to your partition 
#SBATCH --cpus-per-task 48
#SBATCH --mem 100000      
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

eval "$(conda shell.bash hook)"
conda activate fastani-1.34

echo "$(date): Running fastANI..."

fastANI \
    --threads $SLURM_CPUS_PER_TASK \
    --ql data/query.txt \
    --rl data/reference.txt \
    -o data/fastANI.tsv

echo "$(date): Complete!"
