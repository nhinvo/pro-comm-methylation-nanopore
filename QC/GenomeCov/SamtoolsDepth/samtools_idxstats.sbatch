#!/usr/bin/env bash
#SBATCH --job-name=samtools_idxstats
#SBATCH --time 1-0                        
#SBATCH --partition sched_mit_chisholm  
#SBATCH --cpus-per-task 20                                
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

eval "$(conda shell.bash hook)"
conda activate samtools

date

bam_dir=/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/read_mapping
out_dir=data/samtools_idxstats
mkdir -p ${out_dir}

for file_path in ${bam_dir}/*.bam; do 
    # obtain file name without extension 
    file_name=$(basename -s "_sorted.bam" "${file_path}")

    echo "Processing file: ${file_name}"
    echo "File path: ${file_path}"

    samtools idxstats \
        --threads ${SLURM_CPUS_PER_TASK} \
        ${file_path} > ${out_dir}/${file_name}.tsv

done

echo Complete! 
date
