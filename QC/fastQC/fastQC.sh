#!/usr/bin/env bash
#SBATCH --job-name=nanoqc
#SBATCH --time 2-0                       
#SBATCH -p sched_mit_chisholm              
#SBATCH -c 20                               
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH -o logs/fastqc.%j.out
#SBATCH -e logs/fastqc.%j.err

fast5_dir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/fast5_pass
fastq_dir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/demux_fastq

outdir=data/fastQC
mkdir -p ${outdir}

# fastqc \
#     --outdir ${outdir} \
#     --threads ${SLURM_CPUS_PER_TASK} \
#     ${fastq_dir}


for file in ${fastq_dir}/*; do
    fastqc \
        --outdir ${outdir} \
        --threads ${SLURM_CPUS_PER_TASK} \
        "$file"
done

echo Done! 

