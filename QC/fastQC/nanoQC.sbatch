#!/usr/bin/env bash
#SBATCH --job-name=nanoqc
#SBATCH --time 2-0                       
#SBATCH -p sched_mit_chisholm              
#SBATCH -c 5                               
#SBATCH -N 1
#SBATCH --mem 50G
#SBATCH --array=1-24%10
#SBATCH -o logs/nanoqc.%a.%j.out
#SBATCH -e logs/nanoqc.%a.%j.err

eval "$(conda shell.bash hook)"

raw_main_dir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1
demux_dir=${raw_main_dir}/demux_modified_basecall
bam_fpath=$(realpath ${demux_dir}/SQK*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
fname=$(basename "$bam_fpath" .bam)

conda activate samtools 

fastq_outdir=${raw_main_dir}/demux_fastq
mkdir -p ${fastq_outdir}

samtools bam2fq "${bam_fpath}" | gzip > "${fastq_outdir}/${fname}.fastq.gz"

conda deactivate 

conda activate nanoqc

nanoqc_outdir=data/nanoQC/${fname}
mkdir -p ${nanoqc_outdir}

nanoQC "${fastq_outdir}/${fname}.fastq.gz" --outdir ${nanoqc_outdir}

echo Done! 