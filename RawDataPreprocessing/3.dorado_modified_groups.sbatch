#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --partition=mit_normal_gpu 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --time=1-0
#SBATCH --mem 10G 
#SBATCH --gres=gpu:1
#SBATCH --array=31-45%4 # 11-30%4  # 1-10%4
#SBATCH -o logs/dorado-%a-%j.out
#SBATCH -e logs/dorado-%a-%j.err

echo "$(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Number of GPUs: $SLURM_GPUS"
echo -e "Number of GPUs allocated to the batch step: $SLURM_GPUS_ON_NODE\n"

### 1. dorado basecall prep ### 
# obtain path to list of pod5 tsv file 
pod5_tsv_path=$(realpath data/pod5_groups/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
group_name=$(basename "$pod5_tsv_path")  # e.g. pod5_group_1
group_name="${group_name%.*}"
echo -e "Importing pod5 group list at path: ${pod5_tsv_path}"

# output basecall dir 
basecall_outdir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/pod5_modified_basecall/${group_name}
mkdir -p ${basecall_outdir}

# clear basecall dir (in case leftover from previous runs)
# rm -r ${basecall_outdir}/*

### 2. Run dorado basecall on all files in pod5 list ### 
# total number of lines in file - for tracking progress 
total_lines=$(wc -l < "$pod5_tsv_path")
line_count=0
echo -e "Group ${group_name} has ${total_lines} pod5 files to process\n"

while IFS=$'\t' read -r pod5_fpath fname; do 
    # increment line counter
    line_count=$((line_count + 1))

    echo "__________________________________________________________________________________________"
    echo "$(date)"
    echo "Processing file number $line_count/$total_lines"
    echo "Running dorado basecaller on pod5 file: $fname"
    echo "pod5 file location: $pod5_fpath"
    echo -e "Output .bam file location: ${basecall_outdir}/${fname}.bam\n"

    # run dorado modified basecall 
    dorado/bin/dorado basecaller basecall_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
        ${pod5_fpath} \
        --modified-bases-models basecall_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1,basecall_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1 \
        --kit-name SQK-NBD114-24 \
        --trim all --verbose \
        > ${basecall_outdir}/${fname}.bam

done < ${pod5_tsv_path}

echo Completed dorado basecaller on all files in group ${group_name}!
echo "$(date)"