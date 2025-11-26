#!/usr/bin/env bash
#SBATCH --job-name STREME
#SBATCH --time 12:00:00
#SBATCH --partition sched_mit_hill
#SBATCH --cpus-per-task 1  
#SBATCH --nodes 1  
#SBATCH --ntasks 1  
#SBATCH --mem 5G
#SBATCH --array=2-101%100
#SBATCH -o logs/%a-%j.out
#SBATCH -e logs/%a-%j.err


### PREP ### 
echo "$(date): SLURM Job ID: ${SLURM_JOB_ID}"
echo "$(date): Array index number: ${SLURM_ARRAY_TASK_ID} (out of: ${SLURM_ARRAY_TASK_MAX})"
echo "$(date): Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "$(date): CPU Requested/Allocated for task: ${SLURM_CPUS_PER_TASK}"

eval "$(conda shell.bash hook)"
conda activate meme-5.5.8

# input var/paths
FILENAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/STREME_input.tsv | cut -f1 | tr -d '\r')
POS_PATH=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/STREME_input.tsv | cut -f2 | tr -d '\r')
CONTROL_PATH=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/STREME_input.tsv | cut -f3 | tr -d '\r')

# output var/paths
OUTDIR="data/STREME_Outputs/${FILENAME}"
mkdir -p ${OUTDIR}

# STREME
echo "$(date): Running streme v $(streme --version) on: ${FILENAME}"

streme \
    --dna \
    --p ${POS_PATH} \
    --n ${CONTROL_PATH} \
    --oc ${OUTDIR}

echo "$(date): Competed STREME!"