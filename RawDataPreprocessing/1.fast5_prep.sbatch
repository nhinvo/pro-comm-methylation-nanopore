#!/bin/bash
#SBATCH --job-name=pod5
#SBATCH --partition=sched_mit_chisholm
#SBATCH --time 5-0                     
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20 
#SBATCH --mem 250G
#SBATCH -o logs/pod5-%j.out
#SBATCH -e logs/pod5-%j.err

# purpose: to convert all fast5 files to pod5

echo "$(date)"

eval "$(conda shell.bash hook)"
conda activate pod5

fast5_dir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/fast5_pass
cd ${fast5_dir}

### 1. Convert fast5 to pod5 in the same data dir ### 
echo Beginning fast5 to pod5 conversion... 

convert fast5 to pod5 in the same dir 
pod5 convert fast5 \
    *.fast5 \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output . --one-to-one .

echo Completed fast5 to pod5 conversion. 


### 2. move pod5 to different dir ### 
echo Moving pod5 files to separate dir...
cd ../
mkdir -p pod5_pass
mv fast5_pass/*pod5 pod5_pass/

# note: below not needed as pod 5 are split by groups instead
# ### 3. obtain pod5 channel information ###
# echo Obtaining channel information from pod5 files... 
# pod5 view pod5_pass/ --include "read_id, channel" --output pod5_summary.tsv

# ### 4. subset pod5 files - so Dorado runs faster ### 
# echo Splitting pod5 files into different groupb by channels...
# mkdir -p pod5_split_by_channel

# pod5 subset \
#     pod5_pass/ \
#     --threads ${SLURM_CPUS_PER_TASK} \
#     --summary pod5_summary.tsv \
#     --columns channel \
#     --output pod5_split_by_channel

echo Done! 
echo "$(date)"
