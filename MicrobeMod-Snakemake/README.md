# MicrobeMod-Snakemake
Snakemake Pipeline for Methylation Detection of Nanopore Sequencing Data using MicrobeMod.  

## Setup
### 1. Install Snakemake and Conda/Mamba  
Install Snakemake and Conda/Mamba following the instructions at this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#:~:text=for%20installing%20Snakemake.-,Installation%20via%20Conda/Mamba,-This%20is%20the). 

### 2. Install dorado
Install the dorado following instructions on their [Github](https://github.com/nanoporetech/dorado). Edit the path to your dorado install in "inputs/config.yaml" on line 4 (dorado: installation path).

### 3. Install MicrobeMod and MicrobeMod's Dependencies 
Install the MicrobeMod and its Dependencies following the instructions on their [Github](https://github.com/cultivarium/MicrobeMod). Edit the path to your MicrobeMod installation in "inputs/config.yaml" on line 12 (MicrobeMod: installation path), and the conda environment name on line 10 (MicrobeMod: conda env name).

### 3.Set up Snakemake Pipeline
#### 1. Edit Experimental Configurations
Edit **config.yaml** file in the `inputs/` Directory:  
    -  [add more]

Create **samples.tsv** file in the `inputs/` Directory: 
  - Create samples.tsv file for your samples with the following required columns (and any other columns for your samples): 
  - Required columns: ["read path", "sample"]
    - "read path": path to .fastq file
    - "sample": unique sample identifier
  - Note: an example samples.tsv files is included in the "inputs/" directory 

#### 2. Edit Resource Specifications 
Edit **config.yaml** file in the `profile/` Directory:
  - Edit number of jobs (samples/processes/rules) to run at once on line 5
  - Edit "default-resources" specifications on line 26
    - Tips: obtain information on comptational resources for config file by: 
      - `sinfo`: shows partitions you have access to, node time limit, and list of node names. 
      - `scontrol show node [node_name]`: displays information about a node in partition. 
        - "CPUTot=[int]": shows how many CPUs (cores) in this node; "cpus_per_task" on line 30 should not exceed this value. 
        - "RealMemory=[int]": shows memory available in thie node; "mem" on line 29 should not exceed this value. 

## Running Pipeline 
### 1. Edit main "run_assembly.sbatch" Script
- Edit name of partition on line 4
- Edit name of conda environment with Snakemake installed on line 10 (if env name is other than "snakemake")

### 2. Submit Script to Cluster
- Submit job to cluster by:  
  ```
  sbatch run_assembly.sbatch
  ```  

## Troubleshooting 
As the pipeline runs, log messages will be saved into file named "main.[slurm_job_ID].err" in the "logs" folder. Here are some tips on debugging: 
- Each rule/job in the pipeline will get its own .err log file. When a rule fails, check the subfolder with the rule's name and sample where it failed. 
- Locked directory: if you get the error message "Error: Directory cannot be locked":
  - Make sure that all jobs from your previous run of the pipeline have completed/cancelled
  - Uncomment line 10 "unlock: True" in `profile/config.yaml` file
  - Run the pipeline following previous step: [Submit Script to Cluster](#2-submit-script-to-cluster)
  - Wait for snakemake to complete running
  - Comment line 10
  - Re-run pipeline, the lock should now be removed  


## Results
**RENAME_LATER.tsv:** tab-separated file with the following columns:  
- `sample`: unique sample name

## Workflow
- Read Mapping: 
- MicrobeMod call_methylation: 