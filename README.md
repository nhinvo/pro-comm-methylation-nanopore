# Repository for codes in DNA Modificiation Prochlorococcus and co-occurring Heterotrophs  
This repository contains scripts, data, pipelines, outputs and plots for DNA Modification/Methylation analysis of Prochlorococcus and its co-occurring Heterotrophs.  

## Pre-Processing pf Nanopore Raw Data  
Located in the `RawDataPreprocessing/` directory.  

**Scripts:**   
  - `1.fast5_prep.sbatch`: convert fast5 file from sequencer to pod5 format.     
    - Dependencies/requirements: [POD5](https://github.com/nanoporetech/pod5-file-format) v0.3.34.  
    - Raw file path on Engaging cluster: /orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/fast5_pass  
  - `2.dorado_prep_groups.py`: prepares input for the next step. 
    - Note: this is not required, but splitting inputs into groups allows for parallel processing using SLURM arrays. 
    - Dependencies/requirements: Python Pandas.     
  - `3.dorado_modified_groups.sbatch`: performs basecalling on POD5 files to generate .bam outputs.  
    - Dependencies/requirements: 
      - [Dorado](https://github.com/nanoporetech/dorado) - Oxford Nanopore's Basecaller.  
      - Basecall model files (in `basecall_models/` directory).  
      - Access to HPC GPU.  
  - `4.dorado_demux.sbatch`: demultiplexes file into respective samples/barcodes.  
    - Dependencies/requirements:  
      - [Dorado](https://github.com/nanoporetech/dorado) - Oxford Nanopore's Basecaller.  

## Assemble Genomes from Nanopore long reads   
Located in the `GenomeAssembly-Snakemake/` directory. Reference gesnomes are assembled using a Snakemake pipeline.  

**Assembly Workflow:**  
  1. Filter raw reads by length and mean quality using [filtlong](https://github.com/rrwick/Filtlong).  
     1. Params: `--min_length 1000 --min_mean_q 70`.  
  2. Assemble filtered reads using [flye](https://github.com/mikolmogorov/Flye).
     1. Params: `--meta --nano-corr`.  
  3. Map reads to assembled contigs using [minimap2](https://github.com/lh3/minimap2).  
     1. Params: `-ax map-ont`  
  4. BAM/SAM file conversion and inexing using [samtools](https://github.com/samtools/samtools).  
  5. Contig binning using [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) with default parameters.  
  6. Bin quality assessment using [checkM2](https://github.com/chklovski/CheckM2) with default parameters.  
  7. Bin classification using [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) command `classify_wf` with default parameters.    
  8. Combine all bins and their metadata.  

**Outputs:**  
Bins assembled from the pipeline and their associated metadata are located in the `results/` directory.   
 - `aggregate_table.tsv`: table with all bins assembled and their metadata.  
   - Colunms: [sample, bin_id, bin_contig_num, bin_length, Completeness, Contamination, kingdom, phylum, class, order, family, genus, species, closest_genome_reference, closest_genome_ani, warnings]  
 - `HighQualAssemblies/`: folder containing high quality bins, filtered from pipeline outputs.
   - Generated from Python scripts `PostSnakemake/ExtractAssemblies/ExtractAssemblies.py`.     
   - Thresholds: CheckM2 Completeness > 80 & CheckM2 Contamination < 10.  

Additional QC in `PostSnakemake/ANI_ProGenomeDB/` directory: 
  - Performs fastANI on the assembled Prochlorococcus genomes against previously published Prochlorococcus genomes for quality control.  

**Running The Pipeline:**  
To replicate the assembly workflow, follow the steps below:  
  - Requirements/dependencies: 
    - Install Conda/Mamba 
    - Create conda/mamba environment named "snakemake-7.32.4" with Snakemake v.7.32.4.  
      - `conda create -c conda-forge -c bioconda -c nodefaults -n snakemake-7.32.4 snakemake=7.32.4`
    - Install databases: 
      - checkM2 database  
      - GTDB database  
  - Set-up:  
    - Ensure paths to raw read files in `inputs/samples.tsv` are accessible.  
    - Ensure paths are accessible in `inputs/config.yaml`:  
      - `checkM database`, `GTDB database`, and `scratch directory`.  
    - Edit parameters in `profile/config.yaml` if needed: 
      - Ensure you have access to partitions listed in `profile/config.yaml`.  
      - Ensure resource requests matches with your partition/node.  
      - Edit number of `jobs: ` to run in parallel (submitted to the cluster at once).  
    - Create a `logs/` folder to store SLURM and Snakemake log files.  
  - Run pipeline:  
    - Run pipeline using command: `sbatch run_assembly.sbatch`  
      - Progress will be tracked in the `logs/` directory.  


## Genome Annotation using Anvi'o
Assembled reference genomes were annotated using Anvi'o. Scripts are located in `AnvioGeneAnnotation/` directory.  

**Scripts:**  
  - `01.AnvioGenContigDb.sbatch`: creates Anvi'o contigs DB for genomes from     
  - `02.parse_anvio_genes.py`: 