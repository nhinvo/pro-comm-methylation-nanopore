#!/usr/bin/env bash
#SBATCH --job-name anvio_pangenomics
#SBATCH --time 5-0
#SBATCH --partition sched_mit_chisholm
#SBATCH --cpus-per-task 1  
#SBATCH --nodes 1  
#SBATCH --ntasks 1  
#SBATCH --mem 12G
#SBATCH --array=2-57%40
#SBATCH -o logs/%a-%j.out
#SBATCH -e logs/%a-%j.err

########################################################################
## SET UP ## 
date

echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array index number: ${SLURM_ARRAY_TASK_ID} (out of: ${SLURM_ARRAY_TASK_MAX})"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "CPU Requested/Allocated for task: ${SLURM_CPUS_PER_TASK}"

eval "$(conda shell.bash hook)"

## INPUT PATHS 
# annotation database dir
KEGG_DIR=/nobackup1b/users/chisholmlab/databases/anvio_kegg
COG_DIR=/nobackup1b/users/chisholmlab/databases/anvio_cogs

# path to file with paths to all genome 
genome_paths=data/HighQualAssemblies.tsv

# path to genome based on SLURM array ID 
genome_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${genome_paths} | cut -f9 | tr -d '\r')
genome_name=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${genome_paths} | cut -f10 | tr -d '\r')

## OUTPUT PATHs
# folder for fasta with correct format
CORRECTED_FASTA_DIR=data/CorrectedFasta  
mkdir -p ${CORRECTED_FASTA_DIR}

# folder for Anvio Contig DB
OUTDIR=data/AnvioContigDb
mkdir -p ${OUTDIR}

# folder for Cycog data 
CYCOGANNOTDIR=data/CyCogAnnotDir
mkdir -p ${CYCOGANNOTDIR}

# folder for extracted gene call data 
GENECALL_DIR=data/GeneCalls
mkdir -p ${GENECALL_DIR}


########################################################################
## 1. CREATE CONTIGS DB AND ANNOTATE ##
conda activate anvio

echo "$(date): Generating Anvio contig database for ${genome_name}"

# if GENOME.db file doesn't already exists, execute command 
if ! [ -e "${OUTDIR}/${genome_name}.db" ]; then
    #REFORMAT GENOME
    anvi-script-reformat-fasta \
        ${genome_path} \
        --simplify-names \
        --output-file ${CORRECTED_FASTA_DIR}/${genome_name}.fasta

    #MAKE CONTIGS DATABASES 
    echo Generating anvio contigs-db for: ${genome_name}...
    anvi-gen-contigs-database \
        --contigs-fasta ${CORRECTED_FASTA_DIR}/${genome_name}.fasta \
        --project-name MesopelagicPangenomics \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --output-db-path ${OUTDIR}/${genome_name}.db

    #ANNOTATE CONTIGS DATABASES WITH HMM SCAN
    echo 'Running HMM scan'
    anvi-run-hmms \
        --contigs-db ${OUTDIR}/${genome_name}.db \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --add-to-functions-table

    #SCAN FOR TRNAS IN CONTIGS DATABASES
    echo 'Scanning for tRNAs'
    anvi-scan-trnas \
        --contigs-db ${OUTDIR}/${genome_name}.db \
        --num-threads ${SLURM_CPUS_PER_TASK}

    #ANNOTATE KEGG
    echo Performing KEGG annotation for: ${genome_name}...
    anvi-run-kegg-kofams \
        --contigs-db ${OUTDIR}/${genome_name}.db \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --kegg-data-dir ${KEGG_DIR}

    #ANNOTATE NCBI COGS
    echo Performing COG annotation for: ${genome_name}...
    anvi-run-ncbi-cogs \
        --contigs-db ${OUTDIR}/${genome_name}.db \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --cog-data-dir ${COG_DIR}

    #EXPORT GENE CALLS AS PROTEINS IN PREPARATION FOR ANNOTATING WITH EXTERNAL CYCOG6 DATABASE
    echo 'Exporting proteins as faa file'
    anvi-get-sequences-for-gene-calls \
        --contigs-db ${OUTDIR}/${genome_name}.db \
        --get-aa-sequences \
        --output-file ${CYCOGANNOTDIR}/${genome_name}_genecalls.faa 
fi

conda deactivate
echo "$(date): Completed Anvio contig-db generation and annotation!"


###################################################################
## 2. CYCOG ANNOTATION ##
conda activate blast-2.16.0

echo 'Annotating CyCOGv6 functions and making table for anvio import'

# if final output file doesn't already exists, execute commands 
if ! [ -e "${CYCOGANNOTDIR}/${genome_name}_CyCOGv6_functions.txt" ]; then
    # PERFORM BLAST 
    echo "$(date): Performing BLASTP on ${genome_name}..."
    blastp \
        -db /orcd/data/chisholm/001/pmberube/databases002/blast/cycog6/CyCOG6 \
        -query ${CYCOGANNOTDIR}/${genome_name}_genecalls.faa \
        -out ${CYCOGANNOTDIR}/${genome_name}_blast.out \
        -evalue 0.000001 \
        -num_threads ${SLURM_CPUS_PER_TASK} \
        -outfmt "6 qseqid sseqid evalue bitscore pident qcovs qcovhsp"

    # FILTER 
    echo "$(date): Generating ${genome_name}_blast_tophit.out ..."
    cat ${CYCOGANNOTDIR}/${genome_name}_blast.out | awk -F "\t" '!seen[$1]++' > ${CYCOGANNOTDIR}/${genome_name}_blast_tophit.out

    echo "$(date): Generating ${genome_name}_functionsfile1.txt ..."
    cat ${CYCOGANNOTDIR}/${genome_name}_blast_tophit.out | \
        awk 'BEGIN {FS="\t"} ; {OFS="\t"} ; {print $1,$2,$3}' | \
        sed -r 's/\t.+[|]/\tCyCOGv6\t/' | \
        awk 'BEGIN {FS="\t"} ; {OFS="\t"} ; {print $1,$2,$3,$4}' > ${CYCOGANNOTDIR}/${genome_name}_functionsfile1.txt

    echo "$(date): Generating ${genome_name}_functionsfile2.txt ..."
    echo 'gene_callers_id source accession function e_value' | sed -r 's/ /\t/g' > ${CYCOGANNOTDIR}/${genome_name}_functionsfile2.txt
    cat ${CYCOGANNOTDIR}/${genome_name}_functionsfile1.txt >> ${CYCOGANNOTDIR}/${genome_name}_functionsfile2.txt

    for CYCOGIID in `cat ${CYCOGANNOTDIR}/${genome_name}_functionsfile1.txt | awk 'BEGIN {FS="\t"} ; {OFS="\t"} ; {print $3}' | sort -u`
    do
    CYCOGPRODUCT=$(grep "${CYCOGIID}" data/cycog6_product_lookup.tsv)
    sed -i -r "s|${CYCOGIID}|${CYCOGPRODUCT}|" ${CYCOGANNOTDIR}/${genome_name}_functionsfile2.txt
    done

    echo "$(date): Generating ${genome_name}_CyCOGv6_functions.txt ... "
    cat ${CYCOGANNOTDIR}/${genome_name}_functionsfile2.txt > ${CYCOGANNOTDIR}/${genome_name}_CyCOGv6_functions.txt
    echo "$(date): Completed Blast step for ${genome_name}!"

fi 

conda deactivate
echo "$(date): Completed all CyCOG Blast and parsing! "


###################################################################
## IMPORT CYCOG FUNCTIONS INTO THE CONTIGS DATABASES ## 
conda activate anvio

echo "$(date): Importing the annotated CyCOGv6 functions into the contigs databases"

anvi-import-functions \
    --contigs-db ${OUTDIR}/${genome_name}.db \
    --input-files ${CYCOGANNOTDIR}/${genome_name}_CyCOGv6_functions.txt


###################################################################
## EXPORT GENE CALLS AND FUNCTIONS ## 
if ! [ -e "${GENECALL_DIR}/${genome_name}_genecalls.tsv" ]; then
    echo "$(date): Exporting gene calls..."

    anvi-export-gene-calls \
        -c ${OUTDIR}/${genome_name}.db \
        --gene-caller prodigal \
        --skip-sequence-reporting \
        -o ${GENECALL_DIR}/${genome_name}_genecalls.tsv
fi 

if ! [ -e "${GENECALL_DIR}/${genome_name}_functions.tsv" ]; then
    echo "$(date): Exporting gene functions..."
  
    anvi-export-functions \
        -c ${OUTDIR}/${genome_name}.db \
        -o ${GENECALL_DIR}/${genome_name}_functions.tsv
fi 


conda deactivate
echo "$(date): Complete!"