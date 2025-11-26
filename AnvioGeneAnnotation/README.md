# Anvio Gene Annotation 
This folder contains Annotations (COG20, KEGG, CyCOG6) for each assembled genome used in the MicrobeMod Analysis (assembled from GenomeAssembly-Snakemake folder). 

# Files
## data/GeneCalls Folder
For each genome:  
* *_genecalls.tsv: contains genes that were called using Prodigal  
    * Columns: [gene_callers_id,contig,start,stop,direction,partial,call_type,source,version]  
* *_functions.tsv: contains each gene called and their annotations  
    * Columns: [gene_callers_id source  accession       function        e_value]  

Files in this can be used to map methylation sites to genes and their annotations 

## data/samples.tsv
Came from GenomeAssembly-Snakemake folder: 
* PostSnakemake/PrepMicrobeMod/samples.tsv
* Contains paths to genome assemblies used in MicrobeMod

## data/CorrectedFasta
Contains FASTA files of genome assembies in the correct format for Anvio 

## data/AnvioContigDb
Contains Anvio Contigs Db files - used for annotation. 

## data/CyCogAnnotDir + data/cycog6_product_lookup.tsv
Contains data from annotation using CyCOG6 - used for Anvio annotation. 

# Steps 
Script 01_AnvioGenContigDb.sbatch is used to generate all the data in `data/` folder 