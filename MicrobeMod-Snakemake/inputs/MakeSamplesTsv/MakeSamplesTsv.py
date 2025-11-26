"""
Create samples.tsv file for MicrobeMod. 
    - Replace reference assembly/genome paths to the ones corrected by Anvio. 
"""
import pandas as pd 
from pathlib import Path 

ANVIO_CORRECTED_ASSEMBLY_DIR = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/AnvioGeneAnnotation/data/CorrectedFasta"
ASSEMBLY_DF = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/PostSnakemake/ExtractAssemblies/HighQualAssemblies.tsv"

SAMPLE_LOG_SHEET = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/inputs/NanoporeDNAisolations.xlsx"
BAM_DIR = "/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/demux_modified_basecall"

def obtain_assembly_paths(assembly_dir):
    """
    Returns df with path to Anvio corrected assemblies 
    and their original bam file names. 
    """
    assembly_paths = [str(fpath) for fpath in Path(assembly_dir).glob('*')]
    df = pd.DataFrame({'reference_path': assembly_paths})

    # obtain assembly name
    df['assembly_name'] = df['reference_path'].str.split('/').str[-1]  # assembly file name with extension 
    df['assembly_name'] = df['assembly_name'].str.replace('.fasta', '')  # assembly file name with extension 
    
    return df

def process_dfs(df, assembly_df):
    """
    """
    # data merge
    df = pd.merge(df, assembly_df, on='assembly_name', how='inner')

    # rename cols 
    df = df.rename(columns={
        'bam_filename': 'raw_read_fname', 
        'bam_fpath': 'read_path', 
        'assembly_name': 'sample'
    })

    df.to_csv('../samples.tsv', sep='\t', index=False)


def main():
    df = pd.read_table(ASSEMBLY_DF)
    assembly_df = obtain_assembly_paths(ANVIO_CORRECTED_ASSEMBLY_DIR)

    # merge, process, and save df 
    process_dfs(df, assembly_df)





main()