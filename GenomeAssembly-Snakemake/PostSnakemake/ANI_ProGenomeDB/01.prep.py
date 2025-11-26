import pandas as pd 
from pathlib import Path 

ASSEMBLIES_DF = "../ExtractAssemblies/HighQualAssemblies.tsv"
PRO_GENOME_DIR = "/orcd/data/chisholm/002/hstor001/nvo/data/GenomeDB/data/Genomes"

def obtain_query(ASSEMBLIES_DF):
    """
    Creates & saves a list of fasta file paths for Pro assemblies
    from Nanopore reads. 
    """
    df = pd.read_table(ASSEMBLIES_DF)

    df = df[df['genus'].str.contains('Prochlorococcus')]

    path_list = df['bin_fpath'].values.tolist()

    with open('data/query.txt', 'w') as file:
        for path in path_list:
            file.write(f'{path}\n')

def obtain_reference(PRO_GENOME_DIR):
    """
    Creates & saves a list of fasta file paths for Pro genomes
    in the genome database. 
    """
    with open('data/reference.txt', 'w') as file:
        for genome_path in Path(PRO_GENOME_DIR).glob('*'):
            file.write(f'{genome_path}\n')

def main():
    Path('data').mkdir(exist_ok=True)

    obtain_query(ASSEMBLIES_DF)
    obtain_reference(PRO_GENOME_DIR)
    

main()