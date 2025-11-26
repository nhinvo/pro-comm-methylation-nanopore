"""
Purpose: to map each methylated sites to genes and their annotations. 
Note: MicrobeMod already performed coverage filtering (min cov = 10). 

Output: 
    - a .tsv file for each library/file

Nhi N Vo - 06/30/25
"""
import os
import pandas as pd 
from pathlib import Path 

ANNOTATION_DIR = "../../AnvioGeneAnnotation/data/ParsedAnvioAnnot"
METHYLATION_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/microbemod_call_methylation"

def process_data(df, annot_df, fpath):
    """
    """
    # check if the output file already exists
    if os.path.exists(fpath):
        print('Already parsed data! Skipping assembly...')

    # else perform parsing 
    df = df.rename(columns={'Contig': 'contig'})
    df = pd.merge(df, annot_df, how='cross')

    # filter to keep rows where mutation is WITHIN the gene/intergenic regions and has the same chromosome 
    df = df[
        (df['contig_x'] == df['contig_y'])  # methylation and gene on the same contig 
        & (df['Position'] >= df['start']) & (df['Position'] <= df['stop'])
    ]

    # remove 1 of the contig col 
    df = df.drop(columns=['contig_y'])
    df = df.rename(columns={'contig_x': 'contig'})

    # save data as this step is time consuming
    df.to_csv(fpath, sep='\t', index=False)


def process_genome(methylation_fpath, outdir):
    """
    """
    # parse file path 
    genome_name = methylation_fpath.stem.replace('_methylated_sites', '')

    df = pd.read_table(methylation_fpath, usecols=['Contig', 'Position', 'Modification', 'Total_coverage', 'Percent_modified'])

    # obtain path to annotation and parse 
    annot_fath = f'{ANNOTATION_DIR}/{genome_name}.tsv'
    annot_df = pd.read_table(annot_fath)

    # map methylated sites to annotation 
    fpath = f'{outdir}/raw_{genome_name}.tsv'
    process_data(df, annot_df, fpath)


def main():
    outdir = 'data/MethylatedPosAnnot'
    Path(outdir).mkdir(exist_ok=True, parents=True)

    for methylation_fpath in Path(METHYLATION_DIR).glob('*_methylated_sites.tsv'):
        process_genome(methylation_fpath, outdir)

main()