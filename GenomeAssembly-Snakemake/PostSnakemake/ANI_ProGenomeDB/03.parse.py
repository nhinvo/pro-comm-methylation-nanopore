import pandas as pd 
from pathlib import Path 

ASSEMBLIES_DF = "../ExtractAssemblies/HighQualAssemblies.tsv"
ANI_OUT = 'data/fastANI.tsv'

def parse_ANI(ANI_OUT):
    """
    """
    df = pd.read_table(ANI_OUT, usecols=[0,1,2], names=['query', 'reference', 'ANI'])
    df['ANI'] = df['ANI'].astype(float)

    # obtain genome name from full path 
    for col in ['query', 'reference']:
        df[col] = df[col].str.split('/').str[-1]
        df[col] = df[col].str.replace('.fasta', '')\
            .str.replace('.fa', '')\
            .str.replace('.', '_')

    # filter for top hit 
    df = df.sort_values(by=['query', 'ANI'], ascending=[True, False])
    df = df.groupby('query').head(1)

    # rename cols 
    df = df.rename(columns={
        'query': 'assembly_name', 
        'reference': 'topANI_hit', 
    })


    return df 

def main():
    sdf = pd.read_table(ASSEMBLIES_DF)
    df = parse_ANI(ANI_OUT)

    # merge the 2 & save 
    df = pd.merge(sdf, df, on='assembly_name', how='left')
    df.to_csv('data/parsed_fastANI.tsv', index=False, sep='\t')

main()