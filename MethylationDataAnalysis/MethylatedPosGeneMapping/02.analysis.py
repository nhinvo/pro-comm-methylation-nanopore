"""
Purpose: to perform different analyses on the raw methylated-gene data, including: 
    - Coverage filtering

"""
import os
import pandas as pd 
from pathlib import Path 

ASSEMBLY_DF = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/PostSnakemake/ExtractAssemblies/HighQualAssemblies.tsv"
PARSED_DIR = "data/MethylatedPosAnnot"
ANNOT_DIR = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/AnvioGeneAnnotation/data/ParsedAnvioAnnot"

def summary_stats(fpath, annot_fpath):
    """
    Returns summary statistics 
    """
    # import 
    df = pd.read_table(fpath)
    annot_df = pd.read_table(annot_fpath)

    # remove duplicated genes/intergenic (a gene can have multiple methylated sites)
    df = df.drop_duplicates(subset=['gene_callers_id'])

    # obtain total gene count 
    gene_count = len(annot_df[annot_df['type'] == 'gene'])
    intergenic_count = len(annot_df[annot_df['type'] == 'intergenic'])

    data = {
        'Methylated_positions': len(df), 
        'Methylated_genes': len(df[df['type'] == 'gene']), 
        'Methylated_intergenic': len(df[df['type'] == 'intergenic']), 
        'Total_gene_count': gene_count, 
        'Total_intergenic_count': intergenic_count, 
    }

    for key, value in data.items():
        data[key] = [value]
    
    df = pd.DataFrame(data)
    df['Percent_gene_methylated'] = (df['Methylated_genes'] / df['Total_gene_count']) * 100 
    df['Percent_intergenic_methylated'] = (df['Methylated_intergenic'] / df['Total_intergenic_count']) * 100 

    return df

def obtain_summary_data(sdf, parsed_dir, annot_dir):
    """
    Obtain and save df with: gene count, methylated genes, methylated intergenic region, 
    pergcent of genes methylated. 
    """
    sdf = sdf[['assembly_name', 'Sample name',	'rep', 'genus']]
    dfs = []

    for fpath in Path(parsed_dir).glob('*.tsv'):
        # obtain path to annot_df for this assembly 
        assembly_name = fpath.stem.replace('raw_', '')
        annot_fpath = f'{annot_dir}/{assembly_name}.tsv'

        df = summary_stats(fpath, annot_fpath)
        df['assembly_name'] = assembly_name

        dfs.append(df)

    df = pd.concat(dfs)
    df = pd.merge(df, sdf, on=['assembly_name'], how='outer')
    df.to_csv('data/MethylationSummary.tsv', sep='\t', index=False)
    df.to_excel('data/MethylationSummary.xlsx', index=False)

def obtain_methylation_gene_pos(sdf, parsed_dir):
    """
    """
    dfs = []
    for fpath in Path(parsed_dir).glob('*.tsv'):
        assembly_name = fpath.stem.replace('raw_', '')

        df = pd.read_table(fpath, usecols=['Position', 'start', 'stop'])
        total_methylated_pos = len(df)

        # scale values to 0-1
        df['range'] = df['stop'] - df['start']
        df['Position_scaled'] = (df['Position'] - df['start']) / df['range']
        df['Position_scaled'] = df['Position_scaled'].round(1)

        df = df[['Position_scaled']]

        df = df.value_counts().reset_index()
        df['percentage'] = (df['count'] / total_methylated_pos) * 100

        df['assembly_name'] = assembly_name

        dfs.append(df)

    df = pd.concat(dfs)

    # sort byassembly &  scaled position 
    df = df.sort_values(['assembly_name', 'Position_scaled'], ascending=[True, True])
    # print(df)

    # long to wide (might not be the best format for plotting)
    # df = df.pivot(
    #     index=['assembly_name'], 
    #     columns=['Position_scaled'], 
    #     values=['percentage']
    # ).reset_index()

    # df.columns = df.columns.droplevel()  # remove multi-index 
    # df = df.rename(columns={'':'assembly_name'})  # rename col 

    # # add metadata 
    df = pd.merge(df, sdf[['assembly_name', 'Sample name', 'rep', 'genus']])

    df.to_csv('data/MethylatedGenePos.tsv', sep='\t', index=False)


def obtain_COG_old(PARSED_DIR):
    """
    """
    outdir = 'data/Methylated_COG'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for fpath in Path(PARSED_DIR).glob('*'):
        assembly_name = fpath.stem.replace('raw_', '')
        df = pd.read_table(fpath)

        df = df.fillna({'COG20_CATEGORY_ACC': 'No Annotation'})
        df['COG20_CATEGORY_ACC'] = df['COG20_CATEGORY_ACC'].str.split('!!!')
        df = df.explode('COG20_CATEGORY_ACC')

        df = df['COG20_CATEGORY_ACC'].value_counts().reset_index()

        df.to_csv(f'{outdir}/{assembly_name}.tsv', sep='\t', index=False)

def obtain_COG(PARSED_DIR, ANNOT_DIR):
    """
    """
    outdir = 'data/Methylated_COG'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for fpath in Path(PARSED_DIR).glob('*'):
        assembly_name = fpath.stem.replace('raw_', '')
        annot_fpath = f'{ANNOT_DIR}/{assembly_name}.tsv'

        df = pd.read_table(fpath)
        annot_df = pd.read_table(annot_fpath)

        dfs = {'methylated': df, 'annot': annot_df}

        # parse both methylation and annot df 
        for name, df in dfs.items():
            df = df[df['type'] == 'gene'].copy()
            df = df.drop_duplicates(subset=['gene_callers_id'])

            # replace NA with "No Annotation"
            df = df.fillna({'COG20_CATEGORY_ACC': 'No Annotation'})

            # split up accessions 
            df['COG20_CATEGORY_ACC'] = df['COG20_CATEGORY_ACC'].str.split('!!!')
            df = df.explode('COG20_CATEGORY_ACC')

            # count accessions 
            df = df['COG20_CATEGORY_ACC'].value_counts().reset_index()

            df = df.rename(columns={'count': f'{name}_count'})

            dfs[name] = df

        # merge methylated COGs and annotated COGs 
        df = pd.merge(dfs['methylated'], dfs['annot'], on=['COG20_CATEGORY_ACC'], how='outer')
        df = df.fillna(0)

        df['methylated_percent'] = (df['methylated_count'] / df['annot_count']) * 100

        df.to_csv(f'{outdir}/{assembly_name}.tsv', sep='\t', index=False)

def obtain_genes_near_methylated_intergenic(PARSED_DIR, ANNOT_DIR):
    """
    """
    outdir = 'data/GenesNearMethylatedIntergenic'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for fpath in Path(PARSED_DIR).glob('*'):
        assembly_name = fpath.stem.replace('raw_', '')
        annot_fpath = f'{ANNOT_DIR}/{assembly_name}.tsv'

        # obtain gene_callers_id before and after intergenic region 
        df = pd.read_table(fpath)
        df = df[df['type'] == 'intergenic']
        df = df.drop_duplicates(subset=['gene_callers_id'])
        df = df[['gene_callers_id']]
        df['gene_callers_id_before'] = df['gene_callers_id'].astype(int)
        df['gene_callers_id_after'] = df['gene_callers_id'].astype(int) + 1
        genes = df['gene_callers_id_before'].values.tolist() + df['gene_callers_id_after'].values.tolist()

        annot_df = pd.read_table(annot_fpath)
        annot_df = annot_df[annot_df['gene_callers_id'].isin(genes)]

        annot_df.to_csv(f'{outdir}/{assembly_name}.tsv', sep='\t', index=False)


def main():
    # sample metadata 
    sdf = pd.read_table(ASSEMBLY_DF)

    # obtain percent gene/intergenic methylated
    # obtain_summary_data(sdf, PARSED_DIR, ANNOT_DIR)
    
    # # obtain methylated position on gene
    # obtain_methylation_gene_pos(sdf, PARSED_DIR)

    # obtain COG categories of genes methylated 
    obtain_COG(PARSED_DIR, ANNOT_DIR)

    # # obtain genes near the methylated intrgenic regions 
    # obtain_genes_near_methylated_intergenic(PARSED_DIR, ANNOT_DIR)



main()