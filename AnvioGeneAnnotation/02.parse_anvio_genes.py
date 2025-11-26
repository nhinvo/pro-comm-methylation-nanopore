"""
Combine gene call data (coordinates) with gene function data. 
    - Also add intergenic regions. 
"""
import pandas as pd 
from pathlib import Path 

GENE_CALL_DIR = "data/GeneCalls"
OUTDIR = "data/ParsedAnvioAnnot"

def import_func_data(gene_fpath, assembly_name):
    """
    Returns df of functional annotation for genes. 
    """
    fpath = f'{gene_fpath.parent}/{assembly_name}_functions.tsv'

    df = pd.read_table(fpath)

    # obtain top annotation (smallest e-value)
    df = df.sort_values(by=['gene_callers_id', 'e_value'], ascending=[True, True])  
    df = df.groupby(['gene_callers_id', 'source']).head(1)  

    # long to wide pivot - functions 
    function_df = df.pivot(
        index=['gene_callers_id'], 
        columns='source', values='function'
    ).reset_index()

    # long to wide pivot - accession 
    df['accession_type'] = df['source'] + '_ACC'
    accession_df = df.pivot(
        index=['gene_callers_id'], 
        columns='accession_type', values='accession'
    ).reset_index()

    df = pd.merge(function_df, accession_df, on=['gene_callers_id'])

    return df

def add_intergenic(df):
    """
    Create rows for intergenic regions
    """
    columns = df.columns
    # list to store new rows for intergenic regions
    intergenic_rows = []  

    # group by contigs
    contig_groups = df.groupby('contig')

    for index, cdf in contig_groups:
        for i in range(len(cdf) - 1):
            # obtain rows 
            current_row = cdf.iloc[i]
            next_row = cdf.iloc[i + 1]

            # check if there is a gap (intergenic) between rows/genes
            if next_row['start'] > current_row['stop'] + 1:
                # create dict for new row with all columns - empty values 
                new_row = {
                    col_name: '' for col_name in columns
                }

                # obtain intergenic region start and end 
                start = current_row['stop'] + 1
                stop = next_row['start'] - 1

                # add in values 
                new_row['gene_callers_id'] = f'{current_row['gene_callers_id']}.1'
                new_row['contig'] = current_row['contig']
                new_row['start'] = start  # Start right after the 'end' of the current row
                new_row['stop'] = stop  # End right before the 'start' of the next row

                # Append the new row to the list
                intergenic_rows.append(new_row)

    intergenic_df = pd.DataFrame(intergenic_rows)

    # add intergenic to original df 
    intergenic_df['type'] = 'intergenic'
    df['type'] = 'gene'
    df = pd.concat([df, intergenic_df]).sort_values(by='start').reset_index(drop=True)

    return df


def process_sample(gene_fpath, OUTDIR):
    """
    """
    # obtain assembly data from path to gene calls 
    assembly_name = gene_fpath.stem.replace('_genecalls', '')

    # import genes and their coordinates
    gene_df = pd.read_table(gene_fpath, usecols=['gene_callers_id', 'contig', 'start', 'stop'])

    # import genes and their functions/accesssions
    func_df = import_func_data(gene_fpath, assembly_name)

    df = pd.merge(gene_df, func_df, on=['gene_callers_id'], how='inner')

    # add annot rows for intergenic regions 
    df = add_intergenic(df)
    
    df.to_csv(f'{OUTDIR}/{assembly_name}.tsv', sep='\t', index=False)
    

def main():
    Path(OUTDIR).mkdir(exist_ok=True, parents=True)

    for gene_fpath in Path(GENE_CALL_DIR).glob('*_genecalls.tsv'):
        process_sample(gene_fpath, OUTDIR)


main()