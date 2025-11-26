"""
Purpose: to parse all outputs from MicrobeMod annotate_rm and call_methylation. 
"""
import pandas as pd 
from pathlib import Path 

def parse_output(input, extension, output, samples_df):
    """
    """
    # convert to list of path if input is string (1 file)
    paths = [input] if isinstance(input, str) else input

    dfs = []  # list to store motifs.tsv from all samples 
    for fpath in paths:
        fpath = Path(fpath)
        fname = fpath.stem.replace(extension, '')
        print(f"Importing file {fname} at path: {fpath}\n")

        # import table 
        try:
            df = pd.read_table(fpath)
        except pd.errors.EmptyDataError:
            print(f'No results for {fname}, skipping import.\n')
            continue 
        
        df['sample'] = fname
        dfs.append(df)

    if dfs:
        df = pd.concat(dfs)
        df = pd.merge(samples_df, df, on='sample', how='outer')
    else:
        df = pd.DataFrame()

    output = Path(output)
    df.to_csv(output, sep='\t', index=False)
    excel_output = output.parent / f"{output.stem}.xlsx"
    df.to_excel(excel_output, index=False)

def main():
    # import sample table and metadata 
    samples_df = pd.read_table(snakemake.config['input']['sample table'])

    # parse call methylation output and save as tsv + excel 
    parse_output(
        snakemake.input['call_methylation_output'], 
        '_motifs', 
        snakemake.output['parsed_call_methylation'], 
        samples_df
    )

    # parse annotate rm output and save as tsv + excel 
    parse_output(
        snakemake.input['annotate_rm_output'], 
        '.rm.genes', 
        snakemake.output['parsed_annotate_rm'], 
        samples_df,
    )

main()