import os 
import pandas as pd
from pathlib import Path 

DEPTH_DIR = "../SamtoolsDepth/data/samtools_depth"
ANNOT_DIR = "../../../AnvioGeneAnnotation/data/ParsedAnvioAnnot"

def obtain_depth_zscore(DEPTH_DIR):
    """
    """
    output_fpath = 'data/depth_zscore.tsv'

    if not os.path.exists(output_fpath):
        dfs = []

        for fpath in Path(DEPTH_DIR).glob('*'):
            sname = fpath.stem 

            df = pd.read_table(fpath, names=['contig', 'pos', 'count'])
            
            mean = df['count'].mean()
            std = df['count'].std()

            df['count z-score'] = (df['count'] - mean) / std

            # obtain values > 7 std 
            df = df[df['count z-score'] > 7]
            df['assembly_name'] = sname 

            dfs.append(df)

        df = pd.concat(dfs)
        df.to_csv(output_fpath, sep='\t', index=False)

        return df

    else:
        df = pd.read_table(output_fpath)
        return df

def map_peaks(df, ANNOT_DIR):
    """
    """
    dfs = []
    for fpath in Path(ANNOT_DIR).glob('*'):
        assembly_name = fpath.stem 

        # obtain rows of assembly 
        peak_df = df[df['assembly_name'] == assembly_name]
        # peak_df = peak_df.head(500)  # TEMP 

        # skip samples without peaks 
        if len(peak_df) == 0:
            continue 

        # import annotations 
        annot_df = pd.read_table(fpath)

        # merge so all peaks map to all annotations 
        merged_df = pd.merge(peak_df, annot_df, how='cross')

        # filter to keep rows where mutations are within the gene/intergenic regions
        merged_df = merged_df[
            (merged_df['contig_x'] == merged_df['contig_y'])  # same contig 
            & (merged_df['pos'] >= merged_df['start']) 
            & (merged_df['pos'] <= merged_df['stop'])
        ]

        merged_df = merged_df.drop_duplicates(subset=['gene_callers_id'])

        # col edits
        merged_df = merged_df.drop(columns=['pos', 'count', 'count z-score', 'contig_y'])
        merged_df = merged_df.rename(columns={'contig_x': 'contig'})

        dfs.append(merged_df)

    df = pd.concat(dfs)
    df.to_csv('data/peak_annot.tsv', sep='\t', index=False)


def main():
    Path('data').mkdir(parents=True, exist_ok=True)

    # depth zscored - filtered for > 7 (peaks)
    df = obtain_depth_zscore(DEPTH_DIR)

    # map peaks to regions
    map_peaks(df, ANNOT_DIR)

        
main()