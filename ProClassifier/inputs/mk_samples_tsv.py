"""
"""
import pandas as pd 
from pathlib import Path 

def main():
    trimmed_dir = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/MAG_Snakemake/filtered_reads"
    df_data = []
    for fpath in Path(trimmed_dir).glob('*fastq'):
        sname = fpath.stem.replace('_filtered', '')
        data = {
            'sample': sname, 
            'read_path': fpath, 
        }

        df_data.append(data)

    df = pd.DataFrame(df_data)
    df.to_csv('samples.tsv', sep='\t', index=False)


main()