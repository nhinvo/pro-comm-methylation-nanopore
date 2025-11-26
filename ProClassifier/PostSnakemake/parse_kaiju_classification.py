import pandas as pd 
from pathlib import Path 

def main():
    dir = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/ProClassifier/prosyn_reads/read_name_classification"
    for fpath in Path(dir).glob('*'):
        if 'barcode22' not in str(fpath):
            continue

        df = pd.read_table(fpath, names=['name', 'classification'])[['classification']]
        df['genus'] = df['classification'].str.split('; ').str[8]
        df['clade'] = df['classification'].str.split('; ').str[10]
        # print(df['classification'].iloc[0])

        print(fpath)
        print(df['clade'].value_counts())
        print()
        print()
        

        # break

main()