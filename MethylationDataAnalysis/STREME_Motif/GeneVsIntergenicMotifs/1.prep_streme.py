"""
Purpose: to create input files for STREME: 
    1. *_pos.fasta: 
    2. *_control.fasta: control (negative sequences) 
"""
import shutil 
import pandas as pd 
from pathlib import Path 

MICROBEMOD_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/microbemod_call_methylation"
ANNOT_DIR = "../../MethylatedPosGeneMapping/data/MethylatedPosAnnot"

OUTDIR = 'data/STREME_Inputs'

def main():
    Path(OUTDIR).mkdir(exist_ok=True, parents=True)

    # data to make df for STREME 
    STREME_data = {
        'file_name': [], 
        'pos_fpath': [], 
        'control_fpath': [], 
    }

    # cycle through each *_methylated_sites.tsv file
    for methylated_fpath in Path(MICROBEMOD_DIR).glob('*_methylated_sites.tsv'):
        # obtain genome name 
        genome_name = methylated_fpath.stem.replace('_methylated_sites', '')

        # import methylated sites
        methylated_df = pd.read_table(methylated_fpath, usecols=['Contig', 'Position', 'Sequence'])
        methylated_df = methylated_df.rename(columns={'Contig': 'contig'})

        # map pos with gene/intergenic 
        annot_fpath = f'{ANNOT_DIR}/raw_{genome_name}.tsv'
        annot_df = pd.read_table(annot_fpath, usecols=['contig', 'Position', 'type'])

        # merge
        df = pd.merge(methylated_df, annot_df, on=['contig', 'Position'], how='inner')
        
        # split methylated sites into gene/intergenic 
        data = {
            'gene': df[df['type'] == 'gene'],
            'intergenic': df[df['type'] == 'intergenic'],
        }

        # write sequences into files 
        pos_fpaths = []
        file_names = []
        for type, df in data.items():
            i = 0

            out_fpath = f'{OUTDIR}/{genome_name}_{type}_pos.fasta'
            pos_fpaths.append(out_fpath)
            file_names.append(f'{genome_name}_{type}')

            with open(out_fpath, 'w') as f:
                for _, row in df.iterrows():
                    i += 1
                    f.write(f'>{str(i)}\n')
                    f.write(f'{row['Sequence']}\n')

        # copy control sequences over 
        control_fpath = f'{MICROBEMOD_DIR}/{genome_name}_a_control.fasta'
        destination = f'{OUTDIR}/{genome_name}_control.fasta'
        shutil.copy(control_fpath, destination)

        # add data to STREME df 
        STREME_data['sample'].extend(file_names)
        STREME_data['pos_fpath'].extend(pos_fpaths)
        STREME_data['control_fpath'].extend([destination, destination])    

    df = pd.DataFrame(STREME_data)
    df.to_csv('data/STREME_input.tsv', sep='\t', index=False)


main()