"""
Purpose: to parse methylated sites and plot box plot of: 
    - samples with motifs and number of "high quality" sites (>90 percent modified)
    - samples without motifs and number of "high quality" sites (>90 percent modified)

Reason: Do samples with motifs called have higher “high quality” methylated sites? 
    - There are samples without any motifs detected. Their motif calling results had very low
    e-values (MicrobeMod cutoff of 0.1)
"""
import numpy as np
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt

def process_df(fpath):
    """
    Returns df of sample with cols: 
        - Type: whether or not motif was detected in sample 
        - High Quality Count: counts of high quality (>90% methylated) sites
    """
    df = pd.read_table(fpath)[['Percent_modified', 'motif']]

    # obtain number of high-quality sites 
    high_qual_count = len(df[df['Percent_modified'] > 0.9])
    
    # obtain "type" (whether or not there were motifs called)
    # no motifs if all rows are NA
    type = "No Motif Detected" if df['motif'].isna().sum() == len(df) else "Motif Detected"

    # create df for data 
    df = pd.DataFrame({
        'Type': [type], 
        'High Quality Count': [high_qual_count], 
    })

    return df

def main():
    # make output dir 
    outdir = "data/motif_v_modification"
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # metadata df 
    sdf = pd.read_table("../../MicrobeMod-Snakemake/inputs/samples.tsv")[['Sample name', 'sample', 'genus']]

    # cycle through each file and parse
    call_methylation_dir = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/microbemod_call_methylation"

    dfs = []
    for fpath in Path(call_methylation_dir).glob('*_methylated_sites.tsv'):
        fname = fpath.stem.replace('_methylated_sites', '')  # file name

        # process df
        df = process_df(fpath)

        dfs.append(df)

    df = pd.concat(dfs)
    df = df.sort_values(by=['Type'])

    # save df
    df.to_excel(f'{outdir}/motif_v_modification.xlsx', index=False)

    # log10 data (no motif detected has really low counts - skewed results)
    df['log10 High Quality Count'] = df['High Quality Count'].apply(np.log10)
    df = df[df['High Quality Count'] < 10000]

    # obtain lists for plotting
    motif = df[df['Type'] == 'Motif Detected']['log10 High Quality Count'].values.tolist()
    no_motif = df[df['Type'] == 'No Motif Detected']['log10 High Quality Count'].values.tolist()

    fig, ax = plt.subplots()

    ax.boxplot(
        [motif, no_motif], 
        labels=['Motif Detected', 'No Motif Detected']
    )

    ax.set_ylabel('log10 High Quality Count')
    fig.suptitle('Motif Detection vs. Counts of High Quality Sites')

    plt.savefig(f'{outdir}/motif_v_modification.png')
    plt.close()




main()