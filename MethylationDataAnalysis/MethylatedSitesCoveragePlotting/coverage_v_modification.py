"""
Purpose: to parse mythylated site data & plot [coverage vs. percent modified]. 
"""
import pandas as pd 
import numpy as np
from pathlib import Path 
import seaborn as sns 
import matplotlib.pyplot as plt 

def plot(df, assembly_name, sdf, cov_df):
    """
    Plot Coverage vs. Methylation Percent for all methylated positions. 
    """
    # custom colors
    palette = {
        'Motif Detected': '#70cc73', 
        'No Motif Detected': '#d46161', 
    }

    scatter_fig = sns.jointplot(
        data=df, 
        x="Total_coverage", 
        y="Percent_modified", 
        hue="motif", 
        s=10, alpha=0.5,  # linewidth=0.2, edgecolor="black", 
        palette=palette, 
    )

    # plot average coverage line
    cov = cov_df.loc[cov_df['assembly_name'] == assembly_name, 'average_genome_cov'].item()
    plt.axvline(x=cov, color='grey', linestyle='--')

    # obtain sample data
    org = sdf.loc[sdf['assembly_name'] == assembly_name, 'genus'].item()
    sample = sdf.loc[sdf['assembly_name'] == assembly_name, 'Sample name'].item()
    rep = sdf.loc[sdf['assembly_name'] == assembly_name, 'rep'].item()

    scatter_fig.fig.suptitle(f"Coverage vs. Percent Methylation\nSample: {sample} ({rep})\nOrganism: {org}", fontsize=13)
    scatter_fig.fig.subplots_adjust(top=0.90)  

    plt.savefig(f'data/coverage_v_modification_plots/{assembly_name}.png')
    plt.close()

def process_df(fpath):
    """
    """
    df = pd.read_table(fpath)[['Total_coverage', 'Percent_modified', 'motif']]
    df['motif'] = df['motif'].apply(lambda x: "Motif Detected" if pd.notna(x) and x != '' else "No Motif Detected")

    return df

def main():
    Path('data/coverage_v_modification_plots').mkdir(exist_ok=True, parents=True)

    # metadata df 
    sdf = pd.read_table("../../GenomeAssembly-Snakemake/PostSnakemake/ExtractAssemblies/HighQualAssemblies.tsv")[['Sample name', 'assembly_name', 'genus', 'rep']]

    # coverage df 
    cov_df = pd.read_table("../../QC/GenomeCov/data/average_cov.tsv")

    # cycle through each file, process, and plot 
    call_methylation_dir = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/microbemod_call_methylation"
    for fpath in Path(call_methylation_dir).glob('*_methylated_sites.tsv'):
        fname = fpath.stem.replace('_methylated_sites', '')  # file name (assembly name)

        # process df
        df = process_df(fpath)

        # plot df
        plot(df, fname, sdf, cov_df)


main()