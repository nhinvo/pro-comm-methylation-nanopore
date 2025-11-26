"""
Purpose: to plot distribution of percent methylated value. 

Reason: to determine a good cutoff for STREME in MicrobeMod. 
"""
import pandas as pd 
from pathlib import Path 
import seaborn as sns
import matplotlib.pyplot as plt

def plot_dist(data_dict, fname, sdf, outdir):
    """
    Plot distribution of percent methylation for input sample. 
    """
    df = data_dict['df']
    
    # plot distribution 
    dist_fig = sns.displot(
        df['Percent_modified'], #  kind='kde', 
    )

    # plot mean and median
    plt.axvline(x=data_dict['mean'], color='red')
    plt.axvline(x=data_dict['median'], color='blue')

    # obtain sample data
    org = sdf.loc[sdf['sample'] == fname, 'genus'].values[0]
    sample = sdf.loc[sdf['sample'] == fname, 'Sample name'].values[0]

    dist_fig.fig.suptitle(f"Percent Methylation Distribution\nSample: {sample}\nOrganism: {org}", fontsize=13)
    dist_fig.fig.subplots_adjust(top=0.90)  

    plt.savefig(f'{outdir}/plots/{fname}.png')
    plt.close()

def process_df(fpath):
    """
    """
    df = pd.read_table(fpath)[['Total_coverage', 'Percent_modified', 'motif']]

    data_dict = {
        'df': df[['Percent_modified']], 
        'mean': df['Percent_modified'].mean(), 
        'median': df['Percent_modified'].median()
    }

    return data_dict

def main():
    # make output dir 
    outdir = "data/modification_distribution"
    Path(outdir).mkdir(parents=True, exist_ok=True)
    Path(f"{outdir}/plots").mkdir(parents=True, exist_ok=True)
    
    # metadata df 
    sdf = pd.read_table("../../MicrobeMod-Snakemake/inputs/samples.tsv")[['Sample name', 'sample', 'genus']]

    # cycle through each file, process, and plot 
    call_methylation_dir = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/microbemod_call_methylation"

    dfs = []
    for fpath in Path(call_methylation_dir).glob('*_methylated_sites.tsv'):
        fname = fpath.stem.replace('_methylated_sites', '')  # file name

        # process df
        data_dict = process_df(fpath)

        # plot distribution
        plot_dist(data_dict, fname, sdf, outdir)

        # make df of median and mean for sample 
        df = pd.DataFrame({
            'sample': [fname], 
            'mean': data_dict['mean'], 
            'median': data_dict['median'], 
        })

        dfs.append(df)

    df = pd.concat(dfs)
    df = df.sort_values(by=['sample'])
    df.to_excel(f'{outdir}/percent_modification_summary.xlsx', index=False)
    


main()