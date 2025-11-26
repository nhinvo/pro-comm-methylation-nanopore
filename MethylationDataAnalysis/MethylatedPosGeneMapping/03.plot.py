"""
Purpose: to plot data from 02.analysis.py

Nhi N Vo - 06/30/25
"""
from pathlib import Path 
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
from COG_labels import *

ASSEMBLY_DF = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/PostSnakemake/ExtractAssemblies/HighQualAssemblies.tsv"
SUMMARY_TSV = "data/MethylationSummary.tsv"
GENE_POS_TSV = "data/MethylatedGenePos.tsv"
COG_DIR = 'data/Methylated_COG'

TREATMENT_ORDER = [
    'MED4ax', 
    'MED4ax+A1A', 
    'MED4ax+1907',
    'MED4ax+1940', 
    'MED4ax+1943', 
    'MED4Xe', 
    'NATL2Aax',
    'NATL2Aax+A1A',
]

def plot_summary_df(name, df, outdir, genus):
    """
    """
    # obtain mean and std across reps 
    col = 'Percent_gene_methylated' if name == 'gene' else 'Percent_intergenic_methylated'
    summary_df = df.groupby('Sample name').agg({
        col: ['mean', 'std'], 
    })
    summary_df = summary_df.fillna(0)

    # obtain bar_df for plotting bars
    mean_df = summary_df.loc[:, [(col, 'mean')]]
    std_df = summary_df.loc[:, [(col, 'std')]]

    # plotting 
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))
    # bar_df.plot.bar(ax=ax, width=0.7, edgecolor='black', linewidth=1)

    bar_color = 'tab:blue' if name == 'gene' else 'tab:orange'
    ax.bar(
        mean_df.index.values,  # sample name 
        mean_df.values.flatten(),  # mean values 
        yerr=std_df.values.flatten(),  # std 
        capsize=4, edgecolor='black', 
        color = bar_color, 
    )

    for i in range(len(summary_df)):
        bar_sname = summary_df.iloc[i].name

        scatter_y = df[df['Sample name'] == bar_sname][col].values
        scatter_x = [bar_sname] * len(scatter_y)

        sns.stripplot(x=scatter_x, y=scatter_y, color='red', jitter=True, edgecolor='black', size=7, linewidth=1)

    # plot edits 
    ax.set_xlabel('', fontsize=1)  
    ax.set_ylabel('Percent (%)', fontsize=1)  # remove y-axis title ("genus")
    ax.set_ylim(0, 100)  # edit x-axis limits 
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_title(f'{genus}: Percent {name.title()} Methylated', fontsize=10, loc='center')
    
    plt.tight_layout()
    plt.savefig(f'{outdir}/{genus}_{name}.png')

    plt.close()

def plot_summary(SUMMARY_TSV):
    """
    """
    outdir = 'data/plots/MethylationSummary'
    Path(outdir).mkdir(exist_ok=True, parents=True)

    df = pd.read_table(SUMMARY_TSV)

    # treatment name edits 
    df['Sample name'] = df['Sample name'].str.strip().str.replace(' ', '')

    # group by genus 
    genus_groups = df.groupby(['genus'])

    for index, df in genus_groups:
        genus = index[0]
        gene_df = df[['genus', 'Sample name', 'rep', 'Percent_gene_methylated']]
        intergenic_df = df[['genus', 'Sample name', 'rep', 'Percent_intergenic_methylated']]

        for name, df in {'gene':gene_df, 'intergenic': intergenic_df}.items():
            plot_summary_df(name, df, outdir, genus)


def plot_summary_v1(SUMMARY_TSV):
    """
    Plot bar of percent genes/intergenic methylated. 
    """
    outdir = 'data/plots/MethylationSummary_v1'
    Path(outdir).mkdir(exist_ok=True, parents=True)

    df = pd.read_table(SUMMARY_TSV)
    groups = df.groupby(['Sample name'])
    
    for index, df in groups:
        treatment_name = index[0]
        df = df[['genus', 'rep', 'Percent_gene_methylated', 'Percent_intergenic_methylated']]
        df = df.fillna(0)

        # obtain mean and std 
        summary_df = df.groupby('genus').agg({
            'Percent_gene_methylated': ['mean', 'std'], 
            'Percent_intergenic_methylated': ['mean', 'std'], 
        })

        # obtain bar_df for plotting bars
        bar_df = summary_df.loc[:, [
            ('Percent_gene_methylated', 'mean'), 
            ('Percent_intergenic_methylated', 'mean')
        ]]

        # plotting 
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        bar_df.plot.barh(rot=0, ax=ax, width=0.8, edgecolor='black', linewidth=0.7)


        # plot edits 
        ax.set_xlabel('Percent (%)', fontsize=10)  
        ax.set_ylabel('', fontsize=1)  # remove y-axis title ("genus")
        ax.set_xlim(0, 100)  # edit x-axis limits 


        # legend 
        legend = ax.legend(
            ['Percent Gene Methylated', 'Percent Intergenic Methylated'], 
            loc='center', bbox_to_anchor=(0.5, 1.10),  # position of legend
            title=f'{treatment_name} Methylation', title_fontsize=13, 
            alignment='center',  # center align texts 
            fontsize=10, 
            frameon=False,  # remove square frame
            handleheight=1, handlelength=1,  # make markers square instead of rectangles
            handletextpad=0.5, labelspacing=0.25,  # no spacing between markers vertically 
            ncol=2, 
        )  
        
        plt.tight_layout()
        plt.savefig(f'{outdir}/{treatment_name}.png')


def plot_methylated_pos_all(df, outdir):
    """
    """
    # get data 
    x = df['Position_scaled']
    y = df['percentage']

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    sns.boxplot(x=x, y=y)
    # sns.stripplot(x=x, y=y, color='blue', jitter=True, edgecolor='black', size=7, linewidth=1, ax=ax)

    # plot edits 
    ax.set_xlabel('Gene Position', fontsize=10)  
    ax.set_ylabel('Percent (%) of methylated sites within sample', fontsize=10)
    # ax.tick_params(axis='x', labelrotation=90)
    ax.set_title(f'Percent of Methylated sites at each gene position', fontsize=10, loc='center')

    plt.tight_layout()
    plt.savefig(f'{outdir}/AllSamples.png')
    plt.close()


def plot_methylated_pos_genus(df, outdir):
    """
    """
    # sample name edit 
    df['Sample name'] = df['Sample name'].str.strip().str.replace(' ', '')

    # group by genus 
    genus_groups = df.groupby(['genus'])

    for index, gdf in genus_groups:
        genus = index[0]

        # plot each subplot for an exp 
        treatment_groups = gdf.groupby(['Sample name'])
        treatment_groups = [(index[0], tdf) for index, tdf in treatment_groups]

        # determine plot layout
        num_subplots = len(treatment_groups)
        num_cols = int(num_subplots ** 0.5 + 1)
        num_rows = (num_subplots + num_cols - 1) // num_cols
        figsize = (num_cols * 3, num_rows * 3)
        fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)

        for count, ax in enumerate(fig.axes):
            try: 
                treatment_data = treatment_groups[count]
            except:
                ax.set_axis_off()
                continue

            # plot subplot 
            treatment = treatment_data[0]
            tdf = treatment_data[1]

            # scatter plotting 
            sns.stripplot(
                x=tdf['Position_scaled'], y=df['percentage'], 
                color='blue', 
                jitter=True, 
                edgecolor='black', 
                size=7, linewidth=1, 
                ax=ax
            )

            # edits 
            ax.set_title(f"{treatment}")
            ax.tick_params(axis='x', labelrotation=90)

        fig.suptitle(f'{genus}')

        # save plot for genus 
        plt.tight_layout()
        plt.savefig(f'{outdir}/{genus}.png')
        plt.close()


def plot_methylated_pos(GENE_POS_TSV):
    """
    """
    outdir = "data/plots/GeneMethylationDistribution"
    Path(outdir).mkdir(parents=True, exist_ok=True)

    df = pd.read_table(GENE_POS_TSV)

    # plot all pos and percentages (across all samples)
    plot_methylated_pos_all(df, outdir)

    # plot each genus separately 
    plot_methylated_pos_genus(df, outdir)


def plot_COG_old(COG_DIR, ASSEMBLY_DF):
    """
    """
    # input dir 
    COG_DIR = Path(COG_DIR)
    
    # output dir 
    outdir = Path('data/plots/Methylated_COG')
    outdir.mkdir(parents=True, exist_ok=True)

    # samples tsv metadata 
    sdf = pd.read_table(ASSEMBLY_DF)[['Sample name', 'rep', 'assembly_name', 'genus']]

    for fpath in COG_DIR.glob('*'):
        # obtain metadata 
        assembly_name = fpath.stem
        treatment = sdf[sdf['assembly_name'] == assembly_name]['Sample name'].item()
        rep = sdf[sdf['assembly_name'] == assembly_name]['rep'].item()
        genus = sdf[sdf['assembly_name'] == assembly_name]['genus'].item()

        df = pd.read_table(fpath)
        df = df[~df['COG20_CATEGORY_ACC'].isin(['A', 'B'])]  # remove A 

        # df = df[~df['COG20_CATEGORY_ACC'].str.contains('No Annotation')]

        # obtain color list based on order of row 
        colors = [COG_COLORS[cog] for cog in df['COG20_CATEGORY_ACC'].values.tolist()]

        # convert cog labels (e.g "C") into full descriptions 
        df['COG20_CATEGORY_ACC'] = df['COG20_CATEGORY_ACC'].replace(COG_LABELS)

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        df.plot.barh(x='COG20_CATEGORY_ACC', y='count', ax=ax, color=colors)

        fig.suptitle(f'COG of Methylated genes - {treatment} ({rep}) - {genus}')
        ax.set_xlabel('Gene Count', fontsize=10) 
        ax.set_ylabel('', fontsize=1)  # remove y axis label 

        ax.legend().set_visible(False)

        fig.tight_layout()  # Adjust layout to prevent clipping of labels
        fig.savefig(f'{outdir}/{assembly_name}.png')
        plt.close()


def plot_COG(COG_DIR, ASSEMBLY_DF):
    """
    """
    # input dir 
    COG_DIR = Path(COG_DIR)
    
    # output dir 
    outdir = Path('data/plots/Methylated_COG')
    outdir.mkdir(parents=True, exist_ok=True)

    # samples tsv metadata 
    sdf = pd.read_table(ASSEMBLY_DF)[['Sample name', 'rep', 'assembly_name', 'genus']]

    for fpath in COG_DIR.glob('*'):
        # obtain metadata 
        assembly_name = fpath.stem
        treatment = sdf[sdf['assembly_name'] == assembly_name]['Sample name'].item()
        rep = sdf[sdf['assembly_name'] == assembly_name]['rep'].item()
        genus = sdf[sdf['assembly_name'] == assembly_name]['genus'].item()

        df = pd.read_table(fpath)
        df = df.sort_values(by=['methylated_percent'])
        df = df[~df['COG20_CATEGORY_ACC'].isin(['A', 'B'])]  # remove A 

        # df = df[~df['COG20_CATEGORY_ACC'].str.contains('No Annotation')]

        # obtain color list based on order of row 
        colors = [COG_COLORS[cog] for cog in df['COG20_CATEGORY_ACC'].values.tolist()]

        # convert cog labels (e.g "C") into full descriptions 
        df['COG20_CATEGORY_ACC'] = df['COG20_CATEGORY_ACC'].replace(COG_LABELS)

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        df.plot.barh(x='COG20_CATEGORY_ACC', y='methylated_percent', ax=ax, color=colors)

        fig.suptitle(f'COG of Methylated genes - {treatment} ({rep}) - {genus}')
        ax.set_xlabel('Percent Gene Methylated', fontsize=10) 
        ax.set_ylabel('', fontsize=1)  # remove y axis label 

        ax.legend().set_visible(False)

        fig.tight_layout()  # Adjust layout to prevent clipping of labels
        fig.savefig(f'{outdir}/{assembly_name}.png')
        plt.close()

def main():
    # # plot percent gene/intergenic methylated 
    # plot_summary(SUMMARY_TSV)
    # # plot_summary_v1(SUMMARY_TSV)  # old version, not ideal but keep in case

    # # plot distribution of methylated position on gene 
    # plot_methylated_pos(GENE_POS_TSV)

    # plot COG count 
    plot_COG(COG_DIR, ASSEMBLY_DF)

main()