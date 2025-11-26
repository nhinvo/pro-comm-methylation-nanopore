"""
Parse samtools depth output to obtain average coverage. 
"""
import os
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt 

DEPTH_DIR = "data/samtools_depth"
SAMPLES_DF = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/PostSnakemake/ExtractAssemblies/HighQualAssemblies.tsv"

def parse_depth(depth_dir, sdf):
    """
    """
    out_fpath = "data/average_cov.tsv"

    # check if the output file already exists
    if os.path.exists(out_fpath):
        print(f"The file '{out_fpath}' already exists. Importing already-made table...")

        df = pd.read_table(out_fpath)

        return df

    # else run code 
    else:
        data = {
            'assembly_name': [], 
            'average_genome_cov': [], 
        }

        for fpath in Path(depth_dir).glob('*tsv'):
            fname = fpath.stem
            df = pd.read_table(fpath, names=['contig', 'pos', 'cov'])

            data['assembly_name'].append(fname)
            data['average_genome_cov'].append(df['cov'].mean())

        df = pd.DataFrame(data)

        df = pd.merge(df, sdf, on=['assembly_name'], how='inner')

        df = df.sort_values(by=['bam_filename'], ascending=[True])

        df.to_csv('data/average_cov.tsv', sep='\t', index=False)

        # excel for dropbox 
        cols = ['assembly_name', 'average_genome_cov', 'Sample name', 'rep', 'genus', 'Completeness', 'Contamination', 'Contig_N50', 'bin_contig_num']
        df[cols].to_excel('data/average_cov.xlsx', index=False)

        return df 

def edit_pos(df):
    """
    Some samples have multiple contigs, edit position to be continuous. 
    """
    contig_groups = df.groupby(['contig'])

    i = 0
    dfs = []
    for index, cdf in contig_groups:
        if i == 0:
            # first contig, obtain the last pos 
            contig_last_pos = int(cdf.tail(1)['pos'].item())

            # do nothing to the pos in this contig 
            dfs.append(cdf)

            i += 1

            continue 

        else:
            # add the previous contig last position to pos in current contig 
            cdf['pos'] = cdf['pos'] + contig_last_pos

            # obtain the last position after addition 
            contig_last_pos = int(cdf.tail(1)['pos'].item())

            dfs.append(cdf)

    df = pd.concat(dfs)
    return df


def plot_genome_cov(DEPTH_DIR, sdf):
    """
    For each genome, plot the coverage
    """
    plot_outdir = 'data/plots/genome_cov'
    Path(plot_outdir).mkdir(parents=True, exist_ok=True)

    for fpath in Path(DEPTH_DIR).glob('*tsv'):
        # obtain assembly metadata 
        assembly_name = fpath.stem
        treatment = sdf.loc[sdf['assembly_name'] == assembly_name, 'Sample name'].item()
        rep = sdf.loc[sdf['assembly_name'] == assembly_name, 'rep'].item()
        genus = sdf.loc[sdf['assembly_name'] == assembly_name, 'genus'].item()
        N50 = sdf.loc[sdf['assembly_name'] == assembly_name, 'Contig_N50'].item()
        contig_num = sdf.loc[sdf['assembly_name'] == assembly_name, 'bin_contig_num'].item()

        # obtain plot data 
        df = pd.read_table(fpath, names=['contig', 'pos', 'cov'])
        contig_num = df['contig'].nunique()

        if contig_num > 1:
            df = edit_pos(df)
        
        x = df['pos']
        y = df['cov']

        fig, ax = plt.subplots(1, 1, figsize=(7, 2))

        # plot cov 
        ax.plot(x, y, linewidth=0.5)
        
        # plot average cov 
        avg_cov = (sum(y)/len(y))
        ax.axhline(y=avg_cov, color='gray', linestyle='--',)

        ax.set_title(
            f'{treatment} - {rep} - {genus} - mean cov: {round(avg_cov,2)}, N50: {N50}, contig_num: {contig_num}', 
            fontsize=10
        )

        fig.savefig(f'{plot_outdir}/{assembly_name}.png')
        plt.close()


def main():
    sdf = pd.read_table(SAMPLES_DF)

    # parse depth and obtain average
    avg_df = parse_depth(DEPTH_DIR, sdf)

    # plot coverage across genome 
    plot_genome_cov(DEPTH_DIR, sdf)


main()