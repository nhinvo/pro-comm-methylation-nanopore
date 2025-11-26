"""
To extract high quality assemblies from the pipeline and save metadata. 

Outputs: 
    - HighQualAssemblies.tsv: high quality assemblies. 
    - ../../results/HighQualAssemblies directory: contains high quality assemblies. 

Inputs: 
    - BIN_DIR: path to directory with metabat bins 
    - ASSEMBLY_QUALITY_DIR: path to directory with checkm2 quality results 
    - CLASSIFICATION_DIR: path to GTDB-tk classification 
"""
import shutil
import pandas as pd 
import numpy as np
from pathlib import Path 

# INPUT PATHS #
SAMPLE_LOG_SHEET = "/orcd/data/chisholm/002/hstor001/nvo/projects/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/inputs/NanoporeDNAisolations.xlsx"
BAM_DIR = "/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/demux_modified_basecall"
BIN_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/metabat_binning/bins"
ASSEMBLY_QUALITY_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/checkm_bin_quality"
CLASSIFICATION_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/GenomeAssembly-Snakemake/gtdb_classification"

# OUTPUT PATHS #
ASSEMBLY_OUTDIR = "../../results/HighQualAssemblies"


def obtain_bam_paths(dir_path):
    """
    Return all Nanopore basecaling .bam paths 
    """
    bam_data = {
        'bam_filename': [], 
        'BMC Sample ID': [], 
        'bam_fpath': [], 
    }

    for fpath in Path(dir_path).glob('SQK*.bam'):
        sname = fpath.stem
        barcode = int(sname.strip().split('barcode')[1])  # e.g. 7 from SQK-NBD114-24_barcode07
        
        bam_data['bam_filename'].append(sname)
        bam_data['BMC Sample ID'].append(barcode)
        bam_data['bam_fpath'].append(fpath)

    df = pd.DataFrame(bam_data)

    return df

def process_bins(bin_dir):
    """
    Returns a df of information on all bins from metabat2. 
    """
    data = []

    # cycle through every bins 
    for bin_fpath in Path(bin_dir).glob('*/*.fa'):
        sample = bin_fpath.parent.name
        bin_id = bin_fpath.stem.split('.')[-1]  # e.g. '1' from full/file/path/Paragon14_34.1.fa
            
        # count number of contigs and length 
        contig_count = 0
        bin_len = 0
        with open(bin_fpath, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    contig_count += 1
                else:
                    line_len = len(line)
                    bin_len += line_len

        bin_data = {
            'bam_filename': sample, 
            'bin_id': bin_id, 
            'bin_contig_num': contig_count, 
            'bin_length': bin_len, 
            'bin_fpath': bin_fpath, 
        }

        data.append(bin_data)

    df = pd.DataFrame(data)

    return df 

def process_quality(quality_dir):
    """
    Returns df of aggregated data from checkM2 quality assessment. 
    """
    df = pd.concat(pd.read_table(fpath) for fpath in Path(quality_dir).glob('*/quality_report.tsv'))

    df['bam_filename'] = df['Name'].str.split('.').str[0]
    df['bin_id'] = df['Name'].str.split('.').str[1]

    return df

def process_classification(classification_dir):
    """
    Returns df of aggregate data from GTDB-tk classification. 
    """
    df = pd.concat(pd.read_table(fpath) for fpath in Path(classification_dir).glob('*/gtdbtk.bac120.summary.tsv'))

    df['bam_filename'] = df['user_genome'].str.split('.').str[0]
    df['bin_id'] = df['user_genome'].str.split('.').str[-1]

    # split up classification 
    taxonomic_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for counter, rank in enumerate(taxonomic_ranks):
        df['classification'] = df['classification'].str.replace('Unclassified Bacteria', 'd__Unclassified Bacteria')
        df[rank] = df['classification'].str.split(pat=';', expand=True)[counter].str[3:]

    return df 

def process_dfs(sample_df, bam_df, bin_df, quality_df, classification_df):
    """
    """
    # merge data 
    df1 = pd.merge(sample_df, bam_df, on='BMC Sample ID', how='inner')

    cols = ['bam_filename', 'bin_id']
    df2 = bin_df.merge(quality_df, on=cols, how='outer')\
        .merge(classification_df, on=cols, how='outer')

    df = pd.merge(df1, df2, on='bam_filename', how='inner')

    df = df.sort_values(by=['bam_filename', 'Completeness'], ascending=[True, False])

    # genome filtering
    df['Completeness'] = df['Completeness'].astype(float)
    df['Contamination'] = df['Contamination'].astype(float)
    df = df[(df['Completeness'] > 80) & (df['Contamination'] < 10)]

    df['assembly_name'] = df['bam_filename'] + '_' + df['bin_id']

    # remove cols 
    df = df.drop(columns=[
        'total vol (ul)',
        'total ng/ul',
        'total ug',
        'Completeness_Model_Used', 
        'Translation_Table_Used', 
        'Coding_Density', 
        'Additional_Notes', 
        'user_genome', 
        'classification', 
        'closest_genome_reference',
        'closest_genome_reference_radius',
        'closest_genome_taxonomy',
        'closest_genome_ani',
        'closest_genome_af',
        'closest_placement_reference',
        'closest_placement_radius',
        'closest_placement_taxonomy',
        'closest_placement_ani',
        'closest_placement_af',
        'pplacer_taxonomy',
        'classification_method',
        'note',
        'other_related_references(genome_id,species_name,radius,ANI,AF)',
        'msa_percent',
        'translation_table',
        'red_value',
        'warnings'

    ])

    df.to_csv('HighQualAssemblies.tsv', sep='\t', index=False)

    return df

def copy_assemblies(df, outdir):
    """
    Copy high quality assemblies to results/ directory 
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)

    bin_fpaths = df['bin_fpath'].values.tolist()
    for bin_fpath in bin_fpaths:
        assembly_fname = Path(bin_fpath).name
        source = bin_fpath 
        destination = f'{outdir}/{assembly_fname}'
        
        shutil.copy(source, destination)

    return True 


def main():
    # import data from assembly pipeline: sample metadata, bins, quality, and classification 
    sample_df = pd.read_excel(SAMPLE_LOG_SHEET, skiprows=1, sheet_name='BMC ready', usecols=range(6))
    bam_df = obtain_bam_paths(BAM_DIR)
    bin_df = process_bins(BIN_DIR)
    quality_df = process_quality(ASSEMBLY_QUALITY_DIR)
    classification_df = process_classification(CLASSIFICATION_DIR)

    # merge and process all data 
    df = process_dfs(
        sample_df, 
        bam_df, 
        bin_df, 
        quality_df, 
        classification_df
    )

    return False

    # copy extracted assemblies to results/ folder 
    copy_assemblies(df, ASSEMBLY_OUTDIR)

    

main()