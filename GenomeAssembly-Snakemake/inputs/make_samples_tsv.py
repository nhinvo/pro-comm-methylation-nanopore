"""
Purpose: example script to make samples.tsv for ONT assembly workflow 
    - DEMUX_FASTQ_DIR: directory with input fastq files 
    - INPUT_EXTENSION: file extension of input fastq file in directory 

    - Note: Nanopore read files are single end 

01/08/25
"""
import pandas as pd 
from pathlib import Path 

DEMUX_FASTQ_DIR = "/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240730Chi/PseudoID-24/20240805_1712_1B_PAW59116_2cfaf0b1/demux_fastq"
INPUT_EXTENSION = ".fastq.gz"

def main():    
    # dict to store file paths 
    data = {
        'sample': [], 
        'read_path': [], 
    }

    # loop through all files that match extension in input directory 
    for fpath in Path(DEMUX_FASTQ_DIR).glob(f'*{INPUT_EXTENSION}'):
        data['sample'].append(fpath.name.replace(INPUT_EXTENSION, ''))
        data['read_path'].append(fpath)

    df = pd.DataFrame(data)
    df.to_csv('samples.tsv', index=False, sep='\t')


main()