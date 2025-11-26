import pandas as pd 
from pathlib import Path 
import xml.etree.ElementTree as ET

STREME_OUTDIR = "data/STREME_Outputs"
MIN_EVALUE = 0.1

def main():
    data = {
        'genome_name': [],
        'type': [], 
        'motif': [], 
        'e-value': [], 
    }

    for fpath in Path(STREME_OUTDIR).glob('*/streme.xml'):
        file_name = fpath.parent.name
        genome_name = '_'.join(file_name.split('_')[:-1])
        type = file_name.split('_')[-1]

        tree = ET.parse(fpath)
        root = tree.getroot()

        for motif in root[1]:
            evalue = float(motif.get("test_evalue"))

            if evalue < MIN_EVALUE:
                ### Evalue cutoff
                motif_old = motif.get("id").split("-")[1]
                
                data['genome_name'].append(genome_name)
                data['type'].append(type)
                data['motif'].append(motif_old)
                data['e-value'].append(evalue)

    df = pd.DataFrame(data)
    df = df.sort_values(['genome_name', 'type'], ascending=[False, False])
    df.to_csv('data/STREME_output.tsv', sep='\t', index=False)

main()