"""
Purpose: to parse the results from STREME. 

Note: based on MicrobeMod xml parser script
    - https://github.com/cultivarium/MicrobeMod/blob/main/MicrobeMod/microbemod.py

04/17/25
"""
import pandas as pd 
from pathlib import Path 
import xml.etree.ElementTree as ET

CALL_METHYLATIOB_DIR = "/nobackup1b/users/chisholmlab/nvo/24_allison_methylation/Nanopore/MicrobeMod-Snakemake/trash/microbemod_call_methylation_streme_0.8"
MIN_EVALUE = 0.1

def parse_streme_xml(fpath):
    """
    Parse STREME .xml output file. 
    Note: script from MicrobeMod (link above).
    """
    tree = ET.parse(fpath)
    root = tree.getroot()

    motifs = []

    for motif in root[1]:
        evalue = float(motif.get("test_evalue"))

        if evalue < MIN_EVALUE:
            motif_old = motif.get("id")
            motif_new = ""

            print(motif_old)

            # convert to N's 
            i = 0
            for pos in motif:
                freqs = [float(x) for x in pos.attrib.values()]
                print(pos)
                print(freqs)

def main():
    for fpath in Path(CALL_METHYLATIOB_DIR).glob('*/*.xml'):
        parse_streme_xml(fpath)

        break


main()