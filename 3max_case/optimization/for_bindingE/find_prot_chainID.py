###########################################################################
# This script will find the chain ID of the complex PDB
#
# Written by Xingcheng Lin, 01/02/2023;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np
import mdtraj as md
import itertools

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################

def find_prot_chainID(pdbfile, chain_ID_protein_file):

    # Load the PDB
    pdb = md.load(pdbfile)
    # Convert into a pandas table
    pdb_table, pdb_bonds = pdb.topology.to_dataframe()

    # Select rows belonging to protein
    prot_protein_table = pdb_table.loc[(pdb_table['resName'] == 'ALA') |  (pdb_table['resName'] == 'ARG') | (pdb_table['resName'] == 'ASN') \
         | (pdb_table['resName'] == 'ASP') | (pdb_table['resName'] == 'CYS') | (pdb_table['resName'] == 'GLU') \
            | (pdb_table['resName'] == 'GLN') | (pdb_table['resName'] == 'GLY') | (pdb_table['resName'] == 'HIS') \
                | (pdb_table['resName'] == 'ILE') | (pdb_table['resName'] == 'LEU') | (pdb_table['resName'] == 'LYS') \
                    | (pdb_table['resName'] == 'MET') | (pdb_table['resName'] == 'PHE') | (pdb_table['resName'] == 'PRO') \
                        | (pdb_table['resName'] == 'SER') | (pdb_table['resName'] == 'THR') | (pdb_table['resName'] == 'TRP') \
                            | (pdb_table['resName'] == 'TYR') | (pdb_table['resName'] == 'VAL')]

    # We just need the first protein atom to tell us the protein chain ID, suppose there is only one protein chain;
    prot_chainID = prot_protein_table.iloc[0]["chainID"]

    chainID_map = {
    'A': 0,
    'B': 1,
    'C': 2,
    'D': 3,
    'E': 4,
    'F': 5,
    'G': 6,
    'H': 7,
    'I': 8,
    'J': 9,
    'K': 10,
    'L': 11,
    'M': 12,
    'N': 13,
    'O': 14,
    'P': 15,
    'Q': 16,
    'R': 17,
    'S': 18,
    'T': 19,
    'U': 20,
    'V': 21,
    'W': 22,
    'X': 23,
    'Y': 24,
    'Z': 25
    }

    value_to_find = prot_chainID
    for key, value in chainID_map.items():
        if value == int(value_to_find):
            prot_chainID_in_letter = key

    # Write the chain ID into a file;
    outfile = open(chain_ID_protein_file, 'w')
    outfile.write(str(prot_chainID_in_letter) + "\n")
    outfile.close()

    return

############################################################################

if __name__ == "__main__":
    pdbfile = sys.argv[1]
    chain_ID_protein_file = sys.argv[2]

    find_prot_chainID(pdbfile, chain_ID_protein_file)
 
 #   print("I slept and dreamt that life was joy. I awoke and saw that life was service. I acted and behold, service was joy.")
