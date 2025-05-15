#!/usr/bin/env python3

"""
This script extracts the chain ID of the protein in a given PDB file.
Assumes there is only one protein chain.
Written by Xingcheng Lin, modified by Yafan Zhang, 2025/05/07
"""

import sys
import mdtraj as md

def find_prot_chainID(pdbfile, output_file):
    # Load the PDB and convert topology to DataFrame
    pdb = md.load(pdbfile)
    df, _ = pdb.topology.to_dataframe()

    # Define standard amino acid residue names
    aa_residues = {
        'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS',
        'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP',
        'TYR','VAL'
    }

    # Filter for protein residues
    protein_df = df[df['resName'].isin(aa_residues)]

    # Get chain ID from first protein atom
    chain_id = protein_df.iloc[0]['chainID']

    # Convert numerical chain ID to letter 
    chain_letter = chr(ord('A') + int(chain_id))

    # Write to output file
    with open(output_file, 'w') as f:
        f.write(f"{chain_letter}\n")

if __name__ == "__main__":
    find_prot_chainID(sys.argv[1], sys.argv[2])
