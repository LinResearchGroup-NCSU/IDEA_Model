###########################################################################
# This script will the residue and nucleotide IDs that are in contacts
#
# Written by Xingcheng Lin, 12/31/2022;
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


def find_cm_residues(pdbfile, cutoff, random_position_file_protein, random_position_file_DNA):

    # Load the PDB
    pdb = md.load_pdb(pdbfile)
    # Convert into a pandas table
    pdb_table, pdb_bonds = pdb.topology.to_dataframe()

 #   print(pdb_table)

    # Select residue indices belonging to DNA
    DNA_resID = pdb_table.loc[(pdb_table['resName'] == 'DA') |  (pdb_table['resName'] == 'DC') | (pdb_table['resName'] == 'DT') | (pdb_table['resName'] == 'DG')].drop_duplicates(subset=['resSeq']).resSeq.to_numpy()
    # Select residue indices belonging to protein
    prot_resID = pdb_table.loc[(pdb_table['resName'] == 'ALA') |  (pdb_table['resName'] == 'ARG') | (pdb_table['resName'] == 'ASN') \
         | (pdb_table['resName'] == 'ASP') | (pdb_table['resName'] == 'CYS') | (pdb_table['resName'] == 'GLU') \
            | (pdb_table['resName'] == 'GLN') | (pdb_table['resName'] == 'GLY') | (pdb_table['resName'] == 'HIS') \
                | (pdb_table['resName'] == 'ILE') | (pdb_table['resName'] == 'LEU') | (pdb_table['resName'] == 'LYS') \
                    | (pdb_table['resName'] == 'MET') | (pdb_table['resName'] == 'PHE') | (pdb_table['resName'] == 'PRO') \
                        | (pdb_table['resName'] == 'SER') | (pdb_table['resName'] == 'THR') | (pdb_table['resName'] == 'TRP') \
                            | (pdb_table['resName'] == 'TYR') | (pdb_table['resName'] == 'VAL')].drop_duplicates(subset=['resSeq']).resSeq.to_numpy()
    # For mdtraj requirement, made it into 0-indexed
    DNA_resID_0_indexed = DNA_resID - 1
    prot_resID_0_indexed = prot_resID - 1

    # Select the atom indices belong to DNA
    DNA_atom_indices = pdb.topology.select("resname =~ 'D[ATCG]'")
    # Select the atom indices belong to protein
    prot_atom_indices = pdb.topology.select("is_protein == True")

    # Use mdtraj to calculate the contacts
    prot_DNA_respairs = list(itertools.product(DNA_resID_0_indexed, prot_resID_0_indexed))
 #   print(prot_DNA_respairs)
    prot_DNA_respair_distances, residue_pairs = md.compute_contacts(pdb, prot_DNA_respairs, scheme='closest-heavy')

  #  print(prot_DNA_respair_distances)
  #  print(residue_pairs)

    # Select those pairs that are closer than the cutoff
    prot_DNA_respair_distances = prot_DNA_respair_distances.flatten()
    cm_pair_idx = np.argwhere(prot_DNA_respair_distances < cutoff).flatten()

    # Recover the 1-indexed protein and DNA residue IDs that are in contact
    prot_cm_resID = []
    DNA_cm_resID = []
    for i in cm_pair_idx:
        DNA_cm_resID = np.append(DNA_cm_resID, prot_DNA_respairs[i][0] + 1)
        prot_cm_resID = np.append(prot_cm_resID, prot_DNA_respairs[i][1] + 1)

    # Remove the duplicated residue IDs
    prot_cm_resID = np.unique(prot_cm_resID).astype(int)
    prot_cm_resID = np.reshape(prot_cm_resID, (1, np.shape(prot_cm_resID)[0]))
    DNA_cm_resID = np.unique(DNA_cm_resID).astype(int)
    DNA_cm_resID = np.reshape(DNA_cm_resID, (1, np.shape(DNA_cm_resID)[0]))
 #   print(prot_cm_resID)
 #   print(DNA_cm_resID)

    np.savetxt(random_position_file_protein, prot_cm_resID, fmt='%d')
    np.savetxt(random_position_file_DNA, DNA_cm_resID, fmt='%d')

    return

############################################################################

if __name__ == "__main__":
    pdbfile = sys.argv[1]
    # Cutoff for determining residue contacts
    cutoff = float(sys.argv[2])
    # files for recording indicies
    random_position_file_protein = sys.argv[3]
    random_position_file_DNA = sys.argv[4]

    find_cm_residues(pdbfile, cutoff, random_position_file_protein, random_position_file_DNA)
 
    print("I slept and dreamt that life was joy. I awoke and saw that life was service. I acted and behold, service was joy.")
