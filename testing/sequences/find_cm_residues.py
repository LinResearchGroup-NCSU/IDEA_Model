#!/usr/bin/env python3

"""
This script finds all protein–DNA residue contacts in a PDB file,
and writes the contacting protein and DNA residue IDs (1-indexed) to two output files.

Usage:
    python find_cm_residues.py <pdb_file> <cutoff> <out_protein.txt> <out_dna.txt>

Written by Xingcheng Lin, refactored by Yafan Zhang, 2025/05/07
"""

import sys
import numpy as np
import mdtraj as md
import itertools

# Define DNA residue names 
DNA_RESNAMES = {'DA', 'DG', 'DC', 'DT'}

def get_residue_indices(topology):
    dna_res = [res.index for res in topology.residues if res.name in DNA_RESNAMES]
    prot_res = [res.index for res in topology.residues if res.is_protein]
    return dna_res, prot_res

def find_contacts(pdbfile, cutoff):
    traj = md.load_pdb(pdbfile)
    dna_res, prot_res = get_residue_indices(traj.topology)

    if not dna_res or not prot_res:
        raise ValueError("No DNA or protein residues found.")

    # Build all DNA–protein residue pairs
    pairs = list(itertools.product(dna_res, prot_res))
    
    # Compute closest-heavy atom distances for each pair
    distances, _ = md.compute_contacts(traj, contacts=pairs, scheme='closest-heavy')
    distances = distances.flatten()

    # Filter by cutoff
    mask = distances < cutoff
    contact_pairs = np.array(pairs)[mask]

    # Convert to 1-indexed and deduplicate
    dna_ids = np.unique(contact_pairs[:, 0] + 1).astype(int)
    prot_ids = np.unique(contact_pairs[:, 1] + 1).astype(int)

    return prot_ids, dna_ids

def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)

    pdbfile = sys.argv[1]
    cutoff = float(sys.argv[2])
    out_protein = sys.argv[3]
    out_dna = sys.argv[4]

    prot_ids, dna_ids = find_contacts(pdbfile, cutoff)

    np.savetxt(out_protein, prot_ids[np.newaxis, :], fmt='%d')
    np.savetxt(out_dna, dna_ids[np.newaxis, :], fmt='%d')

if __name__ == "__main__":
    main()
