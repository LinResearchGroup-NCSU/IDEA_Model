import sys
import re

# DNA mapping
dna_res_to_char = {
    'DA': 'e',
    'DG': 'l',
    'DC': 'j',
    'DT': 't'
}

# Protein mapping
aa_res_to_char = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'HYP': 'P'
}

def extract_sequence_from_pdb(pdb_file):
    seen = set()
    sequence = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain_id = line[21]
            resseq = line[22:26].strip()
            unique_res_id = (chain_id, resseq)

            if unique_res_id in seen:
                continue
            seen.add(unique_res_id)

            if resname in aa_res_to_char:
                sequence.append(aa_res_to_char[resname])
            elif resname in dna_res_to_char:
                sequence.append(dna_res_to_char[resname])
            # Unknown residues are skipped

    return ''.join(sequence)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)

    pdb_code = sys.argv[1]
    pdb_path = pdb_code + ".pdb"
    sequence = extract_sequence_from_pdb(pdb_path)

    with open(pdb_code + ".seq", "w") as f:
        f.write(sequence + "\n")