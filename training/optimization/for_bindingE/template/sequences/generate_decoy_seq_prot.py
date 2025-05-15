import os
import random
import sys
sys.path.append('../../../../common_functions')
from common_function import *

def load_gbinder_sequences(filepath):
    with open(filepath, 'r') as f:
        return set(f.read().splitlines())

def load_resids_to_randomize(filepath):
    with open(filepath, 'r') as f:
        return [int(x)-1 for x in f.readline().split()]

def generate_single_prot_decoy(sequence, positions):
    amino_acids = ["A", "R", "N", "D", "C", "E", "Q", "G",
                   "H", "I", "L", "K", "M", "F", "P", "S",
                   "T", "W", "Y", "V"]
    sequence_list = list(sequence)
    for pos in positions:
        sequence_list[pos] = random.choice(amino_acids)
    return ''.join(sequence_list)

def generate_protein_decoy_sequences(proteins_list_file_name, num_decoys=10000, randomSeed=None):
    random.seed(randomSeed)

    protein_list = read_column_from_file(proteins_list_file_name, 1)

    method = 'prot_randomization'
    os.makedirs(method, exist_ok=True)

    resids_toberandomized = load_resids_to_randomize("randomize_position_prot.txt")
    gBinder_sequences = load_gbinder_sequences("gBinder_sequences.txt")

    for protein in protein_list:
        print(f"Generating {num_decoys} protein decoys")

        with open(f"{method}/{protein}.decoys", 'w') as output_file:
            sequence_path = f"{method}/{protein}.seq"
            with open(sequence_path, "r") as seq_file:
                native_sequence = seq_file.read().replace('\n', '')

            decoys_generated = set()
            attempts = 0
            max_attempts = num_decoys * 10  
            
            while len(decoys_generated) < num_decoys and attempts < max_attempts:
                decoy_seq = generate_single_prot_decoy(native_sequence, resids_toberandomized)
                if decoy_seq not in gBinder_sequences and decoy_seq not in decoys_generated:
                    output_file.write(decoy_seq + '\n')
                    decoys_generated.add(decoy_seq)
                attempts += 1
            
            if len(decoys_generated) < num_decoys:
                print(f"Warning: Only generated {len(decoys_generated)} unique decoys for {protein} after {attempts} attempts.")

if __name__ == "__main__":
    generate_protein_decoy_sequences("proteins_list.txt", num_decoys=10000, randomSeed=0)