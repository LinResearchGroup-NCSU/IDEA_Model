import os
import random
import sys
sys.path.append('../../../../common_functions')
from common_function import *

def shuffle_string(string):
    string_list = list(string)
    random.shuffle(string_list)
    return ''.join(string_list)

def cyclically_permute_string(string, degree):
    return string[degree:] + string[:degree]

def get_sublist(lst, indices):
    return [lst[i] for i in indices]

def load_gbinder_sequences(filepath):
    with open(filepath, 'r') as f:
        return set(f.read().splitlines())

def load_resids_to_randomize(filepath):
    with open(filepath, 'r') as f:
        return [int(x)-1 for x in f.readline().split()]

def generate_single_decoy(sequence, positions, method):
    sequence_list = list(sequence)
    if method == 'DNA_randomization':
        bases = ["e", "l", "j", "t"]
    elif method == 'prot_randomization':
        bases = ["A", "R", "N", "D", "C", "E", "Q", "G",
                 "H", "I", "L", "K", "M", "F", "P", "S",
                 "T", "W", "Y", "V"]
    else:
        raise ValueError("Unsupported method.")
    for pos in positions:
        sequence_list[pos] = random.choice(bases)
    return ''.join(sequence_list)

def generate_decoy_sequences(proteins_list_file_name, methods=['DNA_randomization'], num_decoys=[1000], randomSeed=None):
    random.seed(randomSeed)

    protein_list = read_column_from_file(proteins_list_file_name, 1)

    for method, n_decoys in zip(methods, num_decoys):
        os.makedirs(method, exist_ok=True)

        resids_toberandomized = load_resids_to_randomize("randomize_position_DNA.txt" if method == 'DNA_randomization' else "randomize_position_prot.txt")
        gBinder_sequences = load_gbinder_sequences("gBinder_sequences.txt")

        for protein in protein_list:
            print(f"Generating {n_decoys} DNA decoys")

            with open(f"{method}/{protein}.decoys", 'w') as output_file:
                sequence_path = f"{method}/{protein}.seq"
                with open(sequence_path, "r") as seq_file:
                    native_sequence = seq_file.read().replace('\n', '')

                decoys_generated = set()
                attempts = 0
                max_attempts = n_decoys * 10  
                
                while len(decoys_generated) < n_decoys and attempts < max_attempts:
                    decoy_seq = generate_single_decoy(native_sequence, resids_toberandomized, method)
                    if decoy_seq not in gBinder_sequences and decoy_seq not in decoys_generated:
                        output_file.write(decoy_seq + '\n')
                        decoys_generated.add(decoy_seq)
                    attempts += 1
                
                if len(decoys_generated) < n_decoys:
                    print(f"Warning: Only generated {len(decoys_generated)} unique decoys for {protein} after {attempts} attempts.")

if __name__ == "__main__":
    generate_decoy_sequences("proteins_list.txt", methods=['DNA_randomization'], num_decoys=[1000], randomSeed=0)