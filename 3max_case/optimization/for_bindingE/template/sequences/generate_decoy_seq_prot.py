####################################################################################
# This script will generate decoy sequences based on your given native sequence
#
# Written by Xingcheng Lin, 06/10/2020
####################################################################################

import math
import subprocess
import os
import time
import sys

import numpy as np

sys.path.append('../../../../common_functions')
from common_function import *

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step


def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
###########################################


def generate_decoy_sequences(proteins_list_file_name, methods=['DNA_randomization'], num_decoys=[1000], randomSeed=None):
    protein_list = read_column_from_file(proteins_list_file_name, 1)
    
    decoys_root_directory = "./"
    
    os.chdir(decoys_root_directory)
    for i, method in enumerate(methods):
        if not os.path.exists(method):
            os.makedirs(method)
        else:
            pass

        os.chdir(method)
        for protein in protein_list:
            print(method, protein)
            output_file = open("%s.decoys" % protein, 'w')
            
            # A random seed is provided if necessary for reproductibility of each protein
            random.seed(randomSeed)

            for j in range(num_decoys[i]):
                output_file.write(generate_decoy_sequence(protein, method=method, degree=j) + '\n')
            output_file.close()
        os.chdir('..')

def generate_decoy_sequence(protein, method='DNA_randomization', degree=None):

    sequences_root_directory = "../"

    with open("%s%s.seq" % (sequences_root_directory + method + '/', protein), "r") as sequence_file:
        native_sequence = sequence_file.read().replace('\n', '')

    if method == 'shuffle':
        return shuffle_string(native_sequence)
    elif method == 'cyclic':
        if degree == None:
            print("Must specify degree with method cyclic")
            sys.exit()
        return cyclically_permute_string(native_sequence, degree)
    elif method == 'constrained_shuffle' or method == 'constrained_cyclic':
        native_sequence = list(native_sequence)
        tm = read_column_from_file(os.path.join(
            tm_root_directory, protein + '.tm'), 1)
        outside_indices = [i for i, x in enumerate(tm) if x == '1' or x == '3']
        outside_list = get_sublist(native_sequence, outside_indices)
        outside_string = ''.join(outside_list)
        membrane_indices = [i for i, x in enumerate(tm) if x == '2']
        membrane_list = get_sublist(native_sequence, membrane_indices)
        membrane_string = ''.join(membrane_list)
        new_sequence = []
        if method == 'constrained_shuffle':
            outside_string = shuffle_string(outside_string)
            membrane_string = shuffle_string(membrane_string)
        elif method == 'constrained_cyclic':
            if degree == None:
                print("Must specify degree with method cyclic")
                sys.exit()
            outside_string = cyclically_permute_string(outside_string, degree)
            membrane_string = cyclically_permute_string(
                membrane_string, degree)
        outside_list = list(outside_string)
        membrane_list = list(membrane_string)
        for region in tm:
            if region == '1' or region == '3':
                new_sequence.append(outside_list.pop(0))
            if region == '2':
                new_sequence.append(membrane_list.pop(0))
        return ''.join(new_sequence)
    elif method == 'DNA_randomization':
        sequence_toberandomized = list(native_sequence)
        resids_toberandomized = open(
            "randomize_position_DNA.txt", 'r').readline().split(' ')
        for resid_toberandomized in resids_toberandomized:
            # Residue to be randomly mutated to (randomly selected from 4 nucleotides);
            resAbbr = random.choice(list(["e", "l", "j", "t"]))

            # Replace the corresponding nucleotide with 4 possibilities;
            sequence_toberandomized[int(resid_toberandomized) - 1] = resAbbr

        newsequence = ''.join(sequence_toberandomized)
        # Check and exclude those decoy sequences that coincide with the gBinder sequences, because we don't want good Binders
        # to be treated as weak binders in our training processes;
        gBinder_sequences = open(
            "gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return
        else:
            return newsequence
    elif method == 'prot_randomization':
        sequence_toberandomized = list(native_sequence)
        resids_toberandomized = open(
            "randomize_position_prot.txt", 'r').readline().split(' ')
        for resid_toberandomized in resids_toberandomized:
            # Residue to be randomly mutated to (randomly selected from 20 amino acids);
            resAbbr = random.choice(list(
                ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]))

            # Replace the corresponding amino acid with 20 possibilities;
            sequence_toberandomized[int(resid_toberandomized) - 1] = resAbbr

        newsequence = ''.join(sequence_toberandomized)
        # Check and exclude those decoy sequences that coincide with the gBinder sequences, because we don't want good Binders
        # to be treated as weak binders in our training processes;
        gBinder_sequences = open(
            "gBinder_sequences.txt", 'r').read().splitlines()
        if newsequence in gBinder_sequences:
            return
        else:
            return newsequence


############################################


generate_decoy_sequences("proteins_list.txt", methods=['prot_randomization'], num_decoys=[10000], randomSeed=0)
