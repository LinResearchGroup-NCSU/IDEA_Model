import math
import subprocess
import os
import time
import sys
import functools
import itertools

import numpy as np
import random


# For Biopython
from Bio.PDB import *
#from Bio.PDB.Polypeptide import protein_letters_1to3, protein_letters_3to1
# For converting three-letter codes to one-letter codes
from Bio.Data.IUPACData import protein_letters_3to1

# For converting one-letter codes to three-letter codes
from Bio.Data.IUPACData import protein_letters_1to3
# Create a new dictionary that maps both uppercase and lowercase three-letter codes to one-letter codes
# Create a new dictionary that maps three-letter codes to one-letter codes, regardless of case
insensitive_protein_letters_3to1 = {}

# Populate the dictionary with uppercase, lowercase, and title case keys
for three_letter, one_letter in protein_letters_3to1.items():
    insensitive_protein_letters_3to1[three_letter.upper()] = one_letter
    insensitive_protein_letters_3to1[three_letter.lower()] = one_letter
    insensitive_protein_letters_3to1[three_letter.title()] = one_letter




phis_directory = "./phis/"

res_type_map = {
    'A': 0,
    'C': 4,
    'D': 3,
    'E': 6,
    'F': 13,
    'G': 7,
    'H': 8,
    'I': 9,
    'K': 11,
    'L': 10,
    'M': 12,
    'N': 2,
    'P': 14,
    'Q': 5,
    'R': 1,
    'S': 15,
    'T': 16,
    'V': 19,
    'W': 17,
    'Y': 18,
    'DA': 20,
    'DG': 21,
    'DC': 22,
    'DT': 23
}


res_name_map = {
    'e': ' DA',
    'l': ' DG',
    'j': ' DC',
    't': ' DT'
}


def save_structure(structure, file_name):
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_name)



def add_virtual_glycine_to_residue(residue):
        # get atom coordinates as vectors
    n = residue['N'].get_vector()
    c = residue['C'].get_vector()
    ca = residue['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis(-np.pi * 120.0 / 180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    atom = Atom.Atom("CB", cb, 0, 1, " ", " CB ", 0, element="CB")
    residue.add(atom)
    return residue



def is_hetero(residue):
    if residue.id[0] != ' ':
        return True
    else:
        return False


def get_res_list(structure, tm_only=False):
    tms_directory = "./tms/"
    pdb_id = structure.get_id().split('/')[-1]
    res_list = Selection.unfold_entities(structure, 'R')

    # Get all residues from a structure
    res_list = [residue for residue in res_list if not is_hetero(residue)]

    if tm_only:
        tm = read_column_from_file(os.path.join(
            tms_directory, pdb_id + '.tm'), 1)
        res_list = [residue for i, residue in enumerate(
            res_list) if tm[i] == '2']

    return res_list


def parse_pdb(pdb_id):
    parser = PDBParser()
    return parser.get_structure(pdb_id, "%s.pdb" % pdb_id)


def read_column_from_file(file_name, column, header_comment_syntax="#", num_header_lines=0, column_delimiter=''):
    list_to_return = []
    for i, line in enumerate(open(file_name, 'r')):
        line = line.strip('\n')
        if len(line) == 0:
            continue
        if i < num_header_lines:
            continue
        if line[0] == header_comment_syntax:
            continue
        if column_delimiter == '':
            line = line.split()
        else:
            line = line.split(column_delimiter)
        list_to_return.append(line[int(column) - 1])
    return list_to_return


def add_virtual_glycines(structure):
    residues = get_res_list(structure)
    for residue in residues:
        if residue.get_resname() == "GLY":
            residue = add_virtual_glycine_to_residue(residue)

    return structure


def add_virtual_glycines_list(proteins_list_file_name):
    proteins_list = read_column_from_file(proteins_list_file_name, 1)
    error_list_file = open("key_errors.dat", 'w')
    for protein in proteins_list:
        structure = parse_pdb(protein)
        try:
            structure = add_virtual_glycines(structure)
        except KeyError:
            error_list_file.write("%s\n" % protein)
            continue
        save_structure(structure, protein + '.pdb')


def read_phi_list(phi_list_file_name, header_comment_syntax="#", num_header_lines=0, column_delimiter=' '):
    input_file = open(phi_list_file_name, 'r')
    phi_list = []
    for i, line in enumerate(input_file):
        line = line.strip()
        if len(line) == 0:
            continue
        if i < num_header_lines:
            continue
        if line[0] == header_comment_syntax:
            continue
        line = line.split(column_delimiter, 1)
        try:
            parameters = line[1].split()
        except IndexError:
            parameters = []
        phi_list.append([line[0], parameters])
    return phi_list

def get_neighbor_list(structure, tm_only=False):
    protein = get_protein_name(structure)
    res_list = get_res_list(structure)
    atom_list = [a for a in get_atom_list(
        structure) if not is_hetero(a.get_parent())]
    if tm_only:
        tm = read_column_from_file(os.path.join(
            tms_directory, protein + '.tm'), 1)
        atom_list = [a for a in atom_list if tm[get_global_index(
            res_list, a.get_parent())] == '2']

    neighbor_list = NeighborSearch(atom_list)
    return neighbor_list

def get_protein_name(structure):
    return structure.get_id().split('/')[-1].split('.')[0]

def get_atom_list(structure):
    atom_list = Selection.unfold_entities(structure, 'A')  # A for atoms
    return atom_list

def get_sequence_from_structure(structure):
    sequence = ""
    ppb = PPBuilder(radius=10.0)
    for pp in ppb.build_peptides(structure, aa_only=False):
        sequence += '%s\n' % pp.get_sequence()
    return sequence.replace('\n', '')


def get_parameters_string(parameters):
    parameter_string = ""
    for parameter in parameters:
        parameter_string += parameter
        if not parameters.index(parameter) + 1 == len(parameters):
            parameter_string += '_'
    return parameter_string

def get_number_of_lines_in_file(filename):
    try:
        return open(filename, 'r').read().count("\n")
    except IOError:
        return 0

def get_local_index(residue):
    return residue.get_id()[1]

def get_chain(residue):
    return residue.get_parent().get_id()

def get_neighbors_within_radius(neighbor_list, residue, radius):
    #print(residue.resname)
    return neighbor_list.search(get_interaction_atom(residue).get_coord(), radius, level='R')

def get_interaction_atom(residue):
    try:
        if (residue.resname.strip() == "DA") or (residue.resname.strip() == "DT") or (residue.resname.strip() == "DC") or (residue.resname.strip() == "DG"):
           # print(residue)
            try:
                residue["C5"]
               # print(residue["P"])
                return residue["C5"]
            except KeyError:
                try: 
                    residue["O5'"]
                   # print(residue["O5'"])
                    return(residue["O5'"])
                except:
                    raise
        else:
            # Only use the CA atom here;
            if residue.resname == "GLY":
                return residue["CA"]
            else:
                return residue["CA"]
    except:
        raise


def get_global_index(residue_list, residue):
    return residue_list.index(residue)

def read_decoy_sequences(sequence_file_name):
    sequences = []
    with open(sequence_file_name, "r") as sequence_file:
        for line in sequence_file:
            line = line.strip()
            sequences.append(line)
    return sequences

def mutate_whole_sequence(res_list, new_sequence):
    for i in range(len(res_list)):
        # If it is a protein sequence:
        if new_sequence[i] in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            res_list[i].resname = protein_letters_1to3[new_sequence[i]]
        # If it is a DNA sequence:
        else:
            res_list[i].resname = res_name_map[new_sequence[i]]
    return res_list

def get_res_type(res_list, residue):
    # If it is DNA:
    if (residue.get_resname().strip() == "DA") or (residue.get_resname().strip() == "DT") or (residue.get_resname().strip() == "DC") or (residue.get_resname().strip() == "DG"):
       # print(residue)
        return res_type_map[residue.get_resname().strip()]
    # If it is protein:
    else: 
        return res_type_map[insensitive_protein_letters_3to1[residue.get_resname()]]

def get_interaction_distance(res1, res2):
    return get_interaction_atom(res1) - get_interaction_atom(res2)


def interaction_well(r, r_min, r_max, kappa):
    return 0.5 * (np.tanh(kappa * (r - r_min)) * np.tanh(kappa * (r_max - r))) + 0.5


def read_native_phi(protein, phi_list, total_phis, jackhmmer=False):
    phi_native = np.zeros(total_phis)
    i_phi = 0
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        parameters = phi_and_parameters[1]
        parameters_string = get_parameters_string(parameters)
        if jackhmmer:
            input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_native_%s" % (
                phi, protein, parameters_string)), 'r')
        else:
            input_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                phi, protein, parameters_string)), 'r')

        for line in input_file:
            line = line.strip().split()
            for i_value, value_i in enumerate(line):
                phi_native[i_phi] = float(line[i_value])
                i_phi += 1
    return phi_native

def read_decoy_phis(protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=False):
    phi_i_decoy = np.zeros((num_decoys, total_phis))

    for i_phi_function, phi_and_parameters in enumerate(phi_list):
        phi = phi_and_parameters[0]
        parameters = phi_and_parameters[1]
        i_phi = phi_list.index(phi_and_parameters)
        parameters_string = get_parameters_string(parameters)
        if jackhmmer:
            input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_decoys_%s" % (
                phi, protein, parameters_string)), 'r')
        else:
            input_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (
                phi, protein, decoy_method, parameters_string)), 'r')
        for i_decoy, line in enumerate(input_file):
            if i_decoy >= num_decoys:
                break
            first_phi = np.cumsum(num_phis)[
                i_phi_function] - num_phis[i_phi_function]
            i_phi = first_phi
            line = line.strip().split()
            for i_value, value_i in enumerate(line):
                phi_i_decoy[i_decoy][i_phi] = float(line[i_value])
                i_phi += 1
    return phi_i_decoy




def get_total_phis_and_parameter_string(phi_list, training_set):

 
    full_parameters_string = ""
    # Find out how many total phi_i there are
    total_phis = 0
    num_phis = []
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        full_parameters_string += phi
        parameters = phi_and_parameters[1]
        i_phi = phi_list.index(phi_and_parameters)
        parameters_string = get_parameters_string(parameters)
        full_parameters_string += parameters_string
        for i_protein, protein in enumerate(training_set):
            if i_protein > 0:
                break
            input_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                phi, protein, parameters_string)), 'r')
            for line in input_file:
                line = line.strip().split()
                num_phis.append(len(line))
                total_phis += len(line)
    return total_phis, full_parameters_string, num_phis



def sort_eigenvalues_and_eigenvectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors

# These three functions provide a way of calling a function multiple times (that run independently)
# on a certain number of processors so that a new function call starts when a processor becomes available


def call_independent_functions_on_n_processors(function, arguments_lists, num_processors):
    from multiprocessing import Pool
    pool = Pool(int(num_processors))
    results = pool.map(universal_worker, pool_args(function, *arguments_lists))


def universal_worker(input_pair):
    function, args = input_pair
    return function(*args)


def pool_args(function, *args):
    return list(zip(itertools.repeat(function), list(zip(*args))))


################################################################################################################

def get_total_phis_and_parameter_string_decoy_structures_provided(phi_list, training_set):
    full_parameters_string = ""
    # Find out how many total phi_i there are
    total_phis = 0
    num_phis = []
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        full_parameters_string += phi
        parameters = phi_and_parameters[1]
        i_phi = phi_list.index(phi_and_parameters)
        parameters_string = get_parameters_string(parameters)
        full_parameters_string += parameters_string
        for i_protein, protein in enumerate(training_set):
            if i_protein > 0:
                break
            input_file = open(os.path.join(phis_directory, "%s_%s_decoy_%s" % (
                phi, protein, parameters_string)), 'r')
            for line in input_file:
                line = line.strip().split()
                num_phis.append(len(line))
                total_phis += len(line)
    return total_phis, full_parameters_string, num_phis


def read_decoy_phi_structures_provided(protein, phi_list, total_phis, jackhmmer=False):
    phi_decoy = np.zeros(total_phis)
    i_phi = 0
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        parameters = phi_and_parameters[1]
        parameters_string = get_parameters_string(parameters)
        if jackhmmer:
            input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_decoy_%s" % (
                phi, protein, parameters_string)), 'r')
        else:
            input_file = open(os.path.join(phis_directory, "%s_%s_decoy_%s" % (
                phi, protein, parameters_string)), 'r')

        for line in input_file:
            line = line.strip().split()
            for i_value, value_i in enumerate(line):
                phi_decoy[i_phi] = float(line[i_value])
                i_phi += 1
    return phi_decoy



def read_all_gammas(phi_list_file_name, training_set_file, training_decoy_method, gamma_file_name=None, noise_filtering=True, read_confidence=False, bootstrapping_confidence=95, bootstrapping_iterations=1000, read_averaged_gammas=False, read_original_phis=False):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)

    # If we need to read in the cases where the decoy structures are explicitly provided, we need to change the name correspondingly;
    if read_original_phis == "decoy" and training_decoy_method == "TCR_modeling":
        total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string_decoy_structures_provided(
            phi_list, training_set)
    else:
        total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
            phi_list, training_set)
    if gamma_file_name == None:
        if noise_filtering:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma_filtered" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_averaged_gammas:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma_averaged" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_original_phis == "native":
            gamma_file_name = os.path.join(phis_directory, "%s_%s_phi_native_summary.txt" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_original_phis == "decoy":
            gamma_file_name = os.path.join(phis_directory, "%s_%s_phi_decoy_summary.txt" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        else:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
    if read_confidence:
        confidence_lower = np.loadtxt(os.path.join(gammas_directory, "%s_%s_confidence_lower_%d_%d" % (training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string, bootstrapping_confidence, bootstrapping_iterations)))
        confidence_upper = np.loadtxt(os.path.join(gammas_directory, "%s_%s_confidence_upper_%d_%d" % (training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string, bootstrapping_confidence, bootstrapping_iterations)))

    # Need to filter out the complex number if in the "filtered" mode;
    if noise_filtering:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    individual_gammas = []
    individual_confidence_lower = []
    individual_confidence_upper = []
    for i, (phi, parameters) in enumerate(phi_list):
        individual_gammas.append(gamma[0:num_phis[i]])
        np.savetxt(gammas_directory + 'individual_gamma' +
                   str(i) + '.txt', individual_gammas[i], fmt='%1.5f')
        gamma = gamma[num_phis[i]:]
        if read_confidence:
            individual_confidence_lower.append(confidence_lower[0:num_phis[i]])
            confidence_lower = confidence_lower[num_phis[i]:]
            individual_confidence_upper.append(confidence_upper[0:num_phis[i]])
            confidence_upper = confidence_upper[num_phis[i]:]

    if read_confidence:
        return individual_gammas, individual_confidence_lower, individual_confidence_upper
    else:
        return individual_gammas

