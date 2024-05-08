####################################################################################
# This script will use the script from the common_function.py to evaluate the
# binding energies of strong and decoy binders
#
#
# Written by Xingcheng Lin, 03/15/2021
####################################################################################

import math
import subprocess
import os
import math
import numpy as np
import sys
import time

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
###########################################


def validate_hamiltonian(hamiltonian, training_set_file, training_decoy_method, num_decoys, test_set_file=None, test_decoy_method=None, use_filtered_gammas=False):
    if test_set_file == None:
        test_set_file = training_set_file
    if test_decoy_method == None:
        test_decoy_method = training_decoy_method
    test_set = read_column_from_file(test_set_file, 1)
    z_scores = []
    e_natives = []
    e_mgs = []
    e_mg_stds = []
    e_decoys = []
    for protein in test_set:
        z, en, emg, emgstd, e_decoy = evaluate_hamiltonian(
            protein, hamiltonian, training_set_file, training_decoy_method, test_decoy_method, num_decoys, use_filtered_gammas)
        if np.isnan(z):
            continue
        z_scores.append(z)
        e_natives.append(en)
        e_mgs.append(emg)
        e_mg_stds.append(emgstd)
        e_decoys.append(e_decoy)

        # Output the randomized sequences decoy energies to file;
        # np.savetxt('randomized_decoy_energies.txt', e_decoys)

    return z_scores, e_natives, e_mgs, e_mg_stds, e_decoys

def evaluate_hamiltonian(protein, hamiltonian, training_set_file, training_decoy_method, test_decoy_method, num_decoys, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    training_set = read_column_from_file(training_set_file, 1)
    # read in Hamiltonian
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set)
    # read in corresponding gammas
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # read in corresponding phis (native and decoys)
    phi_native = read_native_phi(protein, phi_list, total_phis)
    phi_i_decoy = read_decoy_phis(
        protein, phi_list, total_phis, num_phis, num_decoys, test_decoy_method)
    # perform dot products to get energies (native and decoys)
    e_decoy = np.zeros(num_decoys)
    e_native = np.dot(gamma, phi_native)


    for i_decoy in range(num_decoys):
        e_decoy[i_decoy] = np.dot(gamma, phi_i_decoy[i_decoy])
        
    e_mg = np.average(e_decoy)
    e_mg_std = np.std(e_decoy)
    # calculate z-score
    z_score = (e_mg - e_native) / e_mg_std
    return z_score, e_native, e_mg, e_mg_std, e_decoy


############################################################################

if __name__ == "__main__":

    sys.path.append('../../../common_functions')
    from common_function import *

    # Gamma directory
    gammas_directory = "./gammas/randomized_decoy/"
    
    # Read in the target PDB
    infile = open("proteinList.txt", "r")
    
    # Read in lines from the file;
    lines = [line.strip() for line in infile]
    infile.close()
    
    length = len(lines)
    
    for i in my_lt_range(0, length, 1):
    
        # Create the target list for calculation of binding energies
        outfile = open("native_trainSetFiles.txt", "w")
        outfile.write(lines[i] + "\n")
        outfile.close()
    
        z_scores, e_natives, e_mgs, e_mg_stds, e_decoys = validate_hamiltonian(
            hamiltonian='phi1_list.txt', training_set_file ='native_trainSetFiles.txt', training_decoy_method='CPLEX_randomization', 
            num_decoys=1000, use_filtered_gammas=True)
    
    
        # Save the energies
        np.savetxt('enative.' + str(i) + '.txt', np.real(e_natives), fmt='%1.4f')
        np.savetxt('emg.' + str(i) + '.txt', np.real(e_decoys), fmt='%1.4f')
        np.savetxt('zscore.' + str(i) + '.txt', np.real(z_scores), fmt='%1.4f')
        
