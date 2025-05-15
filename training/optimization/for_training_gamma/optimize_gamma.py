import sys

import numpy as np


sys.path.append('../../common_functions')
from common_function import *

def get_filtered_B_inv_lambda_and_P(filtered_lamb, cutoff_mode, P, method='extend_all_after_first_noisy_mode'):
    if method == 'zero_all_after_first_noisy_mode':
        filtered_lamb_inv = 1 / filtered_lamb
        # for "zeroing unreliable eigenvalues"
        filtered_lamb_inv[cutoff_mode:] = 0.0
        filtered_B_inv = np.dot(
            P, np.dot(np.diag(filtered_lamb_inv), np.linalg.inv(P)))
        filtered_lamb = 1 / filtered_lamb_inv
    if method == 'extend_all_after_first_noisy_mode':
        # for "extending lowest reliable eigenvalue"
        filtered_lamb[cutoff_mode:] = filtered_lamb[cutoff_mode - 1]
        filtered_B_inv = np.dot(
            P, np.dot(np.diag(1 / filtered_lamb), np.linalg.inv(P)))

    return filtered_B_inv, filtered_lamb, P

def get_filtered_gamma_B_lamb_P_and_lamb(A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys, noise_iterations=10, relative_error_threshold=0.5):
    lamb, P = np.linalg.eig(B)
    lamb, P = sort_eigenvalues_and_eigenvectors(lamb, P)

    cutoff_modes = []
    for i_noise in range(noise_iterations):
        noisy_B = np.zeros((total_phis, total_phis))
        for i in range(total_phis):
            for j in range(i, total_phis):
                random_B_ij = np.random.normal(
                    loc=half_B[i][j], scale=std_half_B[i][j] / float(num_decoys))
                noisy_B[i][j] = noisy_B[j][i] = random_B_ij - \
                    other_half_B[i][j]

        noisy_lamb, noisy_P = np.linalg.eig(noisy_B)
        noisy_lamb, noisy_P = sort_eigenvalues_and_eigenvectors(
            noisy_lamb, noisy_P)
        # print("np.abs(lamb - noisy_lamb) / lamb: " + str(np.abs(lamb - noisy_lamb) / lamb) )
        
        try:
            cutoff_mode = np.where(np.abs(lamb - noisy_lamb) / lamb > relative_error_threshold)[0][0]
        except IndexError:
            cutoff_mode = len(lamb)
        cutoff_modes.append(cutoff_mode)
    
    cutoff_mode = min(cutoff_modes)
    # Hard set a cutoff_mode by looking at the Lamb file itself;
    cutoff_mode = 70
    print(f"We retain the first {cutoff_mode} eigenvalues to filter gamma file.")

    filtered_lamb = np.copy(lamb)
    print(
    f"The last retained lambda after filtering "
    f"(index {cutoff_mode - 1} in filtered_lamb): {filtered_lamb[cutoff_mode - 1]:.5f}")

    filtered_B_inv, filtered_lamb, P = get_filtered_B_inv_lambda_and_P(
        filtered_lamb, cutoff_mode, P)

    filtered_gamma = np.dot(filtered_B_inv, A)
    filtered_B = np.linalg.inv(filtered_B_inv)
    return filtered_gamma, filtered_B, filtered_lamb, P, lamb

def calculate_A_and_B_wei(average_phi_decoy, phi_native, all_phis):
    # print("calculate_A_and_B")
    # print(datetime.datetime.now())
    A = average_phi_decoy - phi_native
    size_of_training_set, num_decoys, total_phis = all_phis.shape
    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    other_half_B = np.zeros((total_phis, total_phis))

    for p in range(size_of_training_set):
        phis_i = all_phis[p].reshape(num_decoys, total_phis, 1)
        phis_j = all_phis[p].reshape(num_decoys, 1, total_phis)
        half_B += np.average(phis_i * phis_j, axis=0)
        std_half_B += np.std(phis_i * phis_j, axis=0)
    half_B /= size_of_training_set
    std_half_B /= size_of_training_set

    for p in range(size_of_training_set):
        average_phi = np.average(all_phis[p], axis=0)
        other_half_B += average_phi.reshape(total_phis, 1) * average_phi.reshape(1, total_phis)
    other_half_B /= size_of_training_set
    # print("End")
    # print(datetime.datetime.now())

    B = half_B - other_half_B

    return A, B, half_B, other_half_B, std_half_B

def calculate_A_B_and_gamma_xl23(
    training_set_file,
    phi_list_file_name,
    decoy_method,
    num_decoys,
    noise_filtering=True,
    jackhmmer=False,
):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)

    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set
    )
    print(f"training set: {training_set}")

    phi_native_i_protein = np.zeros((len(training_set), total_phis))
    for i, protein in enumerate(training_set):
        phi_native_i_protein[i] = read_native_phi(
            protein, phi_list, total_phis, jackhmmer=jackhmmer
        )

    phi_native = np.average(phi_native_i_protein, axis=0)

    file_prefix = f"{phis_directory}{training_set_file.split('/')[-1].split('.')[0]}_{full_parameters_string}"
    np.savetxt(f"{file_prefix}_phi_native_summary.txt", phi_native, fmt="%1.5f")

    phi_i_protein_i_decoy = np.zeros((len(training_set), num_decoys, total_phis))
    for i, protein in enumerate(training_set):
        phi_i_protein_i_decoy[i] = read_decoy_phis(
            protein,
            phi_list,
            total_phis,
            num_phis,
            num_decoys,
            decoy_method,
            jackhmmer=jackhmmer,
        )

    phi_i_decoy = phi_i_protein_i_decoy.reshape(
        len(training_set) * num_decoys, total_phis
    )
    average_phi_decoy = np.average(phi_i_decoy, axis=0)
    np.savetxt(f"{file_prefix}_phi_decoy_summary.txt", average_phi_decoy, fmt="%1.5f")

    A, B, half_B, other_half_B, std_half_B = calculate_A_and_B_wei(
        average_phi_decoy, phi_native, phi_i_protein_i_decoy
    )
    gamma = np.linalg.pinv(B) @ A

    file_prefix = f"{gammas_directory}{training_set_file.split('/')[-1].split('.')[0]}_{full_parameters_string}"
    np.savetxt(f"{file_prefix}_gamma", gamma, fmt="%1.5f")
    np.savetxt(f"{file_prefix}_A", A, fmt="%1.5f")
    np.savetxt(f"{file_prefix}_B", B, fmt="%1.5f")

    if noise_filtering:
        (
            filtered_gamma,
            filtered_B,
            filtered_lamb,
            P,
            lamb,
        ) = get_filtered_gamma_B_lamb_P_and_lamb(
            A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys
        )

        np.savetxt(f"{file_prefix}_gamma_filtered", filtered_gamma, fmt="%1.5f")
        np.savetxt(f"{file_prefix}_B_filtered", filtered_B, fmt="%1.5f")
        np.savetxt(f"{file_prefix}_lamb_filtered", filtered_lamb, fmt="%1.5f")
        np.savetxt(f"{file_prefix}_P", P, fmt="%1.5f")
        np.savetxt(f"{file_prefix}_lamb", lamb, fmt="%1.5f")
        return

    return A, B, gamma

############################################

gammas_directory = "./gammas/randomized_decoy/"

calculate_A_B_and_gamma_xl23("native_trainSetFiles.txt", "phi1_list.txt", decoy_method='CPLEX_randomization', 
                             num_decoys=11000, noise_filtering=True, jackhmmer=False)