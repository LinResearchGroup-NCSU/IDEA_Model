from visualize import *
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib

sys.path.append('../../common_functions')
from common_function import *

phis_directory = "./phis/"
gammas_directory = "./gammas/randomized_decoy/"

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
            input_file = open(os.path.join(phis_directory,"%s_%s_native_%s" % (
                phi, protein, parameters_string)), 'r')
            for line in input_file:
                line = line.strip().split()
                num_phis.append(len(line))
                total_phis += len(line)
    return total_phis, full_parameters_string, num_phis

def get_parameters_string(parameters):
    parameter_string = ""
    for parameter in parameters:
        parameter_string += parameter
        if not parameters.index(parameter) + 1 == len(parameters):
            parameter_string += '_'
    return parameter_string

def plot_all_gammas_protDNA(phi_list_file_name, individual_gammas, vmin=-0.3, vmax=0.3, invert_sign=False, gammas_to_plot=None, plot_confidence=False, individual_confidence_lower=None, individual_confidence_upper=None, save_path=None, save_prefix="plot"):
    phi_list = read_phi_list(phi_list_file_name)
    if gammas_to_plot == None:
        gammas_to_plot = list(range(len(phi_list)))
    for i, (phi, parameters) in enumerate(phi_list):
        if not i in gammas_to_plot:
            continue
        plot_phi = globals()['plot_protDNA_' + phi]
        # print(phi, parameters)
        if plot_confidence:
            plot_phi(phi, parameters, individual_gammas[i], plot_confidence=plot_confidence,
                     confidence_lower=individual_confidence_lower[i], confidence_upper=individual_confidence_upper[i])
        else:
            plot_phi(individual_gammas[i], vmin=vmin,
                     vmax=vmax, invert_sign=invert_sign, save_path = save_path, save_prefix = save_prefix)

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


def plot_all_gammas(phi_list_file_name, individual_gammas, vmin=-0.3, vmax=0.3, invert_sign=False, gammas_to_plot=None, plot_confidence=False, individual_confidence_lower=None, individual_confidence_upper=None):
    phi_list = read_phi_list(phi_list_file_name)
    if gammas_to_plot == None:
        gammas_to_plot = list(range(len(phi_list)))
    for i, (phi, parameters) in enumerate(phi_list):
        if not i in gammas_to_plot:
            continue
        plot_phi = globals()['plot_' + phi]
        # print(phi, parameters)
        if plot_confidence:
            plot_phi(phi, parameters, individual_gammas[i], plot_confidence=plot_confidence,
                     confidence_lower=individual_confidence_lower[i], confidence_upper=individual_confidence_upper[i])
        else:
            plot_phi(individual_gammas[i], vmin=vmin,
                     vmax=vmax, invert_sign=invert_sign)

            
def plot_protDNA_phi_pairwise_contact_well(gammas, invert_sign=True, fix_colorbar=True, vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0, confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None, save_path=None, save_prefix="plot"):
    size = 24
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            # inverse_res_type_map maps each index to the residue according to the res_type_map, which is the map we used to record interaction matrix when doing the optimization
            # hydrophobicity_map used a new map and map the converted residue to new indices, that is for ease of creating a new heat map ordering residues based on their hydrophobic propensity
            index1 = hydrophobicity_map_protDNA[inverse_res_type_map_protDNA[i]]
            index2 = hydrophobicity_map_protDNA[inverse_res_type_map_protDNA[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111)
    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1
    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="coolwarm")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    # Move the X-axis ticks to the top of the figure, for better visualization 
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters_protDNA)
    ax.set_yticklabels(hydrophobicity_letters_protDNA)
    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)

    # plt.savefig('./plots/direct_contact.pdf')
    if save_path:
        # plt.savefig(os.path.join(save_path, f"{save_prefix}_gamma_plot.png"))
        matplotlib.rcParams['pdf.fonttype'] = 42
        plt.savefig(os.path.join(save_path, f"{save_prefix}.pdf"))
    # plt.show()
    
res_type_map_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

inverse_res_type_map = dict(list(zip(list(range(20)), res_type_map_letters)))

res_type_map_letters_protDNA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'DA', 'DG', 'DC', 'DT']

inverse_res_type_map_protDNA = dict(list(zip(list(range(24)), res_type_map_letters_protDNA)))

hydrophobicity_letters = ['R', 'K', 'N', 'Q', 'D', 'E', 'H', 'Y',
                          'W', 'S', 'T', 'G', 'P', 'A', 'M', 'C', 'F', 'L', 'V', 'I']

hydrophobicity_map = dict(list(zip(hydrophobicity_letters, list(range(20)))))

hydrophobicity_letters_protDNA = ['R', 'K', 'N', 'Q', 'D', 'E', 'H', 'Y',
                          'W', 'S', 'T', 'G', 'P', 'A', 'M', 'C', 'F', 'L', 'V', 'I', 'DA', 'DG', 'DC', 'DT']

hydrophobicity_map_protDNA = dict(list(zip(hydrophobicity_letters_protDNA, list(range(24)))))


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
    'Y': 18
}

res_type_map_protDNA = {
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
    ' DA': 20,
    ' DG': 21,
    ' DC': 22,
    ' DT': 23
}



def main():
    save_path = "./visualize/"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    individual_gammas_randomized_decoy = read_all_gammas("phi1_list.txt", "native_trainSetFiles.txt", training_decoy_method="CPLEX_randomization", noise_filtering=True)
    plot_all_gammas_protDNA('phi1_list.txt', individual_gammas_randomized_decoy, vmin=-10, vmax=10, invert_sign=False, save_path=save_path, save_prefix="trained_gamma")

    original_phis_native=read_all_gammas("phi1_list.txt", "native_trainSetFiles.txt", training_decoy_method="CPLEX_randomization", noise_filtering=False, read_original_phis="native")
    plot_all_gammas_protDNA('phi1_list.txt', original_phis_native, vmin=-3, vmax=3, invert_sign=False, save_path=save_path, save_prefix="native_phi")

    original_phis_decoy=read_all_gammas("phi1_list.txt", "native_trainSetFiles.txt", training_decoy_method="CPLEX_randomization", noise_filtering=False, read_original_phis="decoy")
    plot_all_gammas_protDNA('phi1_list.txt', original_phis_decoy, vmin=-3, vmax=3, invert_sign=False, save_path=save_path, save_prefix="decoy_phi")
    print("The visualized trained energy model is saved in the 'visualize' folder.")
if __name__ == "__main__":
    main()