####################################################################################
# This script will find the structure with the smallest Z-score 
# and output that into the template_list.txt file for the next round of optimization
#
#
# Written by Xingcheng Lin, 03/23/2021
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


############################################################################

if __name__ == "__main__":

    num_of_cases = 60

    # Loop over all the zscore files
    for i in my_lt_range(0, num_of_cases, 1):
        z_score = np.loadtxt("results/zscore." + str(i) + ".txt")

        if i==0:
            z_scores = [z_score]
        else:
            z_scores.append(z_score)

    z_scores = np.asarray(z_scores)

    # Name of TCRs
    TCR_names = ['1ao7', '1bd2', '1lp9', '1oga', '1qrn', '2bnq', '2bnr', '2f53', '2f54', '2gj6', '2j8u',
             '2jcc', '2p5e', '2p5w', '2pye', '2uwe', '2vlj', '2vlk', '2vlr', '3d39', '3d3v', '3gsn',
             '3h9s', '3hg1', '3o4l', '3pwp', '3qdg', '3qdj', '3qdm', '3qeq', '3qfj', '3utt', '4ftv',
             '4jfd', '4jfe', '4jff', '4l3e', '4mnq', '4qok', '5c07', '5c08', '5c09', '5c0a', '5c0b',
             '5c0c', '5d2l', '5d2n', '5e6i', '5e9d', '5eu6', '5euo', '5hhm', '5hho', '5hyj', '5isz',
             '5men', '5nme', '5nmf', '5nmg', '5tez']    

    # We will just flatten the indices array with the minimum Z scores, and get the first element of that (in case there are cases with identical Z scores)
    smallest_z_score_idx = np.argwhere(z_scores == np.amin(z_scores)).flatten().item(0)
    smallest_z_score_label = TCR_names[smallest_z_score_idx]
    print(smallest_z_score_label)
    
    outfile = open("template_list.txt", "a")
    outfile.write(smallest_z_score_label + "\n")
    outfile.close()
    
