###########################################################################
# This script will map DNA sequence from the output of Modeller to real
#
# Written by Xingcheng Lin, 12/12/2016;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################


def mapDNAseq_reverse(DNAseq_file, outputFile):

    # Dictionary between Modeller to Real DNA sequence;
    Real_to_Modeller = {
            "A": "e",
            "G": "l",
            "C": "j",
            "T": "t"
            }

    infile = open(DNAseq_file, 'r')
    outfile = open(outputFile, 'w')

    # Read in lines from the file;

    lines = [line.strip() for line in infile]

    infile.close()

    length = len(lines)

    for i in my_lt_range(0, length, 1):

        seq = list(lines[i])
        real_seq = []

        for idx in my_lt_range(0, len(seq), 1):
            real_seq.append(Real_to_Modeller[seq[idx]])

        converted_seq = ('').join(real_seq)

        outfile.write(converted_seq + "\n")

    return


############################################################################

if __name__ == "__main__":
    dnaSeq_file = sys.argv[1]
    outputFile = sys.argv[2]

    mapDNAseq_reverse(dnaSeq_file, outputFile)
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
