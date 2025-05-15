###########################################################################
# This script will map DNA sequence to its reverse complement
#
# Written by Xingcheng Lin, 06/04/2023;
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

def reverse_complement(DNAseq_file, outputFile):

    # Dictionary between DNA and its complement;
    seq_to_complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
            }

    infile = open(DNAseq_file, 'r')
    outfile = open(outputFile, 'w')

    # Read in lines from the file;

    lines = [line.strip() for line in infile]

    infile.close()

    length = len(lines)

    for i in my_lt_range(0, length, 1):

        seq = list(lines[i])
        reverse_complement_seq = []

        for idx in my_lt_range(0, len(seq), 1):
            # note Python is 0-based;
            reversed_idx = len(seq) - idx - 1

            reverse_complement_seq.append(seq_to_complement[seq[reversed_idx]])

        converted_seq = ('').join(reverse_complement_seq)

        outfile.write(converted_seq + "\n")

    return


############################################################################

if __name__ == "__main__":
    dnaSeq_file = sys.argv[1]
    outputFile = sys.argv[2]

    reverse_complement(dnaSeq_file, outputFile)
    # print("When the voice of the Silent touches my words,")
    # print("I know him and therefore know myself.")
