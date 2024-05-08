###########################################################################
# This script will copy the cleaned template sequence and list the template
# sequences for template_alignment.ali
#
# Written by Xingcheng Lin, 03/13/2021;
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


def copy_seq_fortemplate(inputFile, outputFile, pdb_id):

    infile = open(inputFile, 'r')
    cdr3afile = open('self_aligned_cdr3a', 'r')
    cdr3bfile = open('self_aligned_cdr3b', 'r')
    pepfile = open('self_aligned_peptide', 'r')
    # Can only do append, because there might be multiple templates
    outfile = open(outputFile, 'a')

    # Read in lines from the file;
    lines = [line.strip() for line in cdr3afile]
    cdr3afile.close()

    length = len(lines)

    for i in my_lt_range(0, length, 1):
        line = lines[i].split()
        if (line[0] == pdb_id):
            cdr3a_seq = line[1]
            cdr3a_seq_orig = line[1].replace("-", "")

    # Read in lines from the file;
    lines = [line.strip() for line in cdr3bfile]
    cdr3bfile.close()

    length = len(lines)

    for i in my_lt_range(0, length, 1):
        line = lines[i].split()
        if (line[0] == pdb_id):
            cdr3b_seq = line[1]
            cdr3b_seq_orig = line[1].replace("-", "")

    # Read in lines from the file;
    lines = [line.strip() for line in pepfile]
    pepfile.close()

    length = len(lines)

    for i in my_lt_range(0, length, 1):
        line = lines[i].split()
        if (line[0] == pdb_id):
            pep_seq = line[1]
            pep_seq_orig = line[1].replace("-", "")


    # Find the position of the sequences in the given template string
    lines = [line.strip() for line in infile]
    infile.close()

    # The template sequence should be in the third line
    template_sequence = lines[2]

    # Replace the original CDR3A/CDR3B/peptide sequence with the newly aligned sequences;
    tmp_seq = template_sequence.replace(cdr3a_seq_orig, cdr3a_seq)
    tmp_seq = tmp_seq.replace(cdr3b_seq_orig, cdr3b_seq)
    aligned_template_sequence = tmp_seq.replace(pep_seq_orig, pep_seq)


    length = len(lines)

    for i in my_lt_range(0, length, 1):
        line = lines[i].split()
        try:
            line[0]
        except IndexError:
            continue
        else:
            line_letters = list(line[0])
            if (line_letters[0] == ">"):
                # This is the first line
                outfile.write(lines[i] + "\n")
            elif (line_letters[3] == "u"):
                # This is the 2nd line
                outfile.write(lines[i] + "\n")

    # Output the aligned template sequence
    outfile.write(aligned_template_sequence + "\n")
    outfile.close()
    

############################################################################

if __name__ == "__main__":
    # Cleaned template sequence
    inputFile = sys.argv[1]
    # template_alignment.ali
    outputFile = sys.argv[2]
    # PDB id of the corresponding template
    pdb_id = sys.argv[3]

    copy_seq_fortemplate(inputFile, outputFile, pdb_id)
    
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
