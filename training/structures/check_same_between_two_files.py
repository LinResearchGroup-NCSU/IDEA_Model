###########################################################################
# This script will check if two files have the same column contents and output
# the common contents into a new file

# Written by Xingcheng Lin, 11/06/2021;
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

def check_same_between_two_files(firstFile, secondFile, outputFile):

    first_infile = open(firstFile, 'r')
    second_infile = open(secondFile, 'r')
    outfile = open(outputFile, 'w')

    # Read in lines from the file;

    lines1 = [line.strip() for line in first_infile]
    lines2 = [line.strip() for line in second_infile]
    first_infile.close()
    second_infile.close()

    length1 = len(lines1)
    length2 = len(lines2)

    for i in my_lt_range(0, length1, 1):
        line1 = lines1[i].split()
        for j in my_lt_range(0, length2, 1):
            line2 = lines2[j].split()
            if (line1[0]==line2[0] and line1[1]==line2[1]):
                outfile.write(line1[0] + "\t" + line1[1] + "\n")

    outfile.close()

    return

############################################################################

if __name__ == "__main__":
    firstFile = sys.argv[1]
    secondFile = sys.argv[2]
    outputFile = sys.argv[3]

    check_same_between_two_files(firstFile, secondFile, outputFile)
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
