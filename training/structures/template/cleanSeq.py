###########################################################################
# This script will remove the "\n" of the Modeller generated sequence
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


def cleanSeq(inputFile, outputFile):

    infile = open(inputFile, 'r')
    outfile = open(outputFile, 'w')

    # Read in lines from the file;
    lines = [line.strip() for line in infile]
    infile.close()

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
            else:
                outfile.write(lines[i])


############################################################################

if __name__ == "__main__":

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]

    seq = cleanSeq(inputFile, outputFile)
    
    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
