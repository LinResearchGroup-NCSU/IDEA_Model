import time
import subprocess
import os
import math
import sys
import numpy as np

def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step

def create_tms(random_position_file, tot_resnum, PDBid):

    tms_directory = "./tms/"

    infile = open(random_position_file, 'r')
    lines = [line.strip() for line in infile]
    infile.close()

    random_position = lines[0].split(" ")
    random_position_int = [int(integer) for integer in random_position]
    # print(random_position_int)

    tms_content = []

    for i in my_le_range(1, tot_resnum, 1):
        if i in random_position_int:
            tms_content.append("2")
        else:
            tms_content.append("1")

    with open(tms_directory + PDBid + '.tm', 'w') as outfile:
        outfile.write("\n".join(tms_content) + "\n")

    with open(tms_directory + 'native.tm', 'w') as outfile:
        outfile.write("\n".join(tms_content) + "\n")

    return

if __name__ == "__main__":
    random_position_file = sys.argv[1]
    tot_resnum = int(sys.argv[2])
    PDBid = sys.argv[3]
    print(f"The total number of residues for {PDBid} is {tot_resnum}.")
    create_tms(random_position_file, tot_resnum, PDBid)