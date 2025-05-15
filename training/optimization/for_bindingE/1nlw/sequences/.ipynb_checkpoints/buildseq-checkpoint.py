# This script will get the argument passed by command line
import sys
from modeller import Environ, Model, Alignment

# Get the sequence of the PDB file and write to an alignment file
code = sys.argv[1]

# Initialize Modeller environment
env = Environ()

# Load the model from the PDB file
mdl = Model(env, file=code)

# Create an alignment object and add the model
aln = Alignment(env)
aln.append_model(mdl, align_codes=code)

# Write the sequence to an alignment file
aln.write(file=f"{code}.seq")
