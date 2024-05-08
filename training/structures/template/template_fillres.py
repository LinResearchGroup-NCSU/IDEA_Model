from modeller import *
from modeller.automodel import *    # Load the automodel class
import sys;
import time;

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Read in HETATM records from template PDBs
#env.io.hetatm = True


a = loopmodel(env, 
        alnfile = 'alignment.ali',
        knowns = ('3qeq'), 
        sequence = 'TCRNAME')

a.starting_model= 1
a.ending_model  = 1

a.loop.starting_model = 1
a.loop.ending_model   = 1
a.loop.md_level       = refine.fast

a.make()
