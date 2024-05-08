#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 17 21:56:01 2020
# File Name: cmd.copyFile.sh
# Description: Copy the structure, phis and tms files from the 1ao7_template built each folders into the corresponding folders in this directory 
#########################################################################
#!/bin/bash

while read f
do
    echo $f
    cp ../$f/native_structures_pdbs_with_virtual_cbs/native.pdb native_structures_pdbs_with_virtual_cbs/${f}.pdb
    cp ../$f/phis/phi_pairwise_contact_well_native_native_-8.0_8.0_0.7_10 phis/phi_pairwise_contact_well_${f}_native_-8.0_8.0_0.7_10
    cp ../$f/phis/phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-8.0_8.0_0.7_10 phis/phi_pairwise_contact_well_${f}_decoys_CPLEX_randomization_-8.0_8.0_0.7_10
    cp ../$f/tms/native.tm tms/${f}.tm
done < proteinList.txt
