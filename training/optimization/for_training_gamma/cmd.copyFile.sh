#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 17 21:56:01 2020
# File Name: cmd.copyFile.sh
# Description: Copy the structure, phis and tms files from each folder built in for_bindingE/ 
# into the corresponding folders in this directory for training
#########################################################################
#!/bin/bash

while read f
do
    echo $f
    cp ../for_bindingE/$f/native_structures_pdbs_with_virtual_cbs/native.pdb native_structures_pdbs_with_virtual_cbs/${f}.pdb
    cp ../for_bindingE/$f/phis/phi_pairwise_contact_well_native_native_-8.0_8.0_0.7_10 phis/phi_pairwise_contact_well_${f}_native_-8.0_8.0_0.7_10
    cp ../for_bindingE/$f/phis/phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-8.0_8.0_0.7_10 phis/phi_pairwise_contact_well_${f}_decoys_CPLEX_randomization_-8.0_8.0_0.7_10
    cp ../for_bindingE/$f/tms/native.tm tms/${f}.tm
done < native_trainSetFiles.txt
