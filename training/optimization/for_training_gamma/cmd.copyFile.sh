#!/bin/bash

params=$(grep -v "^#" phi1_list.txt | grep "phi_pairwise_contact_well" | awk '{print $2"_"$3"_"_$4"_"$5}')

for dir in native_structures_pdbs_with_virtual_cbs phis tms; do
  [ -d "$dir" ] && find "$dir" -type f -delete
done

while read f
do
    echo "Copying $f's files"

    cp ../for_bindingE/$f/native_structures_pdbs_with_virtual_cbs/native.pdb native_structures_pdbs_with_virtual_cbs/${f}.pdb

    cp ../for_bindingE/$f/phis/phi_pairwise_contact_well_native_native_${params} phis/phi_pairwise_contact_well_${f}_native_${params}

    cp ../for_bindingE/$f/phis/phi_pairwise_contact_well_native_decoys_CPLEX_randomization_${params} phis/phi_pairwise_contact_well_${f}_decoys_CPLEX_randomization_${params}
    
    cp ../for_bindingE/$f/tms/${f}.tm tms/${f}.tm
    
done < native_trainSetFiles.txt




