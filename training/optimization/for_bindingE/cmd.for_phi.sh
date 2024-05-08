#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Jun 18 22:59:21 2020
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash

while read f
do
    echo $f 
    cp ../../PDBs/${f}_modified.pdb native.pdb

    # Get the chain ID of the protein;
    python find_prot_chainID.py native.pdb chain_ID_protein.txt
    prot_chainID=$(cat chain_ID_protein.txt)

    # Copy the template folder into the PDB folder
    rm -r $f
    cp -r template $f
    # update the cmd.optimization.sh file the residue ID of peptide and PDB id of the corresponding PDB folder;
    gsed "s/PDBID/$f/g; s/PROT_CHAIN_ID/$prot_chainID/g" template_cmd.optimization.sh > $f/cmd.optimization.sh

    cd $f/
    
    bash cmd.optimization.sh
    cd ../
done < proteinList.txt


