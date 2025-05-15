#!/bin/bash

while read f; do
    echo "Processing $f"
    cp ../../PDBs/${f}_modified.pdb native.pdb

    # Get the chain ID of the protein
    python find_prot_chainID.py native.pdb chain_ID_protein.txt || { echo "Failed to get chain ID for $f"; continue; }
    prot_chainID=$(<chain_ID_protein.txt)

    # Prepare working directory
    [ -d "$f" ] && rm -r "$f"
    cp -r template "$f"

    cd "$f" || { echo "Failed to enter folder $f"; continue; }

    # Directly call the preprocessing script with parameters
    bash cmd.preprocessing.sh "$f" "$prot_chainID" || echo "Failed to run preprocessing for $f"

    cd ../
done < proteinList.txt
