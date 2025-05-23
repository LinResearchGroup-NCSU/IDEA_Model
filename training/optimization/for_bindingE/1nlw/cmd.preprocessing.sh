#!/bin/bash

export PDBid=$1
export protChain=$2

cp ../../../PDBs/${PDBid}_modified.pdb native_structures_pdbs_with_virtual_cbs/native.pdb

# Generate randomized sequence for the decoys;
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native

cat native.seq > gBinder_sequences.txt

# Find the indices of contacting protein-DNA residues
# Cutoff for determining contacting residues, unit: nm
# echo "Cutoff for determining contacting residues"
export cutoff=1.2
python find_cm_residues.py native.pdb $cutoff randomize_position_prot.txt randomize_position_DNA.txt

# Generate decoys for the DNA
rm -rf DNA_randomization
mkdir -p DNA_randomization

cp randomize_position_DNA.txt native.seq gBinder_sequences.txt DNA_randomization/

python generate_decoy_seq_DNA.py

# decoys for the protein
rm -rf prot_randomization
mkdir -p prot_randomization

cp randomize_position_prot.txt native.seq gBinder_sequences.txt prot_randomization/

python generate_decoy_seq_prot.py

# Combine the generated DNA and protein decoys together
rm -rf CPLEX_randomization
mkdir -p CPLEX_randomization

cat DNA_randomization/native.decoys prot_randomization/native.decoys > CPLEX_randomization/native.decoys

cd ../

# Create the tms file, where the DNA is labeled as '2', while the protein are labeled as '1';
grep "CA\|O5'" native_structures_pdbs_with_virtual_cbs/native.pdb > tmp.txt

# Get the total number of residues;
tot_resnum=$(grep '^ATOM' tmp.txt | wc -l)
# echo $tot_resnum

python create_tms.py sequences/DNA_randomization/randomize_position_DNA.txt $tot_resnum $PDBid

sed "s/CPLEX_NAME/$PDBid/g; s/PROT_CHAIN/$protChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py