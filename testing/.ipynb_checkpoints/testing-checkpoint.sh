#!/bin/bash

export PDBid=$1

protChain="A"


find . -type f -name "proteinList.txt" -not -path "./proteinList.txt" -exec cp ./proteinList.txt {} \;


cp PDBs/${PDBid}_modified.pdb native_structures_pdbs_with_virtual_cbs/native.pdb

#python add_fakeCB.py

#generate sequence
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native
bash cmd.cleanSequences.sh native.seq 

#dna_half.seq should be generated differently each time

# Use one half to generate the complementary other half
python reverse_complement.py dna_half.seq dna_half_complement.seq

#merge and generate whole sequence 
#paste -d '' dna_half.seq dna_half_complement.seq > dna.seq
python merge.py
#generate modeller whole sequence
python mapDNAseq_reverse.py dna.seq dna_modeller.seq

#combine modeller and original protein sequence 
python combine_DNAPro.py
# Cutoff for determining contacting residues, unit: nm
export cutoff=1.2
python find_cm_residues.py native.pdb $cutoff randomize_position_prot.txt randomize_position_DNA.txt

rm -r DNA_randomization
mkdir -p DNA_randomization 

cp randomize_position_DNA.txt native.seq native.decoys DNA_randomization/

# Generate decoys for the protein
#rm -r prot_randomization
#mkdir -p prot_randomization 

rm -r CPLEX_randomization
mkdir -p CPLEX_randomization

cat DNA_randomization/native.decoys prot_randomization/native.decoys > CPLEX_randomization/native.decoys


cd ../

# Create the tms file, where the DNA is labeled as '2', while the protein are labeled as '1';
grep "CA\|O5'" native_structures_pdbs_with_virtual_cbs/native.pdb > tmp.txt

# Get the total number of residues;
tot_resnum=`cat tmp.txt | awk 'END{print $6}'`
python create_tms.py sequences/DNA_randomization/randomize_position_DNA.txt $tot_resnum

sed "s/CPLEX_NAME/$PDBid/g; s/PROT_CHAIN/$protChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py

