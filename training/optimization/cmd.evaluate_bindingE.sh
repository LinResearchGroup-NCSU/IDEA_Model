#########################################################################
# Author: Xingcheng Lin
# Created Time: Sat Feb 13 19:56:33 2021
# File Name: cmd.TCRtraining.sh
# Description: Training the TCR with each protein in the training set,
# for separate gamma
#########################################################################
#!/bin/bash

# Create evaluateE_use_same_gammas.py
cp evaluateE_use_same_gammas.py for_bindingE/loocv/evaluateE_use_same_gammas.py

# Copy the trained gammma file in
cp for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_gamma_filtered for_bindingE/loocv/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_gamma_filtered

cd for_bindingE/loocv/

# Evaluate the binding energies for native and decoy binders
python evaluateE_use_same_gammas.py

mkdir -p results/
mv enative.*.txt emg.*.txt zscore.*.txt results/

cd ../../

