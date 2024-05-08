#########################################################################
# Author: Xingcheng Lin
# Created Time: Sat Feb 13 19:56:33 2021
# File Name: cmd.TCRtraining.sh
# Description: Training the TCR with each protein in the training set,
# for separate gamma
#########################################################################
#!/bin/bash

# Copy the list of training proteins in;
cp ./proteinList.txt for_training_gamma/native_trainSetFiles.txt

# Training the gamma
cd for_training_gamma/

# Copy the phi and other files in for training; This is different from the original 60 multi-memory training, because we have new contact maps from new training data
bash cmd.copyFile.sh

# Start the training for gamma

python optimize_gamma.py
cd ../
