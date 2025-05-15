#!/bin/bash

log() {
  printf "[%s] [INFO] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1"
}

cd for_bindingE/

# Generate phi file for each protein
log "Starting generation of native and decoy phi files for each protein..."
bash cmd.for_phi.sh
log "All native and decoy phi files have been successfully generated."

# # Copy the generated phi files into the loocv folder
# log "Running cmd.copyFile.sh to copy phi files into the loocv folder..."
# cd loocv/
# bash cmd.copyFile.sh
# log "File copying completed."
# cd ../../

cd ../

log "Copying phi files into the for_training_gamma folder..."
cp ./proteinList.txt for_training_gamma/native_trainSetFiles.txt
cd for_training_gamma/
bash cmd.copyFile.sh
log "File copying completed."

# Do the training
log "Starting gamma file generation..."
python optimize_gamma.py
log "Training has been completed."
