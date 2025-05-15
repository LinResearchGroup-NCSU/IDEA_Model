#########################################################################
# Author: Xingcheng Lin
# Created Time: Sat Feb 13 19:56:33 2021
# File Name: cmd.do_optimization.sh
# Description: 
#########################################################################
#!/bin/bash

cd for_bindingE/

# Run to generate phi file
bash cmd.for_phi.sh

# Copy the generated phi into the loocv folder
cd loocv/
bash cmd.copyFile.sh
cd ../../

# Do the training;
bash cmd.training.sh

# Delete the target folder for saving space, only delete after training is complete
# bash cmd.delete.sh
cd ../
