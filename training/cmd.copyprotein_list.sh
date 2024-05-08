#########################################################################
# Author: Xingcheng Lin
# Created Time: Sun Dec  5 11:18:36 2021
# File Name: cmd.copyprotein_list.sh
# Description: proteinList.txt includes all the strong binder labels and the template PDB
#########################################################################
#!/bin/bash

tail -n +0 proteinList.txt >  ./optimization/proteinList.txt
tail -n +0 proteinList.txt > ./optimization/for_bindingE/proteinList.txt
tail -n +0 proteinList.txt > ./optimization/for_bindingE/loocv/proteinList.txt
tail -n +0 proteinList.txt > ./optimization/for_training_gamma/proteinList.txt
tail -n +0 proteinList.txt > ./optimization/for_training_gamma/phis/proteinList.txt

# In the structure file, we need to add the template PDB
#cp proteinList.txt structures/modeller/

# We copy the updated phi1 parameters into the protocol
cp phi1_list.txt ./optimization/for_training_gamma/
cp phi1_list.txt ./optimization/for_bindingE/template/
cp phi1_list.txt ./optimization/for_bindingE/loocv/
