#########################################################################
# Author: Xingcheng Lin
# Created Time: Sun Jun 14 16:58:44 2020
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash

# bash cmd.preprocessing.sh 3c6l A B 595 606

echo "Check the TCR chain IDs of the testBinder file (because Modeller changed it when rebuilding...) !!!"
echo "Type the eigen-index of PC1, followed by [ENTER]:"
read -rsp $'Press any key to continue...\n' -n1 key

bash cmd.evaluate_bindingE.sh 3c6l A B
