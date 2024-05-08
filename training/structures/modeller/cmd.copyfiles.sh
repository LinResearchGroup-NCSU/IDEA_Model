#########################################################################
# Author: Xingcheng Lin
# Created Time: Sat Feb 13 19:56:33 2021
# File Name: cmd.copyfiles.sh
# Description: 
#########################################################################
#!/bin/bash

acc=1
while read f
do
    echo $f $acc
    # Because job.1 is for the template, so we start from the job.2
    if [[ "$acc" > 1 ]]
    then
        # If Modeller can build the loop structure, use the optimized structure, otherwise, use the automodel built structure
        if [ -f "job.$acc/$f.BL00010001.pdb" ]
        then
            cp job.$acc/$f.BL00010001.pdb results/
        else
            cp job.$acc/$f.B99990001.pdb results/$f.BL00010001.pdb
        fi
    fi
    acc=$((acc+1))
done < HLA0201_list.txt
