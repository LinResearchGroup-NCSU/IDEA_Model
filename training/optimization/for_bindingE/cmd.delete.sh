#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Jun 18 22:59:21 2020
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash

while read f
do
    echo $f
    rm -r $f
done < proteinList.txt
