#!/bin/sh

sed -i 's/\r//' proteinList.txt
sed -i -e '$a\' proteinList.txt


find . -mindepth 2 -name 'phi1_list.txt' -exec cp ./phi1_list.txt {} \;

find . -type f -name "proteinList.txt" -not -path "./proteinList.txt" -exec cp ./proteinList.txt {} \;

cd optimization/

bash cmd.do_optimization.sh

