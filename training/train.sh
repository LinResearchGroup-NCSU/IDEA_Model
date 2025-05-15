#!/bin/sh

# find . -type f -name '*.sh' -print0 |
#   xargs -0 gsed -i.bak 's/\bsed\b/gsed/g'


#sed -i 's/\r//' proteinList.txt
#sed -i -e '$a\' proteinList.txt

#set -x

printf "[%s] [INFO] We will start training the IDEA energy model. The whole training dataset are listed below:\n" "$(date '+%Y-%m-%d %H:%M:%S')"
cat proteinList.txt

find . -mindepth 2 -name 'phi1_list.txt' -exec cp ./phi1_list.txt {} \;

find . -type f -name "proteinList.txt" -not -path "./proteinList.txt" -exec cp ./proteinList.txt {} \;

cd optimization/

bash cmd.optimization.sh

