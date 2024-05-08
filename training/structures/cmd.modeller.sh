#########################################################################
# Author: Xingcheng Lin
# Created Time: Tue Dec 24 17:54:26 2019
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash


export Num_lines=`wc -l self_aligned_cdr3a | awk '{print $1}'`

for ((i=1; i<=$Num_lines; i++))
do
    echo $i
    rm -r modeller/job.$i
    mkdir -p modeller/job.$i
    cp template/* modeller/job.$i/
    TCR_name=`cat self_aligned_cdr3a | awk -v var="$i" '{if(NR==var)print $1}'`
    CDR3A_SEQ=`cat self_aligned_cdr3a | awk -v var="$i" '{if(NR==var)print $2}'`
    CDR3B_SEQ=`cat self_aligned_cdr3b | awk -v var="$i" '{if(NR==var)print $2}'`
    PEP_SEQ=`cat self_aligned_peptide | awk -v var="$i" '{if(NR==var)print $2}'`
    cd modeller/job.$i/
    sed "s/CDR3ASEQ/$CDR3A_SEQ/g; s/CDR3BSEQ/$CDR3B_SEQ/g; s/PEPTIDESEQ/$PEP_SEQ/g; s/TCRNAME/$TCR_name/g" template_alignment.ali > alignment.ali
    # Remove the carriage return
    sed -i 's/^M//g' alignment.ali
    sed "s/TCRNAME/$TCR_name/g" template_fillres.py > fillres.py

    python fillres.py

    # cp $TCR_name.BL00010001.pdb test.$i.pdb

    cd ../../
done
