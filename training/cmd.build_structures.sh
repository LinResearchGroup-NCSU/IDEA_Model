#########################################################################
# Author: Xingcheng Lin
# Created Time: Sat Feb 13 19:56:33 2021
# File Name: cmd.create_template.sh
# Description: 
#########################################################################
#!/bin/bash

# We shall start with template_PDB, then gradually adding in new templates
template_PDB=$1

num_lines=$(wc -l template_list.txt | awk '{print $1}')
current_line_number=0

# Remove the existing template_alignment.ali
rm structures/template/template_alignment.ali

while read f
do
    current_line_number=$((current_line_number+1))
    echo $f
    
    # Copy the template structure 
    cp PDBs/$f.trunc.fitRenumberResidues.pdb structures/template/$f.pdb
    
    cd structures/template/


    # Get the template sequence
    python ~/Dropbox_work/scripts/python3/buildseq.py $f
    # Clean the sequence
    python cleanSeq.py $f.seq $f\_clean.seq
    # Create the template_alignment.ali
    python copy_seq_fortemplate.py $f\_clean.seq template_alignment.ali $f

    # Create the template_fillres.py, adding in the templates
    if [[ "$current_line_number" < "$num_lines" ]]
    then
        SRC="TEMPLATE_PDBID"
        DST="'$f',TEMPLATE_PDBID"
        gsed "s/$SRC/$DST/g" template_template_fillres.py > tmp.py
        mv tmp.py template_template_fillres.py 
    else
        SRC="TEMPLATE_PDBID"
        DST="'$f'"
        gsed "s/$SRC/$DST/g" template_template_fillres.py > template_fillres.py
    fi

    cd ../../

done < template_list.txt

cd structures/template/

# Create the target sequence, using the framework of 1ao7
python copy_seq_fortarget.py $template_PDB\_clean.seq template_alignment.ali $template_PDB
cd ../

# Do the Modeller
bash cmd.modeller.sh
    
# Copy file into the result/ folder
cd modeller/
bash cmd.copyfiles.sh
cd ../../
