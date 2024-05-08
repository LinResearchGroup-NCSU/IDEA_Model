NoLines=`cat decoy_sequences.txt | wc -l`

for ((i=1; i<=$NoLines; i++))
do
   # Get the sequence for that specific line;
   pSequence=$(gsed "${i}q;d" decoy_sequences.txt)
   echo $pSequence
   # Replace that sequence into the alignment file;
   gsed "s/TOBEREPLACED/$pSequence/g" native_template.seq > wBinder.$i.seq

done
# Combine the sequences together
cat $(ls wBinder.*.seq  | gsort -V) > native_wb.decoys
rm wBinder.*.seq
