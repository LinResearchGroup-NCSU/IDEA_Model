NoLines=`cat native_sequences.txt | wc -l`

for ((i=1; i<=$NoLines; i++))
do
   # Get the sequence for that specific line;
   pSequence=$(gsed "${i}q;d" native_sequences.txt)
   echo $pSequence
   # Replace that sequence into the alignment file;
   gsed "s/TOBEREPLACED/$pSequence/g" native_template.seq > gBinder.$i.seq

done
# Combine the sequences together
cat $(ls gBinder.*.seq  | gsort -V) > native_gb.decoys
rm gBinder.*.seq
