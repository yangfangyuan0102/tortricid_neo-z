conda activate orthofinder

#put all protein file *.fas into a folder

sed -i "s/\./__/g" ./*.fas
orthofinder -f ./ -M msa -t 30 -T raxml-ng -oa

#extract single copy gene sequence
mkdir -p ./SCO
ls ./Single_Copy_Orthologue_Sequences | parallel cp ./MultipleSequenceAlignments/{} ./SCO/{} 


#dn/ds
echo "32" > cpu
ParaAT.pl -h leafroller.paired -n orth.cds -a orth.fas -p cpu -m muscle -f axt -g -k -o dnds_leafroller
ParaAT.pl -h outgroup.paired -n orth.cds -a orth.fas -p proc -m muscle -f axt -g -k -o dnds_outgroup
