conda activate jcvi
cd /mnt/dataset/moth/synteny/MCSCAN
###################gff2mscan
for sp in `cat ../specieslist`
do
species=$sp
genome=/mnt/dataset/moth/species/${species}/2genome/${species}.fasta
gff=/mnt/dataset/moth/species/${species}/3anno/evm/${species}.evm.all.gff3
python3 -m jcvi.formats.gff bed --type=mRNA --key=ID ${gff} -o ${species}.bed
python3 -m jcvi.formats.bed uniq ${species}.bed
bedtools getfasta -fi ${genome} -fo ${species}.cds -bed ${species}.uniq.bed -s -name
sed -i "s/::.*//g" ${species}.cds
done
################mapping1by1
for s1 in `cat ../specieslist`;do
 s2=`grep -A 1 ${s1} ../specieslist | tail -1`
 if [ ${s1} != ${s2} ];then
    python3 -m jcvi.compara.catalog ortholog --no_strip_names --cpus=32 ${s1} ${s2}
    python3 -m jcvi.compara.synteny screen --minspan=50 --simple ${s1}.${s2}.lifted.anchors ${s1}.${s2}.lifted.anchors.new
 fi
done


#plot
#colorhttps://matplotlib.org/3.5.1/tutorials/colors/colors.html
#######################prepare layout file
cd /mnt/dataset/moth/synteny/MCSCAN
echo -e "# y, xstart, xend, rotation, color, label, va, bed, label_va" > ../layout
n=.05
for sp in `cat ../specieslist`;do
echo "${n}, 0.1, 0.8, 0, m, ${sp}, bottom, ${sp}.uniq.bed, center" >> ../layout
n=`echo "$n + 0.05" | bc`
done
echo -e "# edges" >> ../layout
n=0
for s1 in `cat ../specieslist`;do
s2=`grep -A 1 ${s1} ../specieslist | tail -1`
if [ ${s1} != ${s2} ];then
echo "e, ${n}, `echo "$n + 1" | bc`, ${s1}.${s2}.lifted.anchors.simple" >> ../layout
fi
n=`echo "$n + 1" | bc`
done
######################prepare seqid
for species in `cat ../specieslist`;do
genome=/mnt/dataset/moth/species/${species}/2genome/${species}.fasta.chr
echo `seqkit seq -n ${genome} | awk '{print $1}'`| tr " " "," >> ../seqid
done
######################plot
python3 -m jcvi.graphics.karyotype ../seqid ../layout --basepair  --keep-chrlabels --nocircles

