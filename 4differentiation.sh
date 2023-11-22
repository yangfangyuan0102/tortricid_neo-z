
python parseVCF.py --skipIndels -i gmgd.filted.vcf.gz -o gmgd.filted.vcf.gz.geno.gz
python popgenWindows.py -w 100000 -s 100000 -g gmgd.M4.dip.vcf.gz.geno.gz --windType coordinate -o output.100k.csv -f phased -T 32  --popsFile popmap -p gm -p gd --writeFailedWindows

awk -F',' '{print $1"\t"$2"\t"$3"\t"$8}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.dxy.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$9}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.fst.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$6}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.pi_gm.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$7}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.pi_gd.bedgraph

##da
da=dxy-(pi_gm+pi_gd)/2

##resudial of fst against recombination rate
#R script
data <- read.table("/mnt/dataset/moth/species/gmgdcp/fstdxy/output.100k.fst.recombinationRate.bedgraph",header = T,sep = "\t")

fst_reb <- lm(fst ~ reb, data = data)
summary(fst_reb)
fst_reb$residuals
data$fst
data$residu <- fst_reb$residuals

write.table(data,"output.100k.da.fste.bedgraph",sep = "\t",quote = F,row.names = F)
