##################################################smc++
conda activate smc

vcf=/mnt/ssd/gm/3gvcf/gm.auto.snp.vcf.gz
genome=/mnt/dataset/moth/species/gm/2genome/gm.fasta.chr.fai

#create mask file
bedtools complement -i $vcf -g $genome | \
bedtools sort | bedtools merge| bgzip > ${vcf}.mask.bed.gz
tabix -b2 -e3 -f ${vcf}.mask.bed.gz

#vcf2smc
mask=${vcf}.mask.bed.gz

bcftools index -f $vcf
mkdir -p ${vcf}_smc
ind=`bcftools query -l $vcf | tr '\n' ','| sed 's/,$/\n/'`
tabix -l ${vcf} | parallel -j 20 "smc++ vcf2smc --cores 1 -m ${mask} ${vcf} ${vcf}_smc/{}.smc.gz {} gm:${ind}"


#run
smc++ estimate --em-iterations 20 --base gm_ -o ${vcf}_smc/ 2.9e-9 ${vcf}_smc/*.smc.gz --cores 32 -rp 4 --spline cubic  -w 200 --timepoints 4e2 4e7
smc++ plot ${vcf}_smc/gm_png -g 0.25 ${vcf}_smc/gm_.final.json --csv --linear -x 0 1e6


###############################################RNLeRNN
URTR="35"
GENOME="/mnt/ssd/gd/gm.fasta.chr.bed"
DIR="/mnt/ssd/gm/reb"
VCF="/mnt/ssd/gm/gm.dp5gq20miss50.tri.bi.snp.p.vcf"
SMC="/mnt/ssd/gm/gm.dp5gq20miss50.tri.vcf.gz_smc/gm_png.csv"
MASK="/mnt/ssd/gm/gm.dp5gq20miss50.tri.vcf.gz.mask.bed"

#VCF:biallelic, noly variation sites
#--forceDiploid
#--nEpochs 20 --nValSteps 2
ReLERNN_SIMULATE --vcf ${VCF} --genome ${GENOME}  --mask ${MASK} \
    --projectDir ${DIR} --assumedMu 2.9e-9 --upperRhoThetaRatio ${URTR} \
    --seed 42 --nCPU 20 --phased  --assumedGenTime 0.25 \
    --demographicHistory ${SMC} && \
ReLERNN_TRAIN --projectDir ${DIR} --seed 42 --gpuID 0  && \
ReLERNN_PREDICT --vcf ${VCF} --projectDir ${DIR} --seed 42 --unphased --gpuID 0 && \
ReLERNN_BSCORRECT --projectDir ${DIR}  --nCPU 30 --gpuID 0

#plot in R
reb <- read.table("./gm.snp.PREDICT.BSCORRECTED.txt",header = T)
reb$pos <- (reb$start+reb$end)/2

reb_chr <- split(reb,reb[,1])
chr.size <- read.table("/mnt/dataset/moth/species/gm/2genome/gm.fasta.chr.fai",header = F)[,1:2]
for (chr in 1:28) {
  name <- paste("gm_",chr,sep = "")
  gm <- ggplot(reb_chr[[name]]) +  theme_bw()+  xlim(0,chr.size[chr,2]) + ylim(0,8e-8) + ylab("") +xlab("")+
    scale_x_continuous(breaks = c(0, chr.size[chr,3]), labels = c(0, chr.size[chr,3])) +
    theme(axis.text.x = element_text(size = 4),axis.text.y = element_text(size = 4)) +
    geom_point(aes(x=pos,y=recombRate),color="grey",size=0.1) +
    geom_smooth(aes(x=pos,y=recombRate),color="#ef8a62", method = "gam") +
    geom_hline(yintercept = median(reb_chr[[name]][,5]),color="black",linetype="dashed")
  assign(name,gm)
}
pdf("gm_reb.CORRECT_gam.pdf",width = 50,height = 5)
grid.arrange(gm_1, gm_2,gm_3,gm_4,gm_5 ,gm_6,gm_7,gm_8,gm_9,gm_10,gm_11,gm_12,gm_13,
             gm_14,gm_15,gm_16,gm_17,gm_18,gm_19,gm_20,gm_21,gm_22,gm_23,gm_24,gm_25,gm_26,gm_27,gm_28,
             ncol=28,nrow=1)
dev.off()

#recombination rate median for each cheomosome
median(reb_chr[[name]][,5])

#chr length corrected with reb
for (chr in 1:28) {
  name <- paste("gm_",chr,sep = "")
chr.size[chr,3] <- median(reb_chr[[name]][,5])
}

pdf("/mnt/dataset/moth/Analysis/figure/person.pdf", width = 5,height = 3)
ggplot(chr.size[-1,],aes(x=V2,y=V3)) +
  geom_point() +
  scale_x_continuous(limits = c(min(chr.size[,2]),max(chr.size[,2]) ))+
  scale_y_continuous(limits = c(1.256349e-08,max(chr.size[,3])))+
  geom_smooth(method = "lm") +
  theme_bw()
dev.off()
chrlength_lm <- lm(V3 ~ V2,chr.size[-1,])
#presicte recombination rate of neo-Z length
predict(chrlength_lm,chr.size[1,]) 
#=1.256349e-08
#obversed=8.846138e-09

cor(chr.size[,3],chr.size[,2],method = "spearman")


###########################################randomForest
#R script
##random input
#reconbination rate dataset
reb <- read.table("./100kWindowsed.reb.bedgraph",header = F)
reb$pos <- (reb$V2+reb$V3)/2
reb_chr <- split(reb,reb[,1])
chr.size <- read.table("/mnt/ssd/gm/gm.fasta.fai",header = F)[,c(1,2)]
#chromosome length dataset
for (chr in 1:28) {
  chr_name <- paste("gm_",chr,sep = "")
  reb_chr[[paste("gm_",chr,sep = "")]]$distance <- chr.size[chr,2]-reb_chr[[paste("gm_",chr,sep = "")]][,5]
}
reb_chr <- do.call(rbind,reb_chr)
reb_chr$dis_tel <- NA
reb_chr$relative_dis_tel <- NA
for (i in 1:nrow(reb_chr)) {
  if (reb_chr$pos[i] < reb_chr$distance[i]) {
    reb_chr$dis_tel[i] <- reb_chr$pos[i]
  } else {
    reb_chr$dis_tel[i] <- reb_chr$distance[i]
  }
  reb_chr$relative_dis_tel[i] <- reb_chr$dis_tel[i]*2 / chr.size[chr.size[,1]==reb_chr[i,1],2]
}


#train dataset
reb_input <- reb_chr[reb_chr[,1] != "gm_1" & reb_chr[,1] != "gm_2" & reb_chr[,1] != "gm_3" & reb_chr[,1] != "gm_4" & reb_chr[,1] != "gm_17", ] #remove these chromosome
#predict dataset
reb_chr1 <- reb_chr[reb_chr[,1] == "gm_1",]
reb_chr2 <- reb_chr[reb_chr[,1] == "gm_2",]
reb_chr3 <- reb_chr[reb_chr[,1] == "gm_3",]
#model
library(randomForest)
model <- randomForest(V4 ~ dis_tel + relative_dis_tel, data=reb_input, ntree=1000)
importance(model)
pred_chr1 <- predict(model,newdata = reb_chr1)
pred_chr2 <- predict(model,newdata = reb_chr2)
pred_chr3 <- predict(model,newdata = reb_chr3)

pdf("pred_obs.cor.pdf")
plot(pred_chr1 ,reb_chr1[,4]*1.420223,xlim = c(0,5e-8),ylim = c(0,5e-8),pch=19,col="#bdbdbd")
abline(lm(reb_chr1[,4]*1.420223 ~ pred_chr1 ),col="blue")
abline(a=0,b=1,lty=2)
plot(pred_chr2 ,reb_chr2[,4],xlim = c(0,5e-8),ylim = c(0,5e-8),pch=19,col="#bdbdbd")
abline(lm(reb_chr2[,4] ~ pred_chr2 ),col="blue")
abline(a=0,b=1,lty=2)
plot(pred_chr3 ,reb_chr3[,4],xlim = c(0,5e-8),ylim = c(0,5e-8),pch=19,col="#bdbdbd")
abline(lm(reb_chr3[,4] ~ pred_chr3 ),col="blue")
abline(a=0,b=1,lty=2)
dev.off()

confint(lm(reb_chr1[,4]*1.420223 ~ pred_chr1 ))
confint(lm(reb_chr2[,4] ~ pred_chr2 ))
confint(lm(reb_chr3[,4] ~ pred_chr3 ))

wilcox.test(pred_chr1 ,reb_chr1[,4]**1.420223,paired = T)
wilcox.test(pred_chr2 ,reb_chr2[,4],    paired = T)
wilcox.test(pred_chr3 ,reb_chr3[,4],    paired = T)





