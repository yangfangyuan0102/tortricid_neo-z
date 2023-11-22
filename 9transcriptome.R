#run in shell
conda activate tpm
#gff32gtf
conda install agat
agat_convert_sp_gff2gtf.pl --gff *.gff3 -o *.gtf
gtf=/mnt/dataset/moth/species/gm/3anno/evm/gm.evm.all.gtf
TPMCalculator -b ./*.bam -g  $gtf  -e -a -p

#R script
#bootstrap resample
library(boot)
library(dplyr)

trans_data <- read.table("/mnt/dataset/moth/species/gm/1trans/Adult_ZW/gm.trans0_genes.out",header=T)[,c(1:4,7)]
  all <- trans_data[,5]
  M20 <- trans_data[1:164,5]
  M17 <- trans_data[165:666,5]
F20_17<- trans_data[1:666,5]
  MZ  <- trans_data[667:1219,5]
 auto <- trans_data[1219:14513,5]
 NEOZ <- trans_data[trans_data$Chr=="gm_1" ,5] #chromosome names may be different
 chr2 <- trans_data[trans_data$Chr=="gm_2" ,5]
 chr3 <- trans_data[trans_data$Chr=="gm_3" ,5]
 chr4 <- trans_data[trans_data$Chr=="gm_4" ,5]
 chr5 <- trans_data[trans_data$Chr=="gm_5" ,5]
 chr6 <- trans_data[trans_data$Chr=="gm_6" ,5]
 chr7 <- trans_data[trans_data$Chr=="gm_7" ,5]
 chr8 <- trans_data[trans_data$Chr=="gm_8" ,5]
 chr9 <- trans_data[trans_data$Chr=="gm_9" ,5]
chr10 <- trans_data[trans_data$Chr=="gm_10",5]
chr11 <- trans_data[trans_data$Chr=="gm_11",5]
chr12 <- trans_data[trans_data$Chr=="gm_12",5]
chr13 <- trans_data[trans_data$Chr=="gm_13",5]
chr14 <- trans_data[trans_data$Chr=="gm_14",5]
chr15 <- trans_data[trans_data$Chr=="gm_15",5]
chr16 <- trans_data[trans_data$Chr=="gm_16",5]
chr17 <- trans_data[trans_data$Chr=="gm_17",5]
chr18 <- trans_data[trans_data$Chr=="gm_18",5]
chr19 <- trans_data[trans_data$Chr=="gm_19",5]
chr20 <- trans_data[trans_data$Chr=="gm_20",5]
chr21 <- trans_data[trans_data$Chr=="gm_21",5]
chr22 <- trans_data[trans_data$Chr=="gm_22",5]
chr23 <- trans_data[trans_data$Chr=="gm_23",5]
chr24 <- trans_data[trans_data$Chr=="gm_24",5]
chr25 <- trans_data[trans_data$Chr=="gm_25",5]
chr26 <- trans_data[trans_data$Chr=="gm_26",5]
chr27 <- trans_data[trans_data$Chr=="gm_27",5]
chr28 <- trans_data[trans_data$Chr=="gm_28",5]

{ all <-   all[  all>0.1 ]
 NEOZ <-  chr1[ chr1>0.1 ]
  M20 <-   M20[  M20>0.1 ];    M20_boot <- replicate(1000, sample(  M20, size = 100, replace = T)) %>% apply(. , 2, median)
  M17 <-   M17[ M17>0.1 ];     M17_boot <- replicate(1000, sample(  M17, size = 100, replace = T)) %>% apply(. , 2, median)
F20_17 <-F20_17[F20_17>0.1];F20_17_boot <- replicate(1000, sample(F20_17,size = 100, replace = T)) %>% apply(. , 2, median)
   MZ <-    MZ[   MZ>0.1 ];    MZ_boot <- replicate(1000, sample(   MZ, size = 100, replace = T)) %>% apply(. , 2, median)
 auto <-  auto[ auto>0.1 ];  auto_boot <- replicate(1000, sample( auto, size = 100, replace = T)) %>% apply(. , 2, median)
 chr2 <-  chr2[ chr2>0.1 ];  chr2_boot <- replicate(1000, sample( chr2, size = 100, replace = T)) %>% apply(. , 2, median)
 chr3 <-  chr3[ chr3>0.1 ];  chr3_boot <- replicate(1000, sample( chr3, size = 100, replace = T)) %>% apply(. , 2, median)
 chr4 <-  chr4[ chr4>0.1 ];  chr4_boot <- replicate(1000, sample( chr4, size = 100, replace = T)) %>% apply(. , 2, median)
 chr5 <-  chr5[ chr5>0.1 ];  chr5_boot <- replicate(1000, sample( chr5, size = 100, replace = T)) %>% apply(. , 2, median)
 chr6 <-  chr6[ chr6>0.1 ];  chr6_boot <- replicate(1000, sample( chr6, size = 100, replace = T)) %>% apply(. , 2, median)
 chr7 <-  chr7[ chr7>0.1 ];  chr7_boot <- replicate(1000, sample( chr7, size = 100, replace = T)) %>% apply(. , 2, median)
 chr8 <-  chr8[ chr8>0.1 ];  chr8_boot <- replicate(1000, sample( chr8, size = 100, replace = T)) %>% apply(. , 2, median)
 chr9 <-  chr9[ chr9>0.1 ];  chr9_boot <- replicate(1000, sample( chr9, size = 100, replace = T)) %>% apply(. , 2, median)
chr10 <- chr10[chr10>0.1 ]; chr10_boot <- replicate(1000, sample(chr10, size = 100, replace = T)) %>% apply(. , 2, median)
chr11 <- chr11[chr11>0.1 ]; chr11_boot <- replicate(1000, sample(chr11, size = 100, replace = T)) %>% apply(. , 2, median)
chr12 <- chr12[chr12>0.1 ]; chr12_boot <- replicate(1000, sample(chr12, size = 100, replace = T)) %>% apply(. , 2, median)
chr13 <- chr13[chr13>0.1 ]; chr13_boot <- replicate(1000, sample(chr13, size = 100, replace = T)) %>% apply(. , 2, median)
chr14 <- chr14[chr14>0.1 ]; chr14_boot <- replicate(1000, sample(chr14, size = 100, replace = T)) %>% apply(. , 2, median)
chr15 <- chr15[chr15>0.1 ]; chr15_boot <- replicate(1000, sample(chr15, size = 100, replace = T)) %>% apply(. , 2, median)
chr16 <- chr16[chr16>0.1 ]; chr16_boot <- replicate(1000, sample(chr16, size = 100, replace = T)) %>% apply(. , 2, median)
chr17 <- chr17[chr17>0.1 ]; chr17_boot <- replicate(1000, sample(chr17, size = 100, replace = T)) %>% apply(. , 2, median)
chr18 <- chr18[chr18>0.1 ]; chr18_boot <- replicate(1000, sample(chr18, size = 100, replace = T)) %>% apply(. , 2, median)
chr19 <- chr19[chr19>0.1 ]; chr19_boot <- replicate(1000, sample(chr19, size = 100, replace = T)) %>% apply(. , 2, median)
chr20 <- chr20[chr20>0.1 ]; chr20_boot <- replicate(1000, sample(chr20, size = 100, replace = T)) %>% apply(. , 2, median)
chr21 <- chr21[chr21>0.1 ]; chr21_boot <- replicate(1000, sample(chr21, size = 100, replace = T)) %>% apply(. , 2, median)
chr22 <- chr22[chr22>0.1 ]; chr22_boot <- replicate(1000, sample(chr22, size = 100, replace = T)) %>% apply(. , 2, median)
chr23 <- chr23[chr23>0.1 ]; chr23_boot <- replicate(1000, sample(chr23, size = 100, replace = T)) %>% apply(. , 2, median)
chr24 <- chr24[chr24>0.1 ]; chr24_boot <- replicate(1000, sample(chr24, size = 100, replace = T)) %>% apply(. , 2, median)
chr25 <- chr25[chr25>0.1 ]; chr25_boot <- replicate(1000, sample(chr25, size = 100, replace = T)) %>% apply(. , 2, median)
chr26 <- chr26[chr26>0.1 ]; chr26_boot <- replicate(1000, sample(chr26, size = 100, replace = T)) %>% apply(. , 2, median)
chr27 <- chr27[chr27>0.1 ]; chr27_boot <- replicate(1000, sample(chr27, size = 100, replace = T)) %>% apply(. , 2, median)
chr28 <- chr28[chr28>0.1 ]; chr28_boot <- replicate(1000, sample(chr28, size = 100, replace = T)) %>% apply(. , 2, median)
}


pdf("trans.0.1.boot.zwxinxian.pdf",width = 5,height = 5)
par(mfrow=c(32,1),mar=c(0.1,0.1,0.1,0.1))
density(  M20_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(  M20_boot),col="#FFA500");abline(v=median(  M20),col="#BA55D3")
density(  M17_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(  M17_boot),col="#FFA500");abline(v=median(  M17),col="#BA55D3")
density(F20_17_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(F20_17_boot),col="#FFA500");abline(v=median(F20_17),col="#BA55D3")
density(   MZ_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(   MZ_boot),col="#FFA500");abline(v=median(  MZ ),col="#BA55D3")
density( auto_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( auto_boot),col="#FFA500");abline(v=median( auto),col="#BA55D3")
density( chr2_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr2_boot),col="#FFA500");abline(v=median( chr2),col="#BA55D3")
density( chr3_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr3_boot),col="#FFA500");abline(v=median( chr3),col="#BA55D3")
density( chr4_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr4_boot),col="#FFA500");abline(v=median( chr4),col="#BA55D3")
density( chr5_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr5_boot),col="#FFA500");abline(v=median( chr5),col="#BA55D3")
density( chr6_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr6_boot),col="#FFA500");abline(v=median( chr6),col="#BA55D3")
density( chr7_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr7_boot),col="#FFA500");abline(v=median( chr7),col="#BA55D3")
density( chr8_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr8_boot),col="#FFA500");abline(v=median( chr8),col="#BA55D3")
density( chr9_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean( chr9_boot),col="#FFA500");abline(v=median( chr9),col="#BA55D3")
density(chr10_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr10_boot),col="#FFA500");abline(v=median(chr10),col="#BA55D3")
density(chr11_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr11_boot),col="#FFA500");abline(v=median(chr11),col="#BA55D3")
density(chr12_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr12_boot),col="#FFA500");abline(v=median(chr12),col="#BA55D3")
density(chr13_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr13_boot),col="#FFA500");abline(v=median(chr13),col="#BA55D3")
density(chr14_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr14_boot),col="#FFA500");abline(v=median(chr14),col="#BA55D3")
density(chr15_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr15_boot),col="#FFA500");abline(v=median(chr15),col="#BA55D3")
density(chr16_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr16_boot),col="#FFA500");abline(v=median(chr16),col="#BA55D3")
density(chr17_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr17_boot),col="#FFA500");abline(v=median(chr17),col="#BA55D3")
density(chr18_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr18_boot),col="#FFA500");abline(v=median(chr18),col="#BA55D3")
density(chr19_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr19_boot),col="#FFA500");abline(v=median(chr19),col="#BA55D3")
density(chr20_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr20_boot),col="#FFA500");abline(v=median(chr20),col="#BA55D3")
density(chr21_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr21_boot),col="#FFA500");abline(v=median(chr21),col="#BA55D3")
density(chr22_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr22_boot),col="#FFA500");abline(v=median(chr22),col="#BA55D3")
density(chr23_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr23_boot),col="#FFA500");abline(v=median(chr23),col="#BA55D3")
density(chr24_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr24_boot),col="#FFA500");abline(v=median(chr24),col="#BA55D3")
density(chr25_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr25_boot),col="#FFA500");abline(v=median(chr25),col="#BA55D3")
density(chr26_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr26_boot),col="#FFA500");abline(v=median(chr26),col="#BA55D3")
density(chr27_boot) %>% plot(main="",xlab="",ylab="",xaxt="n",bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr27_boot),col="#FFA500");abline(v=median(chr27),col="#BA55D3")
density(chr28_boot) %>% plot(main="",xlab="",ylab="",         bty="n",xlim = c(0,15),ylim=c(0,0.5),col="gray");abline(v=mean(chr28_boot),col="#FFA500");abline(v=median(chr28),col="#BA55D3")
dev.off()


###########cutoff robustness
#5-0.1
ratio <- matrix(nrow=50,ncol=4)
row=1
for (cut in seq(0.1,5,by=0.1)) {
  M20_cut <- M20[M20 > cut]
  M17_cut <- M17[M17 > cut]
  MZ_cut  <- MZ[MZ > cut]
  auto_cut <- auto[auto > cut]
  ratio[row,1] <- cut
  ratio[row,2] <-  median(M20_cut)median(auto_cut)
  ratio[row,3] <-  median(M17_cut)/median(auto_cut)
  ratio[row,4] <-  median(MZ_cut)/median(auto_cut)
  row=row+1
}

pdf("a_z_TPMrate_ZZ.pdf",width = 8,height = 6)
num_intervals <- ceiling((max(all) - min(all)) / 0.1)
hist(all, xlim = c(0,5), breaks = seq(min(all), min(all) + num_intervals*0.1, 0.1),xlab = "TPM", ylab = "Frequency")
plot (ratio[,c(1,2)],type = "line",ylim = c(0,1.5),col="#184492")
lines(ratio[,c(1,3)],type = "line",ylim = c(0,1.5),col="#C3C2E0",add=T)
lines(ratio[,c(1,4)],type = "line",ylim = c(0,1.5),col="#C13328",add=T)
dev.off()

#######mutiple linear model
trans_data1 <- read.table("/mnt/dataset/moth/species/gm/1trans/Adult_ZW/gm.trans0_genes.out",header=T)[,c(1:4,7)]
trans_data2 <- read.table("/mnt/dataset/moth/species/gm/1trans/Adult_ZZ/gm.trans0_genes.out",header=T)[,c(1:4,7)]
gene_trans1 <- trans_data1[,2:5]
gene_trans1 <- gene_trans1[order(gene_trans1[,1]),]
gene_trans2 <- trans_data2[,2:5]
gene_trans2 <- gene_trans2[order(gene_trans2[,1]),]


#calculated pi_n and pi_s for each gene

#load
genepi <- read.csv("gene.pi.csv")[,c(1,2,3,6,7)]
genepi <- genepi[order(genepi[,1]),]
genepi[,1]

piRna <- cbind(genepi,(gene_trans1[,4]+gene_trans2[,4])/2)
colnames(piRna) <- c("chr","start","end","pin","pis","tpm")

piRnalog <- piRna
piRnalog$tpm <- log(piRna$tpm)
piRnalog <- na.omit(piRnalog)
piRnalog <- piRnalog[!is.infinite(piRnalog$tpm),]
piRnalog$type <- c(rep("z",times=sum(piRnalog$chr=="gm_1")),rep("a",times=sum(piRnalog$chr!="gm_1")))

#the lm model
lmpiRna <- lm(pin ~ log(tpm) + pis + type, piRnalog)

summary(lmpiRna)

pdf("model_dig.pdf",width = 10,height = 10)
par(mfrow=c(2,2))
plot(lmpiRna,col="green")
dev.off()





















