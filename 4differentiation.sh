
python parseVCF.py --skipIndels -i gmgd.filted.vcf.gz -o gmgd.filted.vcf.gz.geno.gz
python popgenWindows.py -w 100000 -s 100000 -g gmgd.M4.dip.vcf.gz.geno.gz --windType coordinate -o output.100k.csv -f phased -T 32  --popsFile popmap -p gm -p gd --writeFailedWindows

awk -F',' '{print $1"\t"$2"\t"$3"\t"$8}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.dxy.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$9}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.fst.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$6}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.pi_gm.bedgraph
awk -F',' '{print $1"\t"$2"\t"$3"\t"$7}' output.100k.csv | tail -n +2 | grep -v nan > output.100k.pi_gd.bedgraph

######find divergence peak regions
raw_fst <- read.table("output.100k.fst.bedgraph",header = F)
raw_da <- read.table("output.100k.da.bedgraph",header = F)
raw <- cbind(raw_fst,raw_da[,4])
colnames(raw) <- c("chr","start","end","fst","da")
raw$chr <- as.factor(raw$chr)

pdf("fst_da.pdf",width = 196,height = 14)
par(mfrow=c(2,28),mar=c(1,0,1,0))
chr_outliers <- matrix(nrow = 0,ncol = 0)
for (chr in unique(raw[,1])) {
  chr="gm_1"
  chr_raw <- raw[raw[,1]== chr,]
  x <- 1:nrow(chr_raw)
  loess_fst <- loess(chr_raw[,4] ~ x, span = 0.6) %>% predict()
  loess_da  <- loess(chr_raw[,5] ~ x, span = 0.6) %>% predict()
  sd_fst <- chr_raw[,4]-loess_fst
  sd_da <-  chr_raw[,5]-loess_da
  upper_fst <- 1*sd(chr_raw[,4]) + loess_fst
  upper_da  <- 1*sd(chr_raw[,5]) + loess_fst
  plot(x,chr_raw[,4],pch=20, xlim= c(0,350),ylim=c(0.1,0.7),xlab = "",ylab = "")
     lines(x,loess_fst)
     lines(x,upper_fst,col="red")
     abline(v = min(chr_raw[,5]), col = "blue", lty = 2) 
     abline(v = length(chr_raw[,5]), col = "green", lty = 2)
chr_outliers <- rbind(chr_outliers, chr_raw[chr_raw[,4]> upper_fst & chr_raw[,5]> upper_da,])
}
chr_outliers <- matrix(nrow = 0,ncol = 0)
for (chr in unique(raw[,1])) {
  chr_raw <- raw[raw[,1]== chr,]
  x <- 1:nrow(chr_raw)
  loess_fst <- loess(chr_raw[,4] ~ x, span = 0.6) %>% predict()
  loess_da  <- loess(chr_raw[,5] ~ x, span = 0.6) %>% predict()
  sd_fst <- chr_raw[,4]-loess_fst
  sd_da <-  chr_raw[,5]-loess_da
  upper_fst <- median(sd_fst) +  1*sd(sd_fst) + loess_fst
  upper_da  <- median(sd_da)  +  1*sd(sd_da)  + loess_da
  plot(x,chr_raw[,5],pch=20, xlim= c(0,350),ylim=c(0.02,0.09),xlab = "",ylab = "")
   lines(x,loess_da)
   lines(x,upper_da,col="red")
   abline(v =  min(chr_raw[,5]), col = "blue", lty = 2)
   abline(v = length(chr_raw[,5]), col = "green", lty = 2)
  chr_outliers <- rbind(chr_outliers, chr_raw[chr_raw[,4]> upper_fst & chr_raw[,5]> upper_da,])
}
dev.off()
