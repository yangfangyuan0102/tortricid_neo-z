conda activate hic

species=gm
genome=/mnt/dataset/moth/species/${species}/2genome/${species}.fasta
work_dir=/mnt/dataset/moth/species/${species}/4hic
r1=/mnt/dataset/moth/species/${species}/4hic/0raw/r1.fq.gz
r2=/mnt/dataset/moth/species/${species}/4hic/0raw/r2.fq.gz

cd $work_dir
bwa-mem2 index $genome
bwa-mem2 mem -A1 -B4 -E50 -L0 -t 32 $genome  $r1 | samtools view -Shb - > ./1bam/mate_R1.bam
bwa-mem2 mem -A1 -B4 -E50 -L0 -t 32 $genome  $r2 | samtools view -Shb - > ./1bam/mate_R2.bam

###############Hi-C matrix
hicFindRestSite --fasta $genome --searchPattern GATC -o ${genome}_hic.bed
hicBuildMatrix --samFiles ./1bam/mate_R1.bam ./1bam/mate_R2.bam --binSize 50000 \
                 --restrictionSequence GATC --danglingSequence GATC \
                 --restrictionCutFile ${genome}_hic.bed \
                 --threads 30 --inputBufferSize 20000 -o ./2matrix/50k/50k.h5 --QCfolder ./hicQC \
                 --outBam ./1bam/hic.bam
################A/B compartment
#cscoretools
bedtools bamtobed -bedpe -i ./1bam/hic.bam | awk '{print $7"\t"$1"\t"$3"\t"$9"\t"$4"\t"$6"\t"$10}' > ./1bam/hic.summary
generateEqualLengthBed ${genome}.fai ./2matrix/50k.bed 50000
cscoretool1.1 ./2matrix/50k.bed ./1bam/hic.summary ./2matrix/50k/50k.cscore 10 1000000

#R script
##adjust A/B compartment
library(bedtoolsr)
species='gm'
gff=paste('/mnt/dataset/moth/species/',species,'/3anno/evm/',species,'.evm.all.gff3',sep = '')
genome=paste('/mnt/dataset/moth/species/',species,'/2genome/',species,'.fasta.chr.fai',sep = '')
bed <- read.table(paste('/mnt/dataset/moth/species/',species,'/4hic/2matrix/50k/50k.cscore_cscore.bedgraph',sep = ''))

chr <- split(bed,bed[,1])
for (c in 1:length(chr)) {
  up   <- bt.coverage(a=chr[[c]][chr[[c]]$V4>0,],b=gff)
  down <- bt.coverage(a=chr[[c]][chr[[c]]$V4<0,],b=gff)
  if (mean(up$V8)< mean(down$V8) ) {
    chr[[c]][,4]=chr[[c]][,4]*-1 }
}
write.table(unsplit(chr,bed[,1]),paste('/mnt/dataset/moth/species/',species,'/4hic/2matrix/50k/50k.cscore_cscore.adj.bedgraph',sep = ''),
            quote = F,col.names = F,row.names = F,sep = "\t" )
#end

#############average Cscore within Tortricidae and within outfroups
library(dplyr)
library(stringr)
library(bedtoolsr)
#List of homologous genes
OG_genename <- read.table("/mnt/dataset/moth/Analysis/5species/orthfinder/OrthoFinder/Results_Mar31/Orthogroups/Orthogroups.tsv",
                          header = T,row.names = 1,sep = "\t")
#load the gene-species matrix
ABscore <- read.table("/mnt/dataset/moth/Analysis/5species/orthfinder/OrthoFinder/Results_Mar31/Orthogroups/Orthogroups.GeneCount.tsv",
                      header = T,row.names = 1,sep = "\t")[,-c(5,6)] 
                      #modify [,-c(5,6)] to remove species columns such as outgroups. 
                      #Modify again to calculate average Cscores within outgroups
ABscore[,] <- NA
genecount=nrow(OG_genename) #gene count
sp_list <- colnames(ABscore) #species_list
gene_list <- row.names(ABscore) #gene_list
#get average Cscore for each genes
for (sp in sp_list) { #loop in species
  cscore.bed=paste("/mnt/dataset/moth/species/",sp,"/4hic/2matrix/50k/50k.cscore_cscore.sort.bedgraph",sep = "") #c.score
  gene.gff="/mnt/dataset/moth/species/gm/3anno/evm/gm.evm.all.gff3.mRNA.sort.gff3" #
  gene.cscore <- bt.map(a=gene.gff,b=cscore.bed,c=4,o="mean")[,9:10]
  OG_genename[,sp] <- OG_genename[,sp] %>% gsub("__",".",.)
  for (OG in gene_list) { #loop in homologous genes
    OG_each_gene <- OG_genename[OG,sp] %>% strsplit(.,", ") %>% unlist
    OG_each_gene_score <- apply(matrix(OG_each_gene),1,function(x) gene.cscore[grep(x,gene.cscore[,1]),2] %>% as.numeric() )
    ABscore[OG,sp] <- mean(OG_each_gene_score) #get the average
  }
}
ABscore_ingroups <- ABscore #ABscore_ingroup contain average Cscore for species in sp_list.
ABscore_outgroups <- ABscore
plot(ABscore_ingroups,ABscore_outgroups)


#############TAD
hicFindRestSite --fasta $genome --searchPattern GATC -o ${genome}_hic.bed
hicBuildMatrix --samFiles ./1bam/mate_R1.bam ./1bam/mate_R2.bam --binSize 20000 \
                 --restrictionSequence GATC --danglingSequence GATC \
                 --restrictionCutFile ${genome}_hic.bed \
                 --threads 30 --inputBufferSize 20000 -o ./2matrix/50k/50k.h5 --QCfolder ./hicQC \
                 --outBam ./1bam/hic.bam
#normalization
hicCorrectMatrix correct -m 20k.h5 -o 20k.corrected.h5

#plot contact matrix
hicPlotMatrix -m 50k.corrected.h5 -o 50k.corrected.h5.png --dpi 300 --log1p --perChromosome
#Find TAD
hicFindTADs -m 20k.corrected.h5 --outPrefix 20k.corrected.h5.TAD --correctForMultipleTesting fdr -p 24 \
--minDepth 60000 --maxDepth 200000 --thresholdComparisons 0.05 --delta 0.01 \
--TAD_sep_score_prefix myHiCmatrix_min10000_max40000_step1500_thres0.05_delta0.01_fdr 
make_tracks_file --trackFiles corrected.20k.h5 TAD_domains.bed -o tracks3.ini
hicPlotTADs --tracks tracks3.ini --region gm_1:14000000-20000000 -o TAD.png


