####ATAC peak
mamba activate atac
genome=/mnt/dataset/moth/species/gm/2genome/gm.fasta
r1=/mnt/dataset/moth/species/gm/12atac/0raw/F-1_R1.clean.fq.gz
r2=/mnt/dataset/moth/species/gm/12atac/0raw/F-1_R2.clean.fq.gz
r1=/mnt/dataset/moth/species/gm/12atac/0raw/M-1_R1.clean.fq.gz
r2=/mnt/dataset/moth/species/gm/12atac/0raw/M-1_R2.clean.fq.gz

chromap -i -r $genome -o ${genome}.index
chromap -t 32 -r $genome -x ${genome}.index --preset atac --remove-pcr-duplicates -1 $r1 -2 $r2 --SAM -o ./2bam/Male.sam
samtools sort -m 4g -@5 -O BAM > ./2bam/Male.bam && rm ./2bam/Male.sam
#call peak
macs3 callpeak -f BAMPE -t ./2bam/Male.bam ./2bam/Female.bam -g 5e8 -n ./3peak/FM -B -q 0.01

#calculate density
bedtools makewindows -g ${genome}.fai -w 50000 > 50kb_win.bed
bedtools coverage -a 50kb_win.bed -b ./3peak/FM_accessible_regions.gappedPeak.bed > 50kb_atac.density.bedgraph

####CpG methylation
mamba activate nanopolish
genome=/mnt/dataset/moth/species/gm/2genome/gm.fasta
ont_summary=/mnt/dataset/moth/species/gm/0assembly/0raw/ont/20181029-NPL0393-P2-A3-D3.sequencing_summary.txt
ont_read=/mnt/dataset/moth/species/gm/0assembly/0raw/ont/20181029-NPL0393-P2-A3-D3.fastq.gz
nanopolish index -v \
-d /mnt/dataset/moth/species/gm/9CpG/fast5/0 \
-d /mnt/dataset/moth/species/gm/9CpG/fast5/1 \
-d /mnt/dataset/moth/species/gm/9CpG/fast5/2 \
-d /mnt/dataset/moth/species/gm/9CpG/fast5/3 \
-d /mnt/dataset/moth/species/gm/9CpG/fast5/4 \
.
.
.
-d /mnt/dataset/moth/species/gm/9CpG/fast5/1801 \
-s $ont_summary $ont_read

minimap2 -t 30 -a -x map-ont $genome $ont_read | samtools sort -m 2g -O BAM  > ont.bam
samtools index ont.bam

nanopolish call-methylation -t 32 -r $ont_read -b ont.bam -g $genome -q cpg > ./gm_all/methylation.gm_all.tsv
calculate_methylation_frequency.py -c 10 ./gm_all/methylation.gm_all.tsv > ./gm_all/methylation_frequency.gm_all.tsv

#calculate CpG frequency in 50kb windows
bedtools coverage -a 50kb_win.bed -b ./3peak/FM_accessible_regions.gappedPeak.bed > 50kb_atac.density.bedgraph
tail -n +2 ./gm_all/methylation_frequency.gm_all.tsv > ./gm_all/methylation_frequency.gm_all.tsv2
bedtools sort -i ./gm_all/methylation_frequency.gm_all.tsv2 > ./gm_all/methylation_frequency.gm_all.bed
bedtools sort -i 50kb_win.bed > 50kb_win.sort.bed
bedtools map -o sum -c 5,6 -a win50k.bed -b methylation_frequency.gm_all.tsv2 | awk -F'\t' '{print $1,$2,$3,$5/$4}' >./gm_all/methylation_frequency.gm_all.bedgraph
