cpu=32
genome=/path/to/genome.fasta
#fastq files have suffix like: *.R1.fq.gz, *.R2.fq.gz
#all sample names were list in a file sample.list, one per line
#fastq files were put in folder ./0raw

mkdir -p ./1filter ./2bam ./3gvcf

###################################mapping using bwa-mem2
#read filter
parallel -a sample.list --progress -j $cpu "fastp -i ./0raw/{}.R1.fq.gz -I ./0raw/{}.R2.fq.gz -o ./1filter/{}.R1.fq.gz -O ./1filter/{}.R2.fq.gz -w 1"
#mapping
bwa-mem2 index $genome
parallel -a sample.list -j 1 \
"bwa-mem2 mem -R '@RG\tID:{}\tSM:{}' -t $cpu $genome ./1filter/{}.R1.fq.gz ./1filter/{}.R2.fq.gz | samtools sort -m 1g -O BAM -o ./2bam/{}.bam"
#index
parallel -a list_ind -j5 "samtools index ./2bam/{}.bam"

#####################################call snp using Freebayes
find ./2bam/ -name *bam > bamlist
freebayes -f ${genome} --limit-coverage 50 --min-coverage 10 -k --use-mapping-quality --use-best-n-alleles 4 -F 0.2 -C 2 -0 -A cnv_map --populations pop_map\
--haplotype-length -1 --report-monomorphic --ploidy 2 --bam-list bamlist {} | bgzip > ./3gvcf/gvcf.vcf.gz
bcftools index ./3gvcf/gvcf.vcf.gz
bcftools view -e 'F_MISSING > 0.5' ./3gvcf/gvcf.vcf.gz -Oz > ./3gvcf/gm.gvcf.gz

#####################################VCF processing
#haploid2diploid
bgzip -cd gm.gvcf.gz |  sed -r -e 's@\t([0-9.]):@\t\1/\1:@g' | bcftools view -Oz -o gm.gvcf.dip.gz
bgzip -cd gd.gvcf.gz |  sed -r -e 's@\t([0-9.]):@\t\1/\1:@g' | bcftools view -Oz -o gd.gvcf.dip.gz

#intersection between two species
bcftools isec gm.gvcf.dip.gz gd.gvcf.dip.gz -p gdgm_isec -Oz
bcftools merge ./gdgm_isec/0002.vcf.gz ./gdgm_isec/0003.vcf.gz -Oz > gmgd.gvcf.gz

###################################fst and pi
parseVCF.py  -i gmgd.gvcf.gz -o gmgd.snp.vcf.gz.geno.gz
popgenWindows.py -w 100000 -s 100000 -g gmgd.snp.vcf.gz.geno.gz --windType coordinate -o output.100k.csv -f phased -T 32  --popsFile popmap -p gm -p gd --writeFailedWindows




#filter
bcftools view gmgd.gvcf.gz -v snps -m2 -M2 -Oz > gmgd.snp.vcf.gz




