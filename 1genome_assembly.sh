####################################################
################short read cleanning
read1=/path/to/read1.fastq.gz
read2=/path/to/read2.fastq.gz
fastp -i=${read1} -I=${read1} -o=${read1}.clean.gz -O=${read1}.clean.gz -w=10

####################################################
################contig assembly
conda activate nextdenovo
nextDenovo run.cfg

####################################################
################purge_dups
conda activate purge_dups
###input data, draft genome and reads of PacBio
genome=/path/to/contigs.fasta
reads=/path/to/ont.fastq.gz
type=ont #ont,pb,hifi
cpu=32
minimap2 -t ${cpu} -x map-${type} ${genome} ${reads} | pigz -p 16 -c - > ont_aln.paf.gz
pbcstat ont_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
split_fa ${genome} > asm.split
minimap2 -t 32 -x asm5 -DP asm.split asm.split | pigz -p 4 -c > asm.split.self.ont.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.ont.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed ${genome}

####################################################
################Hi-C
conda activate yahs
contigs=/path/to/purged/contigs/purged.fa
r1=/path/to/hic_R1.fq.gz
r2=/path/to/hic_R2.fq.gz
juicer_jar=/path/to/juicer_tools_1.19.02.jar
chromap -i -r $contigs -o ${contigs}.index
chromap -t 32 -r $contigs -x ${contigs}.index --preset hic --remove-pcr-duplicates -1 $r1 -2 $r2 --trim-adapters --SAM -o hic.sam
samtools view test.sam -O BAM > hic.bam && rm hic.sam
samtools index hic.bam
samtools faidx $contigs
yahs -e GATC $contigs hic.bam
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp ${contigs}.fai >out_JBAT.log 2>&1
(java -jar -Xmx80G ${juicer_jar} pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
#using juicertools for manual correction
#after correct
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp $contigs

####################################################
################Polish
conda activate nextgenome
#long read polish
input=/path/to/chromosome/after/hic/out_JBAT.FINAL.fa
round=2
threads=32
read=/path/to/ont.fastq.gz
read_type=ont #{clr,hifi,ont}, clr=PacBio continuous long read, hifi=PacBio highly accurate long reads, ont=NanoPore 1D reads
read1=/path/to/short_read1.fq.gz
read2=/path/to/short_read2.fq.gz
mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
for ((i=1; i<=${round};i++)); do
    minimap2 -ax ${mapping_option[$read_type]} -t ${threads} ${input} ${read}| samtools sort - -m 1g --threads ${threads} -o lgs.sort.bam;
    samtools index lgs.sort.bam;
    ls `pwd`/lgs.sort.bam > lgs.sort.bam.fofn;
    python /mnt/data/miniconda3/envs/nextgenome/share/nextpolish-1.4.1/lib/nextpolish2.py -g ${input} -l lgs.sort.bam.fofn -r ${read_type} -p 4 -sp -o genome.long_polish.fa;
    if ((i!=${round}));then
        mv genome.long_polish.fa genome.longpolishtmp.fa;
        input=genome.longpolishtmp.fa;
    fi;
done
input=genome.long_polish.fa

#short read polish
input=genome.long_polish.fa
rounds=2
threads=32
for ((i=1; i<=${rounds};i++)); do
#step 1:
   bwa-mem2 index ${input};
   bwa-mem2 mem -t ${threads} ${input} ${read1} ${read2}| samtools view --threads 5 -F 0x4 -b -|samtools fixmate -m --threads 5  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   python /mnt/data/miniconda3/envs/nextgenome/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 1 -p 2 -s sgs.sort.bam > genome.polishtemp.fa;
   input=genome.polishtemp.fa;
#step2:
   bwa-mem2 index ${input};
   bwa-mem2 mem -t ${threads} ${input} ${read1} ${read2}| samtools view --threads 5 -F 0x4 -b -|samtools fixmate -m --threads 5  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   python /mnt/data/miniconda3/envs/nextgenome/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 2 -p 2 -s sgs.sort.bam > genome.nextpolish.fa;
   input=genome.nextpolish.fa;
done

####################################################
################BUSCO
conda activate busco
cpu=32
genome=genome=/path/to/genome.fasta
busco_database=busco_database=/path/to/lepidoptera_odb10
busco -m geno --cpu=$cpu -i $genome -l ${busco_database} -f --offline -o ${genome}_BUSCO












