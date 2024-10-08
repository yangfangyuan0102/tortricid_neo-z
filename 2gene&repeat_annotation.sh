cpu=32
species=gm
work_dir=/path/to/species
genome=/path/to/genome.fasta

####################################################
################Transcriptome assembly
cd $work_dir
hisat2-build $genome $genome
hisat2 --dta -p $cpu -x $genome -1 1read1.fq.gz -2 1read2.fq.gz | samtools sort -@ 10 -m 1g -O BAM > 1trans.bam
hisat2 --dta -p $cpu -x $genome -1 2read1.fq.gz -2 2read2.fq.gz | samtools sort -@ 10 -m 1g -O BAM > 2trans.bam
hisat2 --dta -p $cpu -x $genome -1 3read1.fq.gz -2 3read2.fq.gz | samtools sort -@ 10 -m 1g -O BAM > 3trans.bam

#bam2gtf
stringtie ${species}.1trans.bam -o ${species}.1trans.bam.gtf -p 10
stringtie ${species}.2trans.bam -o ${species}.2trans.bam.gtf -p 10
stringtie ${species}.3trans.bam -o ${species}.3trans.bam.gtf -p 10
#merge bam
find ./1trans/ -name *bam.gtf > gtf.list
stringtie --merge gtf.list -o trans.gtf
find ./1trans/ -name *trans.bam > bam.list
samtools merge -f -b ./bam.list -@ 15 trans.bam
rm gtf.list bam.list
gffread trans.gtf -g $genome -w trans.fasta
gffread trans.gtf > trans.gff3

####################################################
################Repeat annotation
conda activate repeat

species=$1
cpu=10 #cpu/4
work_dir=/mnt/dataset/moth/species/${species}
genome=/path/to/genome.fasta
repeatmasker_lib=/database/repeat_masker_Arthropoda.lib
mkdir -p ${work_dir}/3anno/repeat/modeler
mkdir -p ${work_dir}/3anno/repeat/masker
#repeatmodeler
cd ${work_dir}/3anno/repeat/modeler
BuildDatabase -name ${species}.db ${genome}
RepeatModeler -database ./${species}.db -pa $cpu -LTRStruct
#merge the database of repeatmodeler and repeatmasker
cp ./${species}.db-families.fa ../${species}.db-families.fa
cat ./${species}.db-families.fa ${repeatmasker_lib} > ../${species}.repeat.lib
#repeatmasker
cd ${work_dir}/3anno/repeat/masker
RepeatMasker  -no_is -norna -xsmall -q -lib ../ar.repeat.lib -parallel ${cpu} -gff -dir ./ ${genome}

####################################################
################PASA gene prediction
conda activate pasa
species=gm
cpu=32
pasaHome=/mnt/data/miniconda3/envs/pasa/opt/pasa-2.5.2
work_dir=/path/to/${species}
genome=${work_dir}/2genome/${species}.fasta
vec_data=/path/to/UniVec.fasta #Contamination database
transcripts=/path/to/trans.fasta
trans_gff3=/path/to/trans.gff3

mkdir -p ${work_dir}/3anno/pasa
cd ${work_dir}/3anno/pasa
#filter transcriptome
${pasaHome}/bin/seqclean $transcripts -v $vec_data -o ${transcripts}.clean -r ${transcripts}.cln -c 2
awk -F " " '{print $1}' ${transcripts}.cln > FL_accs.txt
#configure
cp ${pasaHome}/pasa_conf/pasa.alignAssembly.Template.txt ./alignAssembly.config
cp ${pasaHome}/pasa_conf/pasa.annotationCompare.Template.txt ./annotCompare.config
sed -i "s@DATABASE=.*@DATABASE=${work_dir}/3anno/pasa/${species}.pasa.sqlite@" *.config
#run pasa
${pasaHome}/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $genome -t ${transcripts}.clean -T -u ${transcripts} -f FL_accs.txt --stringent_alignment_overlap 30.0 --MAX_INTRON_LENGTH 60000 --TRANSDECODER --ALIGNER blat --CPU ${cpu}
${pasaHome}/scripts/build_comprehensive_transcriptome.dbi -c ./alignAssembly.config -t ${transcripts}.clean
${pasaHome}/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta *.assemblies.fasta --pasa_transcripts_gff3 *.pasa_assemblies.gff3
#extract predected protein sequences longer than >=300，for augustus or SNAP
#${pasaHome}/scripts/pasa_asmbls_to_training_set.extract_reference_orfs.pl *.assemblies.fasta.transdecoder.genome.gff3 300 > ./best_candidates.gff3

####################################################
################DeNovo gene prediction using Helixer

disk_path=/mnt/dataset
sudo docker run --runtime=nvidia -it --name helixer_testing_v0.3.2_cuda_11.2.0-cudnn8 --rm --shm-size=50G \
-v ${disk_path}:${disk_path} gglyptodon/helixer-docker:helixer_v0.3.2_cuda_11.8.0-cudnn8

species=gm
work_dir=/path/to/${species}
genome=${work_dir}/2genome/${species}.fasta
model=/path/to/invertebrate_v0.3_m_0200.h5

Helixer.py --lineage invertebrate --model-filepath $model --subsequence-length 320760 \
--compression lzf --batch-size 6 --fasta-path ${genome} --species $species --gff-output-path ${work_dir}/helixer.gff3

####################################################
####################miniprot
miniprot -t1 --gff-only --gff $genome Plutella_xylostella.faa > miniprot.gff3
miniprot_GFF_2_EVM_GFF3.py miniprot.gff3 > miniport_evm.gff3


####################################################
################Evidencemodeler
conda activae evidence
species=$1
cpu=32
work_dir=/mnt/dataset/moth/species/${species}
genome=${species}.fasta

$EVM_HOME/EVidenceModeler --sample_id ly \
--genome $genome --gene_predictions gene_predictions.gff3 \
--protein_alignments miniport_evm.gff3 \
--transcript_alignmentspasa.sqlite.pasa_assemblies.gff3 \
--segmentSize 600000 --overlapSize 30000 --CPU 32 \
--repeats repeatmasker.gff --weight weight


####################################################
################BUSCO

conda activate busco
cpu=32
gff=/path/to/annoted/species.evm.all.gff3
genome=/path/to/genome.fasta
busco_database=/path/to/lepidoptera_odb10
gffread $gff -g $genome -y ${gff}.protein
busco -m prot --cpu=$cpu -i ${gff}.protein -l ${busco_database} -f --offline -o ${gff}.busco 





