#put all mutiple sequance alignments .fasta in work_dir
work_dir=/mnt/dataset/moth/Analysis/dsdn/all_leafrooler/abcc1_aln

#11codeml.sh
realpath ${work_dir}/*fasta | parallel   -j32 "./codeml.sh {}" 
#11chi2.sh
realpath ${work_dir}/*fasta | parallel   -j32 "./chi2.sh {}"
