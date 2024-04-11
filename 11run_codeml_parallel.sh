#put all mutiple sequance alignments .fasta in work_dir
work_dir=/mnt/dataset/moth/Analysis/dsdn/all_leafrooler/abcc1_aln

#11codeml.sh
find ${work_dir} -type f -name "*fasta" | parallel   -j32 "./codeml.sh {}" 
#11chi2.sh
find ${work_dir} -type f -name "*fasta" | parallel   -j32 "./chi2.sh {}"
