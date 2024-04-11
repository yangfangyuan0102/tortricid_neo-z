
#run the branch-site model
msa=$1
tree=/mnt/dataset/moth/Analysis/dsdn/Orth_all/OrthoFinder/24species/paml/tree.root.nwk.nolength.Labeled.txt
modelA_temp=/mnt/dataset/moth/Analysis/dsdn/Orth_all/OrthoFinder/24species/paml/modelA.ctl
modelAnull_temp=/mnt/dataset/moth/Analysis/dsdn/Orth_all/OrthoFinder/24species/paml/modelAnull.ctl

mkdir -p ${msa}_BSM/modelA ${msa}_BSM/modelAnull

sed "s#seqfile = .*#seqfile = ${msa}#" ${modelA_temp} | sed "s#treefile.*#treefile = ${tree}#"| sed "s#outfile = .*#outfile = ${msa}_BSM/modelA/A_mlc#"    > ${msa}_BSM/modelA/ctl
sed "s#seqfile = .*#seqfile = ${msa}#" ${modelAnull_temp} | sed "s#treefile.*#treefile = ${tree}#"| sed "s#outfile = .*#outfile = ${msa}_BSM/modelAnull/Anull_mlc#"> ${msa}_BSM/modelAnull/ctl

cd ${msa}_BSM/modelA
codeml ctl
cd ${msa}_BSM/modelAnull
codeml ctl

