
msa=$1
v1=$(grep lnL ${msa}_BSM/modelA/A_mlc         | awk '{print $5}')
v2=$(grep lnL ${msa}_BSM/modelAnull/Anull_mlc | awk '{print $5}')
c=$(echo "2 * ($v1 - $v2)" | bc)
lrt=$(echo "if ($c < 0) -($c) else $c" | bc)
chi2 2 $lrt >  ${msa}_BSM/lrt
grep prob ${msa}_BSM/lrt | awk '{print $6}' > ${msa}_BSM/p
