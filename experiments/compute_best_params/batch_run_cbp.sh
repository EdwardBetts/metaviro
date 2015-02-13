#!/bin/bash

script=$1
shift
kvals=$*

output=$(basename $script)_k${kvals// /}.results
echo "Results:" > $output

module load python/2.7.3 

for k in $kvals
do
	echo "k =" $k >> $output
	python $script matrix_k${k}_arch_bact_euk_virus.bdat > matrix_k${k}_arch_bact_euk_virus.bestparams >> $output
done

