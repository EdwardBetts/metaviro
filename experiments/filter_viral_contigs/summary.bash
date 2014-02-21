datadir="$1/"

declare -a orgs=(arch bact euk virus)

echo "summaries for kmers..."

params="k3 k4" #k5 k6 k7 k8"

for dist in rof #js kl
do
	for query in arch bact euk virus
	do
		echo "label $params" > ${datadir}${query}_kmer_${dist}.csv
		declare -a arch=("arch")
		declare -a bact=("bact")
		declare -a euk=("euk")
		declare -a virus=("virus")
		for p in $params
		do
			  arch=("${arch[@]}" `tail -n 1 ${datadir}${query}_cont_${p}.${dist} | cut -f2`)
			  bact=("${bact[@]}" `tail -n 1 ${datadir}${query}_cont_${p}.${dist} | cut -f3`)
			  euk=("${euk[@]}" `tail -n 1 ${datadir}${query}_cont_${p}.${dist} | cut -f4`)
			  virus=("${virus[@]}" `tail -n 1 ${datadir}${query}_cont_${p}.${dist} | cut -f5`)
		done
		echo ${arch[@]} >> ${datadir}${query}_kmer_${dist}.csv
		echo ${bact[@]} >> ${datadir}${query}_kmer_${dist}.csv
		echo ${euk[@]} >> ${datadir}${query}_kmer_${dist}.csv
		echo ${virus[@]} >> ${datadir}${query}_kmer_${dist}.csv
	done
done

echo "done."
pwd
