#!/bin/bash


virusfile=$1
shift
livingcreatures=$*

module load blast

for reffile in $livingcreatures ; do
	formatdb -i $reffile -p F -t title
	megablast -d $reffile -i $virusfile -o $reffile.blastout -e 0.001 -D 2 -W 16
	#~ blastall -e 0.1 -p blastn -d $reffile -i $virusfile -o $reffile.blastout
done

#~ blastall -e 10 -p blastx -d nr -i $virusfile -o $virusfile.nr.blastout
#~ blastall -e 10 -p blastn -d nt -i $virusfile -o $virusfile.nr.blastout


