#!/bin/bash

function usage
{
	echo "Usage:"
	echo "run_all_svm.sh -k <kmer_lengths> [-p <svm_parameters>] [-t <testing_dir>] [-c <cores>]"
	echo "Defaults:"
	echo "p='svm_pca,gamma=0.005,C=4,n=131'"
	echo "t='.'"
	echo "c=10"
}

## Managing command line options.

# Default: parameters that seem to work fine
params='svm_pca,gamma=0.005,C=4,n=131'
testdir='.'
cores=10

# Check to see if at least two arguments were specified
if [ $# -lt 1 ] ; then
	echo "You need at least a list of k values"
	usage
	exit 1
fi

# Get the command line options
while getopts k:p:t:c: opt
do
   case "$opt" in
      k) kvals=$OPTARG;;
      p) params=$OPTARG;;
      t) testdir=$OPTARG;;
      c) cores=$OPTARG;;
      \?) echo "Option(s) wrongly formed"
			usage
			;;
   esac
done

# if no k vals specified abort
if [ -z "$kvals" ] ; then
	echo "You need a list of k values"
	usage
	exit 1
fi

## Run, Forrest

for k in $kvals
do
	#python ../runsvm.py matrix_k${k}_arch_bact_euk_virus.bdat ${testdir}/matrix_k${k}_arch_bact_euk_virus.bdat $params $cores
	#python ../runsvm.py matrix_k${k}_arch_bact_euk_virusfiltered.bdat ${testdir}/matrix_k${k}_arch_bact_euk_virus.bdat $params $cores
	python ../epluche_svm.py ${testdir}/${params}-k${k}_arch_bact_euk_virus.csv
	#python ../epluche_svm.py ${testdir}/${params}-k${k}_arch_bact_euk_virusfiltered.csv
done







