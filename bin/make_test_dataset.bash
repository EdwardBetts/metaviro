#!/bin/bash

# Note: the data directory must contain the reference genomes as follows:
# All the families that should be computed must be concatenated in a
# file that is named after the root directory of said family. For instance
# the family X, contained in a/b/c/X must have a file a/b/c/X/X.fasta
# that contains all the genomes of the family in fasta format.
# On the other hand, the resulting resampling will be named a/b/c/X.fasta
# and if there are subfamilies that must also be computed, they will be
# stored in a a/b/c/X directory, e.g. a/b/c/X/Y.fasta
# Yeah, just follow the first rule and see what happens, should be
# easier to understand that way.

function usage
{
	echo "Usage:"
	echo "make_test_dataset.bash -d <path> [-p <path>] [-n <int>] [-l <int>]"
	echo "d = data root path/ p = prefix / n = number of fake contigs to create per fasta / l = maximum depth level"
	echo "Defaults:"
	echo "p='' (the prefix will be glued as is to the standard filename, so you can give full paths too if you don't want the results in the working directory)"
	echo "n='10000'"
	echo "l='100' (0 -> first level only)"
}

## Managing command line options.

# Default: resamplings of 10k contigs

ncontigs='10000'
depth=100

# Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
	echo "You need at least one argument specified"
   usage
   exit 1
fi

# Get the command line options
while getopts d:p:n:l: opt
do
   case "$opt" in
      d) dataroot=$OPTARG;;
      p) prefix=${OPTARG}_;;
      n) ncontigs=$OPTARG;;
      l) depth=$OPTARG;;
      \?) echo "Wrong options."
			usage
			;;
   esac
done

# if no file specified abort
if [ "$dataroot" = '' ] ; then
	echo "You need to specify a data root directory"
   usage
   exit 1
fi

# If the prefix is a relative path, we need to "absolutify" it before
# changing the current dir
if [[ $prefix != "/"* ]] ; then
	prefix=`pwd`/$prefix
fi

dataroot=`pwd`/${dataroot/%\//} # cleaning data path
datadest=${prefix}`basename $dataroot`_n${ncontigs}_l$depth # creating ref path
bindir=`dirname $(readlink -f $0)`

echo 'Reference data root path:' $dataroot
echo 'Destination root path:' $datadest

# Checking if the dataset already exists
if [ -d $datadest ]
then
	echo "Dataset already exists, continue anyway? (y/N)"
	read ans
	if [ "$ans" = 'y' ]
	then
		echo "Overwriting..."
	else
		echo "Aborting."
		exit 0
	fi
else
	mkdir -p $datadest
fi

module load python/2.7.3

cd $dataroot # just because, it's easier this way (find being a general pain in the butt and me being lazy)

# Look for fasta file within the depth limit (0 -> first level only)
for taxdir in `find * -maxdepth $depth -type d`
do
	if [ -s ${taxdir}/$(basename $taxdir).fasta ]
	then
		echo 'Current file:' `ls ${taxdir}/$(basename $taxdir).fasta`
		if [ ! -d ${datadest}/`dirname $taxdir` ]
		then
			# Make destination dir if needed
			echo '- mkdir' ${datadest}/`dirname $taxdir`
			mkdir -p ${datadest}/`dirname $taxdir`
		fi
		# Resample!
		echo '- input' ${taxdir}/$(basename $taxdir).fasta
		echo '- output' ${datadest}/${taxdir}.fasta
		# (remove previous resamples / fasta indices if needed)
		if [ -s ${datadest}/${taxdir}.fasta ]
		then
			rm ${datadest}/${taxdir}.fasta
			rm ${datadest}/${taxdir}.fasta.fai
		fi
		python $bindir/generate_fake_contigs.py -n $ncontigs -o ${datadest}/${taxdir}.fasta ${taxdir}/$(basename $taxdir).fasta
		# Compute k-mers
		#~ for k in 3 4 5 6 7 8
		#~ do
			#~ echo 'k =' $k
			#~ ~/metaviro2/data_generation/run_kmer.bash  -k $k -f ${datadest}/${taxdir}.fasta
			#~ python ~/metaviro2/data_generation/convert_kmer.py -t 20 ${datadest}/${taxdir}_k${k}.frags
			#~ rm ${datadest}/${taxdir}_k${k}.frags
		#~ done
		#~ # Compute segments
		#~ for bip in 0.0 0.2 0.4
		#~ do
			#~ echo 'bip =' $bip
			#~ ~/metaviro2/data_generation/run_basio.bash  -p $bip -f ${datadest}/${taxdir}.fasta
			#~ python ~/metaviro2/data_generation/convert_basio.py -t 20 ${datadest}/${taxdir}_bip${bip}.sgm
			#~ rm ${datadest}/${taxdir}_bip${bip}.sgm
		#~ done
	fi
done

echo "Resampling done."
