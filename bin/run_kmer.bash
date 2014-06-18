#!/bin/bash

## Some functions. We all like functions.

function usage
{
	echo "Usage:"
	echo "run_kmer.bash [-k <int>] [-d <shared dir>] -f <multifasta file path>"
	echo "k = size of k-mers / b = size of batches of sequences per job"
	echo "Defaults:"
	echo "k='3'"
	echo "b='1000'"
	echo "d='/tmp'"
}

function cleanlogs
{
	if [ -n "`ls $shareddir/$runid/job*.err 2> /dev/null`" ]
	then
		for errfile in $shareddir/$runid/job*.err
		do
			if [ ! -s $errfile ]
			then
				errfileprefix="${errfile%.*}"
				rm $errfileprefix.err $errfileprefix.out
			fi
		done
	fi
}

function waitjobs
{
	qstat -x | grep $runid | wc -m
	while [ `qstat -x | grep $runid | wc -m` -gt 0 ]
	do
		echo "Waiting for jobs to finish..."
		sleep 5
		cleanlogs
	done
}

## Managing command line options.

# Default: no border insertion penalty for basio and 1000 sequences per job
# and shared directory is /tmp (if it runs on a single machine)
k=3
pack='1000'
shareddir='/tmp'

# Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
	echo "You need at least one argument specified"
   usage
   exit 1
fi

# Get the command line options
while getopts k:b:f:d: opt
do
   case "$opt" in
      k) k=$OPTARG;;
      b) pack=$OPTARG;;
      f) fastafile=$OPTARG;;
      d) shareddir=$OPTARG;;
      \?) echo "Option(s) wrongly formed"
			usage
			;;
   esac
done

# if no file specified abort
if [ -z "$fastafile" ] ; then
	echo "You need an input fasta file"
   usage
   exit 1
fi

## Preparing the environment.

# just to have nice variables for the paths (with lame names but meh)
multifasta=`readlink -f $fastafile`
fulldir=`dirname $multifasta`
filename=`basename $fastafile`
extension="${filename##*.}"
filenameprefix="${filename%.*}"
prefix=$fulldir/$filenameprefix
bindir=`dirname $(readlink -f $0)`

# Set a name for the run, using a random number, to name the jobs
runid=${filenameprefix}_$RANDOM
echo "Run ID = $runid"

# and the directory to put the job results in before final merging
mkdir $shareddir/$runid

# create the multifasta index (once and for all, we'll use it a lot)
if [ ! -s $multifasta.fai ]
then
	echo "Indexing multifasta file..."
	qsub <<EOS
#!/bin/sh

#PBS -l nodes=1:ppn=1,walltime=00:05:00
#PBS -N $runid
#PBS -o $shareddir/$runid/job\${PBS_JOBID}.out
#PBS -e $shareddir/$runid/job\${PBS_JOBID}.err
module load samtools
samtools faidx $multifasta
EOS
	waitjobs
fi

## Jobs!

# Create sequence ids list
grep "^>" $multifasta | sed -e 's/^>//' | cut -d" " -f1 > $shareddir/$runid/ids.txt

# ( °_°)/⌒●~*
echo "Launching jobs..."
for ((I=1; I <= `wc -l $shareddir/$runid/ids.txt | cut -d" " -f 1`; I+=$pack))
do
	qsub <<EOS
#!/bin/sh

#PBS -l nodes=1:ppn=1,walltime=00:05:00
#PBS -N $runid
#PBS -o $shareddir/$runid/job\${PBS_JOBID}.out
#PBS -e $shareddir/$runid/job\${PBS_JOBID}.err
module load samtools

# Moving to /tmp so that if we're on a cluster we don't mess with the NFS
if [ "$shareddir" != "/tmp" ]
then
	mkdir /tmp/$runid
fi

# Segmenting!
sed -n ${I},$(( $I + $pack - 1 ))p $shareddir/$runid/ids.txt | while read header
do
	samtools faidx $multifasta \$header > /tmp/$runid/\${header}.fasta
	echo $multifasta \$header
	$bindir/kmer $k /tmp/$runid/\${header}.fasta
	rm "/tmp/$runid/\${header}.fasta"
	echo ">\${header}" >> /tmp/$runid/job\${PBS_JOBID}.frags
	cat /tmp/$runid/\${header}.frags >> /tmp/$runid/job\${PBS_JOBID}.frags
	rm /tmp/$runid/\${header}.frags
	shift
done

if [ "$shareddir" != "/tmp" ]
then
	mv /tmp/$runid/job\${PBS_JOBID}.frags $shareddir/$runid/job\${PBS_JOBID}.frags
	rm -R /tmp/$runid
fi

echo 'Done :)'

EOS
	# Clean up a little if there's no error messages in the logs
	cleanlogs
	# wait a little so PBS doesn't go nuts
	sleep 0.1
done

# wait for all jobs to finish keeping everything tidy
waitjobs

## Post-processing.

# merge results (maybe a bit ugly? too many files?)
if [ -s ${prefix}_k$k.frags  ]
then
	rm ${prefix}_k$k.frags 
fi
for sgmfile in $shareddir/$runid/job*.frags
do
	cat $sgmfile >> ${prefix}_k$k.frags 
	rm $sgmfile 
done

# clean /tmp (or wherever the shared dir is)
cleanlogs
rm $shareddir/$runid/ids.txt
rmdir $shareddir/$runid 
