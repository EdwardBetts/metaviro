#!/bin/bash

## Some functions. We all like functions.

function usage
{
	echo "Usage:"
	echo "run_basio.bash [-p '<float>'] [-k <int>] [-d <shared dir>] -f <multifasta file path>"
	echo "p =  Border Insertion Penalty / k = number of seq per job"
	echo "Defaults:"
	echo "p='0.0'"
	echo "k='1000'"
	echo "d='/tmp'"
	echo "Note:"
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
bip='0.0'
pack='1000'
shareddir='/tmp'

# Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
   usage
   exit 1
fi

# Get the command line options
while getopts p:k:f:d: opt
do
   case "$opt" in
      p) bip=$OPTARG;;
      k) pack=$OPTARG;;
      f) fastafile=$OPTARG;;
      d) shareddir=$OPTARG;;
      \?) usage;;
   esac
done

# if no file specified abort
if [ -z "$fastafile" ] ; then
	echo "You need an input fasta file" \"$fastafile\" not valid
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
bindir=`dirname $(readlink -f $0)`/bin

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
grep "^>" $multifasta | sed -e 's/^>//' > $shareddir/$runid/ids.txt

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
	$bindir/basio_splitb /tmp/$runid/\${header}.fasta
	$bindir/basio_basio --bip=$bip /tmp/$runid/\${header}.sgm
	rm /tmp/$runid/\${header}.fasta
	echo ">\${header}" >> /tmp/$runid/job\${PBS_JOBID}.sgm
	cat /tmp/$runid/\${header}.sgm >> /tmp/$runid/job\${PBS_JOBID}.sgm
	echo >> /tmp/$runid/job\${PBS_JOBID}.sgm
	rm /tmp/$runid/\${header}.sgm
	shift
done

if [ "$shareddir" != "/tmp" ]
then
	mv /tmp/$runid/job\${PBS_JOBID}.sgm $shareddir/$runid/job\${PBS_JOBID}.sgm
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
if [ -s ${prefix}_bip${bip}.sgm ]
then
	rm ${prefix}_bip${bip}.sgm
fi
for sgmfile in $shareddir/$runid/job*.sgm
do
	cat $sgmfile >> ${prefix}_bip${bip}.sgm 
	rm $sgmfile 
done

# clean /tmp (or wherever the shared dir is)
cleanlogs
rm $shareddir/$runid/ids.txt
rmdir $shareddir/$runid 
