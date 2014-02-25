#! /usr/bin/env python
# encoding: utf-8
import itertools
# import prettytable
import sys
import collections
from Bio.Seq import Seq
from Bio.SeqUtils import CheckSum
import scipy
from Bio import SeqUtils
from Bio import SeqIO

from Bio.Alphabet import generic_dna

import sys
import os
import argparse
# from optparse import OptionParser 

import logging
if "logger" not in globals():
	logger = logging.getLogger('kmerize')
	logger.setLevel(logging.DEBUG)

	# while len(logger.handlers()) > 0:
	#  	logger.pop()

	# create console handler and set level to debug
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s',"%Y-%m-%d %H:%M:%S")
	# formatter = logging.Formatter('%(asctime)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)

	# add ch to logger
	logger.addHandler(ch)

help="""Usage 
%(tool)s  FASTAFILE FASTAFILE 

Options
========

Example usage
==============
%(tool)s genome1.fa genome2.fa

Todo
======
* Long/Short format: DONE
* Add wildcard in pos 3 

"""%{"tool":sys.argv[1]}


def find_fasta_files(src_dir):
	all_fastas=[]
	for dirname, dirnames, filenames in os.walk(src_dir):
	 # print path to all filenames.
		for filename in filenames:
			if filename.endswith(".fasta"):
				all_fastas.append(os.path.join(dirname, filename))
	return all_fastas

def generate_wide_table(all_fastas):
	global args
	basis=[['A','T','G','C']]*args.kmer_length
	all_kmers=sorted(["".join(x) for x in tuple(itertools.product(*basis))])

	records_to_kmer={}
	for f in all_fastas: 
		logger.debug("Processing file %s"%(f))
		for record in SeqIO.parse(f, "fasta", generic_dna):
			logger.debug("Processing sequence %s"%(record.description))

			seq=str(record.seq)
			fasta_keys=[f,record.description,str(len(seq))]
			# Add additional features to the sequence
			fasta_keys.append(str(SeqUtils.GC(record.seq)))

			# 
			fasta_keys=tuple(fasta_keys)
			records_to_kmer[fasta_keys]=collections.defaultdict(int)
			for i in range(0,len(seq)-args.kmer_length):
				kmer=seq[i:i+args.kmer_length]
				records_to_kmer[fasta_keys][kmer]+=1

	print >>args.outfile, "\t".join(["path","sequence_description","sequence_length","GC"]+all_kmers)
	for k,kmer_values in records_to_kmer.items():
		all_values = list(k)
		all_values.extend(map(str,[kmer_values.get(x,0) for x in all_kmers]))
		# print len(all_values)
		print >>args.outfile, "\t".join(all_values)

def generate_long_table(all_fastas):
	global args
	print >>args.outfile, "\t".join(["path","sequence_description","sequence_length","GC","kmer","count"])
	for f in all_fastas: 
		logger.debug("Processing file %s"%(f))
		for record in SeqIO.parse(f, "fasta", generic_dna):
			kmer_count=collections.defaultdict(int)

			logger.debug("Processing sequence %s"%(record.description))

			seq=str(record.seq)
			fasta_keys=[f,record.description,str(len(seq))]
			# Add additional features to the sequence
			fasta_keys.append(str(SeqUtils.GC(record.seq)))

			# 
			# fasta_keys=tuple(fasta_keys)
			# kmer_count[fasta_keys]=collections.defaultdict(int)
			for i in range(0,len(seq)-args.kmer_length):
				kmer=seq[i:i+args.kmer_length]
				kmer_count[kmer]+=1
			for kmer, count in kmer_count.items():
				print >>args.outfile,"\t".join(fasta_keys+[kmer,str(count)])

def generate_count_table(all_fastas,prefix):
	logger.info("Genrating counts")
	global args
	print >>args.outfile, "\t".join(["prefix","kmer","count"])
	kmer_count=collections.defaultdict(int)
	try:
		for f in all_fastas: 
			logger.debug("Processing file %s"%(f))
			for record in SeqIO.parse(f, "fasta", generic_dna):

				logger.debug("Processing sequence %s"%(record.description))

				seq=str(record.seq)
				# fasta_keys=[f,record.description,str(len(seq))]
				# # Add additional features to the sequence
				# fasta_keys.append(str(SeqUtils.GC(record.seq)))

				# 
				# fasta_keys=tuple(fasta_keys)
				# kmer_count[fasta_keys]=collections.defaultdict(int)
				for i in range(0,len(seq)-args.kmer_length):
					kmer=seq[i:i+args.kmer_length]
					if not args.keep_IUPAC:
						wrong=False
						for j in range(0, args.kmer_length):
							if not (kmer[j]=="A" or kmer[j]=="C" or kmer[j]=="G" or kmer[j]=="T"):
								wrong=True
								break
						if (wrong):
							continue
					if args.star:
						kmer=list(kmer)
						for i in range(2,args.kmer_length,3):
							kmer[i]="*"
						kmer="".join(kmer)
					kmer_count[kmer]+=1
	except KeyboardInterrupt:
		pass
	logger.info("Finished counting, outputting")
	for kmer, count in kmer_count.items():
		print >>args.outfile,"\t".join([prefix,kmer,str(count)])

	# for k,kmer_values in records_to_kmer.items():
	# 	all_values = list(k)
	# 	all_values.extend(map(str,[kmer_values.get(x,0) for x in all_kmers]))
	# 	# print len(all_values)
	# 	print >>args.outfile, "\t".join(all_values)



def main(argv=None):
	global args
	parser=argparse.ArgumentParser(description="Compute kmer counts for all sequences in the input multi fasta file")
	parser.add_argument('-k',dest="kmer_length",help="mer length",default=3,type=int)
	parser.add_argument('-r',dest="recursive",help="Perform recursive search for fasta files",action="store_true")
	parser.add_argument('-L',dest="long",help="Switch output from wide(default) to Long format(One fasta/kmer per line, sparse)",action="store_true")
	parser.add_argument('-S',dest="star",help="Generated starred k-mers, using * for each third nt",action="store_true")
	parser.add_argument('-I',dest="keep_IUPAC",help="Do not throw k-mer with IUPAC code nt",action="store_true")
	parser.add_argument('-C',default="None",type=str,dest="counts_prefix",help="Indicate prefix for count output, and switch output from wide(default) to count format(One kmer per line, count in all the FATSTAs in the input)")
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")


	parser.add_argument('FASTAFILE',action='append',nargs="+",help='list of fasta files')
	args=parser.parse_args()

	FASTAFILE=args.FASTAFILE[0]
	all_fastas=set()
	if args.recursive:
		for f in FASTAFILE:
			all_fastas.update(find_fasta_files(f))
	else:
		all_fastas=set(FASTAFILE)

	logger.info("Will process %d fasta files"%(len(all_fastas)))
	if(args.kmer_length>=10 and not (args.long or args.counts_prefix)):
		logger.critical("K-mer length %d too large for wide output format"%(args.kmer_length))
		sys.exit(1)

	print args.counts_prefix
	if args.long:
		generate_long_table(all_fastas)
	elif args.counts_prefix !="None":
		generate_count_table(all_fastas,args.counts_prefix)
	else:
		generate_wide_table(all_fastas)

if __name__ == "__main__":
	sys.exit(main())

