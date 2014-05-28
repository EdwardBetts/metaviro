#! /usr/bin/env python
# encoding: utf-8
import itertools
import random
import math
# import prettytable
import sys
import collections
from Bio.Seq import Seq
from Bio.SeqUtils import CheckSum
import scipy
from Bio import SeqUtils
from Bio import SeqIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

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

Receive a set of multi-fasta files (or a root directory to recurse into), and perform stratified sampling of contigs with a specified length distribution. 
Stratification is based on the directory structure. 

TODO: A tree representing the location and length of all sequences is built, then length are totalled to build the samping amount.

Finally, the sampling is performed by 
Options
========

Example usage
==============
%(tool)s genome1.fa genome2.fa

Todo
======
"""%{"tool":sys.argv[0]}


def find_fasta_files(src_dir):
	all_fastas=[]
	for dirname, dirnames, filenames in os.walk(src_dir):
	 # print path to all filenames.
		for filename in filenames:
			if filename.endswith(".fasta") | filename.endswith(".fa"):
				all_fastas.append(os.path.join(dirname, filename))
	return all_fastas

def main(argv=None):
	parser=argparse.ArgumentParser(description="Compute kmer counts for all sequences in the input multi fasta file")
	# parser.add_argument('-k',dest="kmer_length",help="mer length",default=3,type=int)
	parser.add_argument('-r',dest="recursive",help="Perform recursive search for fasta files",action="store_true")


	parser.add_argument('-n',dest="n_contigs",help="Number of total contigs to sample",default=1e4,type=int)
	parser.add_argument('-l',dest="avg_length",help="Length of conting to sample",default=500,type=int)
	parser.add_argument("-s", "--stdev",dest="stdev", help="standard deviation", type=int, default=200)
	parser.add_argument("-p", dest="preview", help="Only print preview of number of contigs output", default=False,action="store_true")

	parser.add_argument('-I',dest="keep_IUPAC",help="Indicate whether to keep ambiguous site",default=False,action="store_true")
	parser.add_argument('-s',dest="simplify_ncbi_fasta_header",help="If set, headers of fasta file will only contain the GI number from the NCBI full description",default=False,action="store_true")
	parser.add_argument('-A',dest="append",help="If set, fasta sequences are appended to the output file",default=False,action="store_true")
	parser.add_argument('-o',dest="output_name",help="Name of output file",default="contigs.fasta",type=str)

	parser.add_argument('FASTAFILE',action='append',nargs="+",help='list of fasta files')
	args=parser.parse_args()


	FASTAFILE=args.FASTAFILE[0]
	all_fastas=set()
	if args.recursive:
		for f in FASTAFILE:
			all_fastas.update(find_fasta_files(f))
		logger.info("Recursively found %d fasta files"%(len(all_fastas)))
	else:
		all_fastas=set(FASTAFILE)

	if len(all_fastas)<1:
		logger.critical("No fasta files found, bailing out")
		sys.exit(1)
	multi_fasta_lengths=collections.defaultdict(int)
	sequence_lengths={} # Assuming description is unique
	fasta_to_sequences=collections.defaultdict(list)

	sequences={}

	for f in all_fastas: 
		for record in SeqIO.parse(f, "fasta", generic_dna):
			fasta_to_sequences[f].append(record.description)
			sequence_lengths[record.description]=len(record)
			multi_fasta_lengths[f]+=len(record)
			sequences[record.description]=record

	# compute the number of samples to take from each fasta 
	total_length=sum(multi_fasta_lengths.values())
	n_samples=dict([(k,math.ceil(float(v)/total_length*args.n_contigs)) for k,v in multi_fasta_lengths.items()])
	all_records=[]
	for fasta,n_samples in n_samples.items():
		# Get the length of these sequences
		these_sequence_length = dict([(k,v) for k,v in sequence_lengths.items() if k in fasta_to_sequences[fasta]])
		total_length=sum(these_sequence_length.values())

		# compute the number of samples per sequences
		n_samples_this_fasta=dict([(k,int(math.ceil(float(v)/total_length*n_samples))) for k,v in these_sequence_length.items()])

		logger.info("Fasta: %s,length:%d, N Sample:%d"%(fasta, multi_fasta_lengths[fasta],n_samples))
		if(args.preview:
			continue

		fasta_name=os.path.split(fasta)[-1]
		for k,v in n_samples_this_fasta.items():
			# logger.info("Fasta:%s, Seq: %s,length:%d, N Sample:%d"%(fasta,k, sequence_lengths[k],v))
			this_contig_length=int(random.normalvariate(args.mean, args.stdev))
			seq=str(sequences[k].seq)

			for i in range(0,v):
				# Sample start position,

				if len(seq)>=this_contig_length:
					start=random.randint(0,len(seq)-this_contig_length)
					end=start+this_contig_length
				else:
					start=0
					end=this_contig_length

				sub_seq=seq[start:start+this_contig_length]

                if args.reverse and bool(random.getrandbits(1)):
                    forward_sequence = Seq(sub_seq, generic_dna)
                    sub_seq=forward_sequence.reverse_complement()


				if (not args.keep_IUPAC) and (set(sub_seq)!=set(['A','C','G','T'])):
					# contains ambiguous 
					continue
				if args.simplify_ncbi_fasta_header:
					seq_id = sequences[k].id.split("|")[1]+"_sample_"+str(i)+"_"+fasta_name
					seq_id = seq_id.replace(".","_")
					seq_name=seq_id
					sequence_description=seq_id

				else:
					seq_id = sequences[k].id+"sample_"+str(i)
					seq_name=sequences[k].name+"sample_"+str(i)
					sequence_description=sequences[k].description

				record = SeqRecord(Seq(sub_seq,generic_dna),id=seq_id, name=seq_name,description=sequence_description)
				# record = SeqRecord(Seq(sub_seq,generic_dna))
				all_records.append(record)
	if args.append:
		output_handle = open(args.output_name, "wa")
	else:
		output_handle = open(args.output_name, "w")
	SeqIO.write(all_records, output_handle, "fasta")
	output_handle.close()





	# 		seq=str(record.seq)
	# 		fasta_keys=(f,record.description,str(len(seq)))
	# 		records_to_kmer[fasta_keys]=collections.defaultdict(int)
	# 		for i in range(0,len(seq)-args.kmer_length):
	# 			kmer=seq[i:i+args.kmer_length]
	# 			records_to_kmer[fasta_keys][kmer]+=1

	# print "\t".join(["path","sequence_description","sequence_length"]+all_kmers)
	# for k,kmer_values in records_to_kmer.items():
	# 	all_values = list(k)
	# 	all_values.extend(map(str,[kmer_values.get(x,0) for x in all_kmers]))
	# 	# print len(all_values)
	# 	print "\t".join(all_values)

if __name__ == "__main__":
	sys.exit(main())

