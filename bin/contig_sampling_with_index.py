#! /usr/bin/env python
# encoding: utf-8
import itertools
import random
import math
import numpy as np

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



def class_balanced_sampling(items,key_idx,n_by_class):
	sample_idx_by_key=collections.defaultdict(set)
	for i in range(len(items)):
		sample_idx_by_key[items[i][0]].add(i)

	selected_indices=set()
	for idx,sample_indices in sample_idx_by_key.items():
		selected_indices.update(random.sample(sample_indices,min(n_by_class,len(sample_indices))))
	final_items=[items[i] for i in selected_indices]
	rejected_items=[items[i] for i in range(len(items)) if i not in selected_indices]
	return final_items,rejected_items

def sample_intervals(interval_lengths,n_contigs,BYSEQUENCE=False,AT_LEAST_ONE=False,length_avg=500,length_sd=200):
	tot_length=sum(interval_lengths)
	interval_lengths_cs=np.hstack(([0],np.cumsum(interval_lengths)))
	sample_points=[]
	if AT_LEAST_ONE:
		for i in range(len(interval_lengths_cs)-1):
			sample_points.append(random.randint(interval_lengths_cs[i],interval_lengths_cs[i+1]))


	if BYSEQUENCE:
		for i in range(len(interval_lengths_cs)-1):
			sample_points.extend(random.sample(xrange(interval_lengths_cs[i],interval_lengths_cs[i+1]),n_contigs))
	else:
		sample_points.extend(random.sample(xrange(0,tot_length),n_contigs))

	sample_points.sort()

	contig_coordinates=[]

	last_start=0
	sequence_assigment=[-1]*len(sample_points)

	# For each sample_point, find the corresponding interval, then generate three sample possibilities 
	# where the sample point can be either start, end, or middle

	for i_idx in xrange(len(sample_points)):
		for j in xrange(last_start,len(interval_lengths_cs)-1):
			if ((sample_points[i_idx] >= interval_lengths_cs[j]) and (sample_points[i_idx] < interval_lengths_cs[j+1])):
				sequence_assigment[i_idx]=j
				# Generate length, then determine if sample_points is the start or the end of the contig 

				this_contig_length=min(int(math.ceil(random.normalvariate(length_avg,length_sd))),interval_lengths[j])

				# print "with length",this_contig_length,
				# possible_samplings=[]

				if sample_points[i_idx]+this_contig_length <= interval_lengths_cs[j+1]: # can be the start 
					contig_coordinates.append((j,sample_points[i_idx]-interval_lengths_cs[j],sample_points[i_idx]+this_contig_length-interval_lengths_cs[j]))

				if sample_points[i_idx]-this_contig_length >= interval_lengths_cs[j]: # can be the end
					contig_coordinates.append((j,sample_points[i_idx]-this_contig_length-interval_lengths_cs[j],sample_points[i_idx]-interval_lengths_cs[j]))

				if this_contig_length>interval_lengths[j]/2:
					# can be the middle 
					left_end=sample_points[i_idx]-this_contig_length/2-interval_lengths_cs[j]
					right_end=sample_points[i_idx]-interval_lengths_cs[j]+this_contig_length/2
					if left_end<=0: 
						left_end=0
						right_end=this_contig_length
					elif right_end>=interval_lengths[j]:
						right_end=interval_lengths[j]
						left_end=interval_lengths[j]-this_contig_length
					contig_coordinates.append((j,left_end,right_end))

				last_start=j
				break


	# reduce the contig_coordinates to match the requested number of contigs

	if AT_LEAST_ONE:
		guaranteed_contigs,remaining=class_balanced_sampling(contig_coordinates,0,1)
		contig_coordinates=remaining
	else:
		guaranteed_contigs=[]

	if BYSEQUENCE:
		# we sample at most n_contigs for each sequence 
		contig_coordinates,foo = class_balanced_sampling(contig_coordinates,0,n_contigs)
	else:
		contig_coordinates=random.sample(contig_coordinates,min(n_contigs,len(contig_coordinates)))

	contig_coordinates=set(contig_coordinates)
	contig_coordinates.update(guaranteed_contigs)
	contig_coordinates=list(contig_coordinates)

	# for idx,start,end in contig_coordinates:
	# 	assert(start<=interval_lengths[idx])
	# 	assert(end<=interval_lengths[idx])
	# 	assert(start>=0)


	# Sort 
	# contig_coordinates.sort(key=lambda x: x[0:2])

	return contig_coordinates


def main(argv=None):
	parser=argparse.ArgumentParser(description="Compute kmer counts for all sequences in the input multi fasta file")
	# parser.add_argument('-k',dest="kmer_length",help="mer length",default=3,type=int)
	parser.add_argument('-r',dest="recursive",help="Perform recursive search for fasta files",action="store_true")


	parser.add_argument('-n',dest="n_contigs",help="Number of total contigs to sample",default=1000,type=int)
	parser.add_argument('-S',dest="by_sequence",help="Sampling of contigs is made for each sequence",default=False,action="store_true")
	parser.add_argument('-l',dest="avg_length",help="Length of conting to sample",default=500,type=int)
	parser.add_argument('-m',dest="min_length",help="Minimal length of conting to sample",default=200,type=int)
	parser.add_argument("-d", "--stdev",dest="stdev", help="standard deviation", type=int, default=200)
	parser.add_argument("-p", dest="preview", help="Only print preview of number of contigs output", default=False,action="store_true")
	parser.add_argument("-R", dest="reverse", help="Randomly perform reverse complement", default=False,action="store_true")

	parser.add_argument('-I',dest="keep_IUPAC",help="Indicate whether to keep ambiguous site",default=False,action="store_true")
	parser.add_argument('-i',dest="index_path",help="MANDATORY path to an SQLite file index for out of memory sampling. If absent, will be generated; if present will be loaded back. Will trigger an error if an existing index do not match the provided input fasta files.",default=None,type=str)
	parser.add_argument('-s',dest="simplify_ncbi_fasta_header",help="If set, headers of fasta file will only contain the GI number from the NCBI full description",default=False,action="store_true")
	parser.add_argument('-A',dest="append",help="If set, fasta sequences are appended to the output file",default=False,action="store_true")
	parser.add_argument('-o',dest="output_name",help="Name of output file",default="contigs.fasta",type=str)
	parser.add_argument('-P',dest="pretty",help="Use (and require) prettytable for summary output",default=False,action="store_true")
	parser.add_argument('-k',dest="key",help="Indicate string key inserted in sampled contig headers",default="")
	parser.add_argument('-M',dest="max_sequences",help="Maximal number of sequences to consider,-1: all",default=-1,type=int)
	parser.add_argument('-F',dest="max_input_fastas",help="Maximal number of fastas files to process; -1: all",default=-1,type=int)
	parser.add_argument('-u',dest="at_least_one",help="Enable subsampling: less than one contig per sequence is possible",default=True,action="store_false")

	parser.add_argument('FASTAFILE',action='append',nargs="+",help='list of fasta files')
	args=parser.parse_args()

	if args.pretty:
		import prettytable


	FASTAFILE=args.FASTAFILE[0]
	all_fastas=set()
	if args.recursive:
		for f in FASTAFILE:
			all_fastas.update(find_fasta_files(f))
	else:
		all_fastas=set(FASTAFILE)

	logger.info("Found %d fasta files"%(len(all_fastas)))
	if len(all_fastas)<1 and not args.index_path:
		logger.critical("No fasta files found, no index provided, bailing out")
		sys.exit(1)
	if args.max_input_fastas != -1:
		logger.info("Downsampling input list of FASTA files down to %d"%(args.max_input_fastas))
		all_fastas=random.sample(list(all_fastas), k=args.max_input_fastas)
	all_fastas=list(all_fastas)


	# load the index if provided
	sequence_index=None
	if args.index_path:
		if len(all_fastas)>0:
			sequence_index=SeqIO.index_db(args.index_path,all_fastas,"fasta")
			logger.info("Built DB index %s"%(args.index_path))
		else:
			sequence_index=SeqIO.index_db(args.index_path)
			logger.info("Reusing DB index %s"%(args.index_path))

	if sequence_index==None:
		logger.critical("No index provided, bailing out")
		sys.exit(1)

	# Build sequence DB

	sequence_lengths=dict([(x,len(sequence_index[x])) for x in sequence_index])
	sequence_length_keys=sequence_lengths.keys()
	sequence_length_values=sequence_lengths.values()

	total_length=sum(sequence_length_values)
	logger.info("Computed total sequence lengths")

	all_records=[]

	if args.pretty:
		table=prettytable.PrettyTable(["Fasta","length(kb)","N sequences","N sampled contig","Avg contig/seq","median","min", "max"])
		table.align["File"] = "l" 

	# these_sequence_length = dict([(k,v) for k,v in sequence_lengths.items() if k in fasta_to_sequences[fasta]])
	
	contig_coordinates=sample_intervals(sequence_length_values,args.n_contigs,BYSEQUENCE=args.by_sequence,AT_LEAST_ONE=args.at_least_one,length_avg=args.avg_length,length_sd=args.stdev)
	
	#index contig coordinates 
	contig_coordinates_by_key=collections.defaultdict(list)
	for c in contig_coordinates:
		contig_coordinates_by_key[c[0]].append(c)

	n_contigs_per_key=[len(v) for v in contig_coordinates_by_key.values()]

	avg_contig_per_sequence=scipy.average(n_contigs_per_key)
	med_contig_per_sequence=scipy.median(n_contigs_per_key)
	min_contig_per_sequence=min(n_contigs_per_key)
	max_contig_per_sequence=max(n_contigs_per_key)

	info=[args.index_path, total_length/1000.0,len(sequence_index),len(contig_coordinates),avg_contig_per_sequence,med_contig_per_sequence,min_contig_per_sequence,max_contig_per_sequence]
	if args.pretty:
		table.add_row(info)
	else:
		logger.info("Fasta: %s,length:%d, N sequences: %d, N sampled contig:%d, Avg contig per seq:%f, med: %d, min seq:%d, max seq:%d "%(tuple(info)))
	

	if(args.preview):
		return 

	fasta_name='/'.join(os.path.split(args.index_path)[-2:])
	last_seq_idx=None
	sample_id=0
	generated_samples=0

	for seq_idx,start,end in contig_coordinates:
		if seq_idx!=last_seq_idx:
			sample_id=0
			fasta_key=sequence_length_keys[seq_idx]
			current_seq=str(sequence_index[fasta_key].seq)

			if args.simplify_ncbi_fasta_header:
				seq_id = sequence_index[fasta_key].id.split("|")[1]+"_"+args.key+"_sample_"+str(sample_id)+"_"+fasta_name
				seq_id = seq_id.replace(".","_")
				seq_name=seq_id
				sequence_description=seq_id
			else:
				seq_id = sequence_index[fasta_key].id+args.key+"_sample_"+str(sample_id)
				seq_name=sequence_index[fasta_key].name+args.key+"_sample_"+str(sample_id)
				sequence_description=sequence_index[fasta_key].description
		generated_samples+=1
		if (generated_samples % 500)==0: 
			logger.info("Generated %d/%d sequences"%(generated_samples,len(contig_coordinates)))

		sample_id+=1
		sub_seq=current_seq[start:end]
		if (len(sub_seq)<= args.min_length):
			continue

		if(len(sub_seq)>=10000):
			assert False

		if (not args.keep_IUPAC) and (set(sub_seq)!=set(['A','C','G','T'])):
			# contains ambiguous 
			continue

		record = SeqRecord(Seq(sub_seq,generic_dna),id=seq_id, name=seq_name,description=sequence_description)

		if args.reverse and bool(random.getrandbits(1)):
			recordRC=record.reverse_complement()
			recordRC.id=record.id+"_rev"
			# recordRC.name=record.name
			recordRC.description=""
			record=recordRC
		# record = SeqRecord(Seq(sub_seq,generic_dna))
		all_records.append(record)


	if args.pretty:
		print table.get_string()
	if args.append:
		logger.info("Appending to file")
		output_handle = open(args.output_name, "a")
	else:
		output_handle = open(args.output_name, "w")

	SeqIO.write(all_records, output_handle, "fasta")
	output_handle.close()


if __name__ == "__main__":
	sys.exit(main())

