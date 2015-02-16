#! /usr/bin/env python
# encoding: utf-8
import pprint
pp = pprint.PrettyPrinter(indent=4)
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
%(tool)s  BLAST_DB_CONTENT.txt 


Reads a tab separated description of sequences from a blast db  and generate coordinates of randomly sampled sub-sequences.


Input format 
============
The input table can be generated with commands such as 
	blastdbcmd -entry all -db {db_name} -outfmt "%i\t%a\t%g\t%t\t%l" > {output}

Options
========
Number of contigs to samples etc. 

"""
# %{"tool":sys.argv[0]}

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

def sample_intervals(interval_lengths,n_contigs,BYSEQUENCE=False,AT_LEAST_ONE=False,length_avg=500,length_sd=200,min_length=200):
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

				this_contig_length=max(abs(min(int(math.ceil(random.normalvariate(length_avg,length_sd))),interval_lengths[j])),min_length)

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
				# for i in range(len(contig_coordinates)-4,len(contig_coordinates)):
				# 	assert contig_coordinates[i][1]<=contig_coordinates[i][2]
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



def parse_blast_content(path):
	all_contents=[]
	with open(path,"r") as f:
		for line in f:
			cols=line.split("\t")
			all_contents.append((cols[0],int(cols[-1])))

	return all_contents


def main(argv=None):
	parser=argparse.ArgumentParser(description="Generate sampling coordinates from a blast DB description.")
	parser.add_argument('-n',dest="n_contigs",help="Number of total contigs to sample",default=1000,type=int)
	parser.add_argument('-S',dest="by_sequence",help="Sampling of contigs is made for each sequence",default=False,action="store_true")
	parser.add_argument('-l',dest="avg_length",help="Length of contig to sample",default=500,type=int)
	parser.add_argument('-m',dest="min_length",help="Minimal length of contig to sample",default=200,type=int)
	parser.add_argument("-d", "--stdev",dest="stdev", help="standard deviation", type=int, default=200)
	parser.add_argument("-p", dest="preview", help="Only print preview of number of contigs output", default=False,action="store_true")
	parser.add_argument('-P',dest="pretty",help="Use (and require) prettytable for summary output",default=False,action="store_true")
	parser.add_argument('-u',dest="at_least_one",help="Enable subsampling: less than one contig per sequence is possible",default=True,action="store_false")
	parser.add_argument('-o',dest="output_name",help="Name of output file",default="sampled_contigs.txt",type=str)

	parser.add_argument('BLASTCONTENTFILE',action='append',nargs="+",help='TSV file with sequence ID (1st column) and sequence length (last column)')
	args=parser.parse_args()

	if args.pretty:
		import prettytable


	if len(args.BLASTCONTENTFILE[0])!=1:
		logger.critical("Not enough arguments, BLASTCONTENTFILE mandatory")
		return -1
	BLASTCONTENTFILE=args.BLASTCONTENTFILE[0][0]
	content=parse_blast_content(BLASTCONTENTFILE)
	logger.info("Parsed blastDB content for %d sequences totalling %d nt"%(len(content),sum([x[1] for x in content])))
	logger.info("First lines are %s"%(pp.pformat(content[0:5])))

	# Build sequence DB

	sequence_lengths=dict(content)
	sequence_length_keys=sequence_lengths.keys()
	sequence_length_values=sequence_lengths.values()

	total_length=sum(sequence_length_values)
	logger.info("Computed total sequence lengths")

	all_records=[]

	if args.pretty:
		table=prettytable.PrettyTable(["Fasta","length(kb)","N sequences","N sampled contig","Avg contig/seq","median","min", "max"])
		table.align["File"] = "l" 

	# these_sequence_length = dict([(k,v) for k,v in sequence_lengths.items() if k in fasta_to_sequences[fasta]])
	
	contig_coordinates=sample_intervals(sequence_length_values,args.n_contigs,BYSEQUENCE=args.by_sequence,AT_LEAST_ONE=args.at_least_one,length_avg=args.avg_length,length_sd=args.stdev,min_length=args.min_length)
	
	#index contig coordinates 
	contig_coordinates_by_key=collections.defaultdict(list)
	for c in contig_coordinates:
		contig_coordinates_by_key[c[0]].append(c)

	n_contigs_per_key=[len(v) for v in contig_coordinates_by_key.values()]

	avg_contig_per_sequence=scipy.average(n_contigs_per_key)
	med_contig_per_sequence=scipy.median(n_contigs_per_key)
	min_contig_per_sequence=min(n_contigs_per_key)
	max_contig_per_sequence=max(n_contigs_per_key)

	info=[BLASTCONTENTFILE, total_length/1000.0,len(content),len(contig_coordinates),avg_contig_per_sequence,med_contig_per_sequence,min_contig_per_sequence,max_contig_per_sequence]
	if args.pretty:
		table.add_row(info)
	else:
		logger.info("DB: %s,length:%d, N sequences: %d, N sampled contig:%d, Avg contig per seq:%f, med: %d, min seq:%d, max seq:%d "%(tuple(info)))
	

	if(args.preview):
		return 

	with open(args.output_name,"w") as f:
		for c in contig_coordinates:
			f.write("%s %d-%d\n"%(sequence_length_keys[c[0]],c[1],c[2]))

	return 
	fasta_name='/'.join(os.path.split(BLASTCONTENTFILE)[-2:])
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

