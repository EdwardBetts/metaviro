#! /usr/bin/env python
# encoding: utf-8

"""
Script to generate and k-merize nt sequences (FASTA) with spaced k-mers
Given a max number of non-starred k-mers, generate at most K random patterns, count their occurences per sequences and save results in a pattern specific file
"""
import random
import pprint
import os,sys
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


pp = pprint.PrettyPrinter(indent=4)



## Logger

import logging
if "logger" not in globals():
	logger = logging.getLogger('spaced_kmer')
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

METAVIROOT=os.path.expanduser("~/metaviro/")
sys.path.append(METAVIROOT+"bin/")

## Set inputs 

all_fastas=["long_kmers_m300_n4.fa"]
split_delim=None #Whether we should split the FASTA id or not 
split_delim="_" #Whether we should split the FASTA id or not 
# output_file="long_kmers_m300_n4_w9_p40.csv"


# Generate list of existing files 
existing_files=os.listdir(".")




### Generate some random patterns 

## User defined parameters 
max_1 = 8
n_patterns = 30



all_patterns=[]
for i in xrange(n_patterns):
	number_zeros=random.choice(xrange(0,max_1+1))
	pattern_size=max_1+number_zeros
	# make a permutation
	this_pattern=random.sample([1]*max_1 + [0] * number_zeros,k=pattern_size)
	# We strip leftmost and rightmost "0" 
	this_pattern= map(int,list("".join(map(str,this_pattern)).strip("0")))
	all_patterns.append(tuple(this_pattern))
all_patterns=set(all_patterns)
logger.info("Will k-merize with patterns:\n%s",pp.pformat(["".join(map(str,x)) for x in all_patterns]))


## Build output file name 

output_files=dict()
patterns_to_skip=set()
for p in all_patterns:
	p_str="".join(map(str,p))
	for f in all_fastas:
		output_file=f+"_"+p_str+".csv"
		if output_file in existing_files:
			logger.info("%s already existing for pattern %s and input %s"%(output_file,p_str,f))
			patterns_to_skip.add(p)
			continue
		output_files[(f,p)]=open(output_file,"w")
		logger.info("Creating file for pattern %s and input %s:%s"%(p_str,f,output_file))


all_patterns=all_patterns.difference(patterns_to_skip)

if len(all_patterns)==0:
	logger.info("No more available patterns, bailing out")
	sys.exit(0)

"""
***TODO***
Add column ID line to output_handle @DONE 
Split the species etc. during this process
"""
# output_handle=open(output_file,"wa")

for a_file in output_files.values():
	print >>a_file, "\t".join(["path","species","sequence_description","sequence_length","GC","pattern","kmer","count"])

for f in all_fastas: 
	logger.debug("Processing file %s"%(f))
	for record in SeqIO.parse(f, "fasta", generic_dna):

		logger.debug("Processing sequence %s"%(record.description))

		seq=str(record.seq)
		if split_delim:
			fasta_keys=[f]+[record.description.split(split_delim)[0]]+[record.description,str(len(seq))]
		else:
			fasta_keys=[f,record.description,str(len(seq))]
		# Add additional features to the sequence
		fasta_keys.append(str(SeqUtils.GC(record.seq)))

		# 
		# fasta_keys=tuple(fasta_keys)
		# kmer_count[fasta_keys]=collections.defaultdict(int)
		for pattern in all_patterns:
			kmer_count=collections.defaultdict(int)
			pattern_length=len(pattern)
			for i in range(0,len(seq)-pattern_length):
				kmer=seq[i:i+pattern_length]
				kmer=list(kmer)
				for i in range(0,pattern_length):
					if pattern[i]==0:
						kmer[i]="*"
				kmer="".join(kmer)
				kmer_count[kmer]+=1
			output_handle=output_files[f,pattern]
			for kmer, count in kmer_count.items():
				print >>output_handle,"\t".join(fasta_keys+["".join(map(str,pattern)),kmer,str(count)])



for a_file in output_files.values():
	a_file.close()



## Parse the sequences and compute the mean distance between sequences of the same taxa 