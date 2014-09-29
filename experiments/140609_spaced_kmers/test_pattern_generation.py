#! /usr/bin/env python
# encoding: utf-8
import random
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

# from kmerize import logger

### Generate some random patterns 

## User defined parameters 
max_1 = 8
n_patterns = 30
output_file="test_w8_30pat.csv"


all_patterns=[]
for i in xrange(n_patterns):
	number_zeros=random.choice(xrange(0,max_1+1))
	pattern_size=max_1+number_zeros
	# make a permutation
	this_pattern=random.sample([1]*max_1 + [0] * number_zeros,k=pattern_size)
	all_patterns.append(tuple(this_pattern))
all_patterns=set(all_patterns)
print ["".join(map(str,x)) for x in all_patterns]


## Set inputs 

all_fastas=[METAVIROOT+"experiments/140530_long_kmer_analysis/long_kmers_b2.fa"]

"""
***TODO***
Add column ID line to output_handle
Split the species etc. during this process
"""
output_handle=open(output_file,"wa")
print >>output_handle, "\t".join(["path","sequence_description","sequence_length","GC","pattern","kmer","count"])

for f in all_fastas: 
	logger.debug("Processing file %s"%(f))
	for record in SeqIO.parse(f, "fasta", generic_dna):

		logger.debug("Processing sequence %s"%(record.description))

		seq=str(record.seq)
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
			for kmer, count in kmer_count.items():
				print >>output_handle,"\t".join(fasta_keys+["".join(map(str,pattern)),kmer,str(count)])





## Parse the sequences and compute the mean distance between sequences of the same taxa 