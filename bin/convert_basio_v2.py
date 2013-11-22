#! /usr/bin/env python
# encoding: utf-8
# Make the basio output into a table by counting the number of occurences of each segment 
import itertools
import sys
import collections
import sys
import os
import argparse

from Bio.Seq import Seq
from Bio.SeqUtils import CheckSum
from Bio import SeqUtils
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import scipy
from scipy.sparse import coo_matrix

import csv


import logging
if "logger" not in globals():
	logger = logging.getLogger('basio_converter')
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
%(tool)s  basio_file_1.sgm original_sequence.fa

Options
========

Example usage
==============
%(tool)s  basio_file_1.sgm original_sequence.fa

Todo
======
* Additional output format: sparse matrix, turtled sequence, splitted sequences
* Allow for fasta id used as filters? 
* Allow for a globl pre-existing Sequence -> ID mapping instead of building a new one
* Allow specification of a domain recognizer in re format : E.g. (viruses|virus)|(archea) etc.

"""%{"tool":sys.argv[0]}


fasta_index={}
basio_breakpoints=collections.defaultdict(list)
all_coordinates=None

# Map tables
all_segments_to_id={}
all_id_to_segments={}
all_fasta_keys_to_id={}
all_id_to_fasta_keys={}
all_fasta_keys_to_sampling_file={}

def index_fasta_files(fasta_files):
	global fasta_index,all_fasta_keys_to_id,all_id_to_fasta_keys
	fasta_keys_in_order=[]
	for a_file in fasta_files:
		logger.debug("Processing fasta file %s"%(a_file))
		for record in SeqIO.parse(a_file, "fasta", generic_dna):
			if record.id in fasta_index:
				if fasta_index[record.id].seq != record.seq:
					logger.debug("Duplicate fasta ID with different sequences:",record.id)
					sys.exit(1)
			fasta_index[record.id]=record
			fasta_keys_in_order.append(record.id)
			all_fasta_keys_to_sampling_file[record.id]=a_file
	all_id_to_fasta_keys=dict(enumerate(fasta_keys_in_order))
	all_fasta_keys_to_id=dict([reversed(i) for i in all_id_to_fasta_keys.items()])


def index_breakpoints(segment_files):
	global fasta_index,basio_breakpoints
	for a_file in segment_files:
		logger.debug("Processing segment file %s"%(a_file))
		last_fast_idx=None
		last_full_fasta_path=None
		for content in open(a_file):
			content=content.strip()
			if content[0]==">": # Start a new sequence
				last_fast_idx=content.strip(">")
				assert last_fast_idx in fasta_index,"FASTA ID %s in segment file but not in FASTA files"%(last_fast_idx)
				continue
			else:
				elements=content.split()
				if elements[0] in ["AL","FN","NR","VD","SC"]:
					continue
				elif elements[0] == 'SC' : #Get the path of the original fasta
					this_path=elements.split()[1]
					this_path_components=os.path.split(this_path)
					if this_path_components[0].startwiths("/tmp/"):
						this_path_components[0]=this_path[4:]

				else:
					assert elements[0]=="BR","syntax error on line:%s"%(content)
					elements[0]=last_fast_idx
					basio_breakpoints[last_fast_idx].append(int(elements[1]))


def build_sparse_matrix():
	global fasta_index,basio_breakpoints
	global all_coordinates,all_segments_to_id,all_id_to_segments,all_values
	logger.info("Indexing for sparse matrix")
	
	next_sequence_uid=0
	rows,cols,values=[],[],[]
	next_insertion=0

	#Preallocation
	all_values=scipy.zeros((30000,3))

	for fasta_key,breakpoints in basio_breakpoints.items():
		this_contig_bow=collections.defaultdict(int)
		# this_contig_sequences=[]
		for i in range(0,len(breakpoints)-1):
			start,end=breakpoints[i],breakpoints[i+1]
			this_contig_bow[str(fasta_index[fasta_key].seq[start:end])]+=1

		# Add to  the sparse matrix representation
		for sequence,multiplicity in this_contig_bow.items():
			if sequence not in all_segments_to_id:
				all_segments_to_id[sequence]=next_sequence_uid
				next_sequence_uid+=1

			if next_insertion>=all_values.shape[0]:
				#Augment array
				all_values=scipy.vstack([all_values,scipy.zeros((30000,3))])
			all_values[next_insertion,0]=all_fasta_keys_to_id[fasta_key]
			all_values[next_insertion,1]=all_segments_to_id[sequence]
			all_values[next_insertion,2]=multiplicity
			next_insertion+=1
			# rows.append(all_fasta_keys_to_id[fasta_key])
			# cols.append(all_segments_to_id[sequence])
			# values.append(multiplicity)

	# trim down
	all_values=all_values[:next_insertion,]

	# logger.info("Building the sparse matrix:%d items"%(len(values)))
	logger.info("Building the sparse matrix:%d items"%(all_values.shape[0]))
	# all_coordinates=coo_matrix((values, (rows,cols)))
	all_coordinates=coo_matrix((all_values[:,2], (all_values[:,0],all_values[:,1])))
	logger.info("Built the sparse matrix:%d x %d items"%(all_coordinates.shape[0],all_coordinates.shape[1]))

	all_id_to_segments=dict([reversed(i) for i in all_segments_to_id.items()])


def print_splitted_sequence(print_restriction):
	global fasta_index,basio_breakpoints,args
	# global all_coordinates

	for fasta_key,breakpoints in basio_breakpoints.items():
		print fasta_key,"\t",
		for i in range(0,len(breakpoints)-1):
			start,end=breakpoints[i],breakpoints[i+1]
			if print_restriction and ((end-start)!=print_restriction):
				continue
			print fasta_index[fasta_key].seq[start:end],
		print ""

	



def save_output_files(output_prefix):
	segment_mult_file=output_prefix+"_mvdist.csv"
	logger.info("Writing sparse matrix")

	with open(segment_mult_file, 'wb') as csvfile:
		csvwriter = csv.writer(csvfile)
		all_coordinates_csr=all_coordinates.tocsr()
		csvwriter.writerow(['path','sequence','segid','count'])
		for contig_id in range(all_coordinates_csr.shape[0]):
			fasta_key=all_id_to_fasta_keys[contig_id]
			fasta_path=all_fasta_keys_to_sampling_file[fasta_key]

			seg_id=all_coordinates_csr[contig_id,:].nonzero()[1]
			counts=all_coordinates_csr[contig_id,:].data
			for component in range(len(seg_id)):
				# csvwriter.writerow(this_contig_meta+[seg_id[component],counts[component]])
				csvwriter.writerow((fasta_path,fasta_key,seg_id[component],counts[component]))

	segment_id_file=output_prefix+"_mvseg.csv"
	logger.info("Writing id to segment mapping")

	with open(segment_id_file, 'wb') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(['segid','segment'])
		csvwriter.writerows(all_id_to_segments.items())



def process(segment_files,fasta_files,output_prefix):
	global args
	index_fasta_files(fasta_files)
	index_breakpoints(segment_files)
	if args.print_only:
		print_splitted_sequence(args.print_restriction)
	else:
		build_sparse_matrix()
		save_output_files(output_prefix)
	pass



def main(argv=None):
	global args
	parser=argparse.ArgumentParser(description="""Make the basio output into a table by counting the number of occurences of each segment """)
	parser.add_argument("-o", "--out", dest="out", help ="Prefix of output file: will generate files OUT.mvdist OUT.mvseg")
	parser.add_argument("-p", "--print", dest="print_only", help ="Don't write files but instead print segmented sequence",action="store_true")
	parser.add_argument("-s", "--print_restriction", dest="print_restriction", help ="Only print parts of segments of this size",type=int)
	parser.add_argument('INPUTFILES',action='append',nargs="+",help='list of fasta (.fa | .fasta) and segment files (.sgm), in any order')
	args=parser.parse_args()

	INPUTFILES=args.INPUTFILES[0]
	SEGMENTFILES=[x for x  in INPUTFILES if x.endswith(".sgm")]
	FASTAFILES=[x for x in INPUTFILES if x.endswith(".fa")|x.endswith(".fasta")]
	assert len(FASTAFILES)!=0,"No fasta file recognized in input %r"%(INPUTFILES)
	assert len(SEGMENTFILES)!=0,"No segment file recognized in input %r"%(INPUTFILES)
	if not args.print_only and not args.out:
		print "Either -p or -o must be specified"
		sys.exit(1)
	process(SEGMENTFILES,FASTAFILES,args.out)

if __name__ == "__main__":
	sys.exit(main())

