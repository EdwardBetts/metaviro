#! /usr/bin/env python
# encoding: utf-8


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
import logging

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('mix')


help="""Usage 
%(tool)s  FASTAFILE FASTAFILE 

Will add all sequences in input FASTAFILE to target directory as multiple unifasta files

Options
========

Example usage
==============
%(tool)s genome1.fa genome2.fa


Todo
======
* Check sequences already present in the directory
"""%{"tool":sys.argv[1]}
		
def main(argv=None):
	parser=argparse.ArgumentParser(description="Split all sequences from the input fasta to single file unifasta")
	parser.add_argument("-o", "--out", dest="out", help ="Target directory of output files: Created if needed",required=True)
	parser.add_argument("-p", "--prefix", dest="prefix", help ="Prefix used to name the output file",default="sequence")
	parser.add_argument('FASTAFILE',action='append',nargs="+",help='List of fasta files to keep. Use "*" to keep them all')
	args=parser.parse_args()
	FASTAFILE=args.FASTAFILE[0]
	last_idx=0

	try:
		os.mkdir(args.out)
	except OSError:
		logger.info("output directory %s already exists"%(args.out))
		pass

	for f in FASTAFILE: 
		logger.info("processing %s"%(f))
		for record in SeqIO.parse(f, "fasta", generic_dna):
			tgt_file = args.out+"/"+args.prefix+"_"+str(last_idx)+".fasta"
			output_handle = open(tgt_file, "w")
			SeqIO.write([record], output_handle, "fasta")
			output_handle.close()
			last_idx+=1
if __name__ == "__main__":
	sys.exit(main())
