# -*- coding: utf-8 -*-

# Generates random sequences in fasta format, of random lengths according to a
# normal distribution. 

import argparse, random

argparser = argparse.ArgumentParser()
argparser.add_argument("-n", help="number of sequences", type=int, default=100000)
argparser.add_argument("-m", help="mean length of random sequences", type=int, default=500)
argparser.add_argument("-s", help="standard deviation", type=int, default=200)
argparser.add_argument("-t", help="maximun length", type=int, default=4000)
argparser.add_argument("-b", help="minimum length", type=int, default=100)
argparser.add_argument("-o", help="output file name", default="randseq.fa")
args = argparser.parse_args()

outputf=open(args.o, "w")
for i in range(args.n):
	while True:
		length = int(random.normalvariate(args.m, args.s))
		if length >= 100:
			seq=''.join(random.choice("ATGC") for i in xrange(length))
			outputf.write(">randomsequence_len%d\n" % length )
			outputf.write("%s\n" % seq)
			break
outputf.close()
