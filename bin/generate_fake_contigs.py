# -*- coding: utf-8 -*-
"""
Generates random contigs in fasta format, of random lengths according to a normal distribution, picked from a a set of sequences in FASTA format.

@author: Louise-Amélie Schmitt
"""

import argparse, random, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser( description = "Generates fake contigs randomly selecting bits from a set of sequences from a multifasta file given as an argument. The sizes of the contigs are randomly chosen, following a normal distribution that can be tuned with the parameters." )
argparser.add_argument("-n", "--number", help="number of sequences", type=int, default=100000)
argparser.add_argument("-m", "--mean", help="mean length of random sequences", type=int, default=500)
argparser.add_argument("-s", "--stdev", help="standard deviation", type=int, default=200)
argparser.add_argument("-t", "--top", help="maximun length", type=int, default=4000)
argparser.add_argument("-b", "--bottom", help="minimum length", type=int, default=100)
argparser.add_argument("-r", "--reverse", help="create reverse complement contigs too (p=0.5)", action="store_true")
argparser.add_argument("-o", "--out", help="output file name", default=None)
argparser.add_argument("input_file", help="Multifasta file to generate the contigs from")
args = argparser.parse_args()

# Loading input file contents
infile = open(args.input_file, "rU")
sequences = []
for record in SeqIO.parse(infile, "fasta"):
    sequences.append(record.seq)
print len(sequences), "sequences in the input file"
infile.close()

# Checking boundaries
longest_sequence_length =  len(max(sequences, key=len))
if args.top > longest_sequence_length:
    maxlength = longest_sequence_length
else:
    maxlength = args.top
if maxlength < args.bottom:
    raise ValueError("max length is lower than min length - most likely, your input file doesn't contain sequences long enough for your min length.")

# Setting output file name
if args.out:
    outfile = args.out
else:
    outfile = os.path.splitext(args.input_file)[0] + "_contigs.fasta"

# Creating cumulative lengths list
cum_lengths = []
cumul = 0
for sequence in sequences:
    cumul = cumul + len(sequence)
    cum_lengths.append(cumul)
#print cum_lengths

contigs = []
left_to_go = args.number

# There is indeed a tiny chance that this will end up as an endless loop...
# If that ever happens, I'm willing to eat my damn keyboard.
while left_to_go > 0:
    print "New loop, contigs left to go =", left_to_go
    # First create the sample of positions we want to process in this loop turn
    sample = []
    for i in range(left_to_go): # Most efficient way I found to sample positions, not pretty but meh.
        sample.append(random.randrange(0, cum_lengths[len(cum_lengths) -1]))
    for value in sample:
        # Set the length once and for all instead of retrying if it doesn't fit
        # hopefully this will reduce the bias induced in the distribution
        length = int(random.normalvariate(args.mean, args.stdev))
        if length >= args.bottom and length < maxlength:
            # Then, check from which sequence this position belongs
            for index in range(len(cum_lengths)):
                if value < cum_lengths[index]:
                    break
            # If it fits I sits :3
            if value + length <= cum_lengths[index]:
                # Calculating the position within the sequence from the cumulative values
                if index == 0:
                    start = value
                    end = value + length + 1 # +1 'cause it cuts before the end index
                else:
                    start = value - cum_lengths[index - 1]
                    end = value - cum_lengths[index - 1] + length + 1
                # If we also want reverse complement contigs, let's make some half of the time
                if args.reverse and bool(random.getrandbits(1)):
                    forward_sequence = Seq(str(sequences[index][start:end]), generic_dna)
                    contigs.append(forward_sequence.reverse_complement())
                else:
                    contigs.append(sequences[index][start:end])
                left_to_go = left_to_go - 1
            else:
#                print "Doesn't fits, no sits :3"
                pass
        else:
#            print "Length out of bounds, skipping"
            pass
# Saving results
print "writing to file..."
outputf=open(outfile, "w")
i = 1
for contig in contigs:
    outputf.write(">%d_%s_len%d\n" % (i, os.path.splitext(os.path.basename(args.input_file))[0], len(contig)) )
    outputf.write("%s\n" % contig)
    i = i+1
outputf.close()

print "ta-da"