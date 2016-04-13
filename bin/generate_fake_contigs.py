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
total_sequences = []
total_ids = []
for record in SeqIO.parse(infile, "fasta"):
    total_sequences.append(record.seq)
    total_ids.append(record.id)
print len(total_sequences), "sequences in the input file"
infile.close()

# Checking boundaries
longest_sequence_length =  len(max(total_sequences, key=len))
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

# Creating total cumulative lengths list
total_cum_lengths = []
cumul = 0
for sequence in total_sequences:
    cumul = cumul + len(sequence)
    total_cum_lengths.append(cumul)
#print cum_lengths

contigs = []
idscontigs = []
used_ids = []
left_to_go = args.number

# Get sequences to resample from (so that a maximum of sequences are represented)
def get_sequences_left():
    if len(used_ids) == 0 or len(total_ids) == len(used_ids):
        return total_ids, total_sequences, total_cum_lengths
    else:
        left_ids = []
        left_sequences = []
        for i in range(len(total_ids)):
            if total_ids[i] not in used_ids:
                left_ids.append(total_ids[i])
                left_sequences.append(total_sequences[i])
        left_cum_lengths = []
        cumul = 0
        for sequence in left_sequences:
            cumul = cumul + len(sequence)
            left_cum_lengths.append(cumul)
        return left_ids, left_sequences, left_cum_lengths
            

# There is indeed a tiny chance that this will end up as an endless loop...
# If that ever happens, I'm willing to eat my damn keyboard.
while left_to_go > 0:
    print "New loop, contigs left to go =", left_to_go
    # Get sequences to work with:
    ids, sequences, cum_lengths = get_sequences_left()
#    ids, sequences, cum_lengths = total_ids, total_sequences, total_cum_lengths
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
                if ids[index] not in used_ids:
                    used_ids.append(ids[index])
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
                idscontigs.append(ids[index])
                left_to_go = left_to_go - 1
            else:
#                print "Doesn't fits, no sits :3"
                pass
        else:
#            print "Length out of bounds, skipping"
            pass

# How many of the original sequences have we kept in the resulting contigs?
print "{:.2%}".format(float(len(used_ids))/len(total_ids)) + " of the original sequences are present in the resampling."

# Saving results
print "writing to file..."
outputf=open(outfile, "w")
i = 1
for i in range(len(contigs)):
    outputf.write(">%d_%s_%s_len%d\n" % (i, os.path.splitext(os.path.basename(args.input_file))[0], idscontigs[i], len(contigs[i])) )
    outputf.write("%s\n" % contigs[i])
    i = i+1
outputf.close()

# Saving used sequences ids
print "writing ids to file..."
outputf=open("%s.used_ids" % outfile, "w")
for header in used_ids:
    outputf.write("%s\n" % header)
outputf.close()

#for i in range(len(ids)):
#    if ids[i] not in used_ids:
#        print sequences[i]

print "ta-da"
