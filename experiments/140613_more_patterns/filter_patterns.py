## Goal is to reduce the space of patterns by removing too frequent or too rare k-mers 
## We do it in python since the file cannot be read in memory by R 

import os,sys
from collections import defaultdict


occurence_counts=defaultdict(set)

# First we count the number of lines 

# fh=open("long_kmers_m300_n4_w9_p40.csv")
# lines = 0
# buf_size = 1024 * 1024
# read_f = fh.read # loop optimization

# buf = read_f(buf_size)
# while buf:
# 	lines += buf.count('\n')
# 	buf = read_f(buf_size)

# print lines

fh=open("long_kmers_m300_n4_w9_p40.csv")



# Skip first line
fh.next()

for line in fh:
	fields=line.split("\t")
	path,species,sequence_description,sequence_length,GC,pattern,kmer,count=fields
	occurence_counts[kmer].add(int(species))
	if len(occurence_counts)>4000000:
		break
	if (len(occurence_counts) % 50000)==0:
		print len(occurence_counts)


# We tally 

most_frequent_kmers=sorted(occurence_counts.items(),key=lambda x:len(x[1]))
print len(most_frequent_kmers)
print most_frequent_kmers[-5:]
print "Kept:",len([x for x in occurence_counts if len(occurence_counts[x])>=5])
