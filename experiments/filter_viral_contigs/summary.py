# -*- coding: utf-8 -*-

'''
Give for one distance measure and one query, all the files for different k values
'''
from collections import defaultdict
import sys, string, os

data = defaultdict(list)
kvals = []

files = sys.argv[1:]


for i in range(20):
	for f in files:
		if string.find(f, "_k%d" % i) >= 0:
			firstline = None
			lastline = None
			with open(f, "r") as fh:
				for line in fh.readlines():
					if not firstline and string.find(line, "contig_id") == 0:
						firstline = line.split()
					if string.find(line, "cumulative_max") == 0:
						lastline = line.split()
			for ref, counts in zip(firstline[1:len(firstline)-1], lastline[1:len(lastline)-1]):
				ref = ref[:string.find(ref, "_k%d" % i)]
				data[ref].append(counts)
			kvals.append(str(i))



with open("%s_kmer_%s.csv" %(files[0][:string.find(files[0], "_k")], os.path.splitext(files[0])[1][1:]), "w") as csvfh:
	csvfh.write("label %s\n" % " ".join(kvals))
	sortedkeys = data.keys()
	sortedkeys.sort()
	for key in sortedkeys:
		csvfh.write("%s %s\n" % (key, " ".join(data[key])) )
