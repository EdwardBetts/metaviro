
import sys
try:
   import cPickle as pickle
except:
   import pickle

# parameters: blastout.virless arch.fasta bact.fasta euk.fasta
print "parsing blast results..."
results = {}
current_contig = ''
prevline = ''
with open(sys.argv[1], 'r') as blastout:
	for line in blastout.readlines():
		if line[:6] != "Query=" and line[0] != ">":
			continue
		if prevline != '' and prevline[:6] == "Query=":
			current_contig = prevline[7:].strip()
			if line[0] == ">":
				results[current_contig] = 1
			else:
				results[current_contig] = 0
		prevline = line

if prevline[:6] == "Query=":
	results[prevline[7:].strip()] = 0

print "adding non viral headers..."
for fastafile in sys.argv[2:]:
    with open(fastafile, 'r') as fh:
        for line in fh.readlines():
            if line[0] == ">":
                results[line[1:].strip()] = -1

with open("filtered_contigs.bdat", "w") as fh:
	pickle.dump(results, fh)

print "done"