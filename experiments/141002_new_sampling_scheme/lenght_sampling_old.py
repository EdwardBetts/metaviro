import math
import random
import scipy 
import numpy as np


# ll=[12,28,30,15,15]
NSEQUENCES=20
MAXSEQLENGTH=90
MINSEQLENGTH=20
n_contigs=15
length_avg=25
length_sd=2
MINCONTIGLENGTH=10
MAXCONTIGLENGTH=MAXSEQLENGTH

SEQUENCE_CHARS=list("ABCDEFGHIJKLMNOPQRSTUVXYZabcdefghijklmnopqrstuvxyz")


ll=np.random.randint(MINSEQLENGTH,MAXSEQLENGTH,NSEQUENCES)
coverages=[]
for i in range(len(ll)):
	coverages.append([0]*ll[i])

tot_length=sum(ll)

ll_cs=np.hstack(([0],np.cumsum(ll)))


n_samples={0:math.ceil(tot_length/tot_length*n_contigs)}

these_sequence_length = dict(zip(range(NSEQUENCES),ll))
total_length=sum(these_sequence_length.values())

n_samples_this_fasta=dict([(k,int(math.floor(float(v)/total_length*n_samples[0]))) for k,v in these_sequence_length.items()])
# In the case where the number of requested contig is less than the number of input sequence, pick at random at most n seq
n_total_samples=sum(n_samples_this_fasta.values())
if n_total_samples<n_contigs:
	selected_sequences=random.sample(n_samples_this_fasta.keys(),k=n_contigs)
	for k in selected_sequences:
		n_samples_this_fasta[k]=1

contig_to_sample=[]
for k,v in n_samples_this_fasta.items():
	seq="X"*these_sequence_length[k]
	for i in range(0,v):				# Sample start position,
		this_contig_length=min(int(random.gauss(length_avg, length_sd)),len(seq))

		if len(seq)>=this_contig_length:
			start=random.randint(0,len(seq)-this_contig_length)
			end=start+this_contig_length
		else:
			start=0
			end=this_contig_length

		sub_seq=seq[start:start+this_contig_length]
		if (len(sub_seq)<= MINCONTIGLENGTH):
			continue
		if(len(sub_seq)>=MAXCONTIGLENGTH):
			assert False

		# Randomly skip some ambiguous sites 
		if random.random()>=0.95:
			# contains ambiguous 
			continue
		contig_to_sample.append((k,start,end))


for idx,start,end in contig_to_sample:
	assert(start<=ll[idx])
	assert(end<=ll[idx])
	assert(start>=0)


# Sort 
contig_to_sample.sort(key=lambda x: x[0:2])

# update coverage 
for c_idx,start,end  in contig_to_sample:
	for i in range(start,end):
		coverages[c_idx][i]+=1


print "--"*12
# print coverage 
for c in coverages:
	print scipy.mean(c),c


# Print contigs 

last_cont_idx=None
for cont_idx,start,end in contig_to_sample:
	if last_cont_idx!=cont_idx:
		print "\n\n"
		print SEQUENCE_CHARS[cont_idx]*ll[cont_idx],len([x for x in sequence_assigment if x==cont_idx]),"bp"
		last_cont_idx=cont_idx
	print " "*(start)+"."*(end-start)+" "*(ll[cont_idx]-end)
# Make statistics 
