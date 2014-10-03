import math
import random
import scipy 
import collections
import numpy as np


# ll=[12,28,30,15,15]
NSEQUENCES=15
MAXSEQLENGTH=90
MINSEQLENGTH=20
n_contigs=5
length_avg=5
length_sd=2
BYSEQUENCE=False
AT_LEAST_ONE=False


SEQUENCE_CHARS=list("ABCDEFGHIJKLMNOPQRSTUVXYZabcdefghijklmnopqrstuvxyz")

# opt: Generate all length at once ? 
# opt: mode at least one contig per sequence
# opt: mode minimal/maximal length boundaries 
# opt: account for any distribution of contig length 
# opt: by sequence 



def class_balanced_sampling(items,key_idx,n_by_class):
	sample_idx_by_key=collections.defaultdict(set)
	for i in range(len(items)):
		sample_idx_by_key[items[i][0]].add(i)

	selected_indices=set()
	for idx,sample_indices in sample_idx_by_key.items():
		selected_indices.update(random.sample(sample_indices,min(n_by_class,len(sample_indices))))
	final_items=[items[i] for i in selected_indices]
	rejected_items=[items[i] for i in range(len(items)) if i not in selected_indices]
	return final_items,rejected_items




ll=np.random.randint(MINSEQLENGTH,MAXSEQLENGTH,NSEQUENCES)
coverages=[]
for i in range(len(ll)):
	coverages.append([0]*ll[i])
tot_length=sum(ll)
ll_cs=np.hstack(([0],np.cumsum(ll)))




sample_points=[]
if AT_LEAST_ONE:
	for i in range(len(ll_cs)-1):
		sample_points.append(random.randint(ll_cs[i],ll_cs[i+1]))


if BYSEQUENCE:
	for i in range(len(ll_cs)-1):
		sample_points.extend(random.sample(xrange(ll_cs[i],ll_cs[i+1]),n_contigs))
else:
	sample_points.extend(random.sample(xrange(0,tot_length),n_contigs))

sample_points.sort()

contig_to_sample=[]

last_start=0
sequence_assigment=[-1]*len(sample_points)


for i_idx in xrange(len(sample_points)):
	# print i_idx,":",sample_points[i_idx],
	for j in xrange(last_start,len(ll_cs)-1):
		if ((sample_points[i_idx] >= ll_cs[j]) and (sample_points[i_idx] < ll_cs[j+1])):
			sequence_assigment[i_idx]=j
			# print "Found in [",ll_cs[j],ll_cs[j+1],"](",ll[j],"nt)",
			# Generate length, determine if sample_points is the start or the end of the contig 

			this_contig_length=min(int(math.ceil(random.normalvariate(length_avg,length_sd))),ll[j])

			# print "with length",this_contig_length,
			# possible_samplings=[]

			if sample_points[i_idx]+this_contig_length <= ll_cs[j+1]: # can be the start 
				contig_to_sample.append((j,sample_points[i_idx]-ll_cs[j],sample_points[i_idx]+this_contig_length-ll_cs[j]))

			if sample_points[i_idx]-this_contig_length >= ll_cs[j]: # can be the end
				contig_to_sample.append((j,sample_points[i_idx]-this_contig_length-ll_cs[j],sample_points[i_idx]-ll_cs[j]))

			if this_contig_length>ll[j]/2:
				# can be the middle 
				left_end=sample_points[i_idx]-this_contig_length/2-ll_cs[j]
				right_end=sample_points[i_idx]-ll_cs[j]+this_contig_length/2
				if left_end<=0: 
					left_end=0
					right_end=this_contig_length
				elif right_end>=ll[j]:
					right_end=ll[j]
					left_end=ll[j]-this_contig_length
				contig_to_sample.append((j,left_end,right_end))

			last_start=j
			break



if AT_LEAST_ONE:
	guaranteed_contigs,remaining=class_balanced_sampling(contig_to_sample,0,1)
	contig_to_sample=remaining
else:
	guaranteed_contigs=[]

if BYSEQUENCE:
	# we sample at most n_contigs for each sequence 
	contig_to_sample,foo = class_balanced_sampling(contig_to_sample,0,n_contigs)
else:
	contig_to_sample=random.sample(contig_to_sample,min(n_contigs,len(contig_to_sample)))

contig_to_sample=set(contig_to_sample)
contig_to_sample.update(guaranteed_contigs)
contig_to_sample=list(contig_to_sample)

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


print "--"*12,len(contig_to_sample),"sampled contigs of mean size",scipy.mean([end-start for idx,start,end in contig_to_sample])
# print coverage 
for c in coverages:
	print scipy.mean(c),c


# Print contigs 

last_cont_idx=None
for cont_idx,start,end in contig_to_sample:
	if last_cont_idx!=cont_idx:
		print "\n\n"
		print SEQUENCE_CHARS[cont_idx]*ll[cont_idx],len([x for x in sequence_assigment if x==cont_idx]),"sample_points"
		last_cont_idx=cont_idx
	print " "*(start)+"."*(end-start)+" "*(ll[cont_idx]-end)
# Make statistics 
