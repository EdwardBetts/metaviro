import math
import random
import scipy 
import collections
import numpy as np


# ll=[12,28,30,15,15]
NSEQUENCES=15000
MAXSEQLENGTH=90000
MINSEQLENGTH=20000
n_contigs=5000
length_avg=500
length_sd=20
BYSEQUENCE=False
AT_LEAST_ONE=False


SEQUENCE_CHARS=list("ABCDEFGHIJKLMNOPQRSTUVXYZabcdefghijklmnopqrstuvxyz")

# opt: Generate all length at once ? 
# opt: mode at least one contig per sequence
# opt: mode minimal/maximal length boundaries 
# opt: account for any distribution of contig length 
# opt: by sequence 




ll=np.random.randint(MINSEQLENGTH,MAXSEQLENGTH,NSEQUENCES)

coverages=[]
for i in range(len(ll)):
	coverages.append([0]*ll[i])





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

def sample_intervals(interval_lengths,n_contigs,BYSEQUENCE=False,AT_LEAST_ONE=False,length_avg=500,length_sd=200):
	tot_length=sum(interval_lengths)
	interval_lengths_cs=np.hstack(([0],np.cumsum(interval_lengths)))
	sample_points=[]
	if AT_LEAST_ONE:
		for i in range(len(interval_lengths_cs)-1):
			sample_points.append(random.randint(interval_lengths_cs[i],interval_lengths_cs[i+1]))


	if BYSEQUENCE:
		for i in range(len(interval_lengths_cs)-1):
			sample_points.extend(random.sample(xrange(interval_lengths_cs[i],interval_lengths_cs[i+1]),n_contigs))
	else:
		sample_points.extend(random.sample(xrange(0,tot_length),n_contigs))

	sample_points.sort()

	contig_to_sample=[]

	last_start=0
	sequence_assigment=[-1]*len(sample_points)

	# For each sample_point, find the corresponding interval, then generate three sample possibilities 
	# where the sample point can be either start, end, or middle

	for i_idx in xrange(len(sample_points)):
		for j in xrange(last_start,len(interval_lengths_cs)-1):
			if ((sample_points[i_idx] >= interval_lengths_cs[j]) and (sample_points[i_idx] < interval_lengths_cs[j+1])):
				sequence_assigment[i_idx]=j
				# Generate length, then determine if sample_points is the start or the end of the contig 

				this_contig_length=min(int(math.ceil(random.normalvariate(length_avg,length_sd))),interval_lengths[j])

				# print "with length",this_contig_length,
				# possible_samplings=[]

				if sample_points[i_idx]+this_contig_length <= interval_lengths_cs[j+1]: # can be the start 
					contig_to_sample.append((j,sample_points[i_idx]-interval_lengths_cs[j],sample_points[i_idx]+this_contig_length-interval_lengths_cs[j]))

				if sample_points[i_idx]-this_contig_length >= interval_lengths_cs[j]: # can be the end
					contig_to_sample.append((j,sample_points[i_idx]-this_contig_length-interval_lengths_cs[j],sample_points[i_idx]-interval_lengths_cs[j]))

				if this_contig_length>interval_lengths[j]/2:
					# can be the middle 
					left_end=sample_points[i_idx]-this_contig_length/2-interval_lengths_cs[j]
					right_end=sample_points[i_idx]-interval_lengths_cs[j]+this_contig_length/2
					if left_end<=0: 
						left_end=0
						right_end=this_contig_length
					elif right_end>=interval_lengths[j]:
						right_end=interval_lengths[j]
						left_end=interval_lengths[j]-this_contig_length
					contig_to_sample.append((j,left_end,right_end))

				last_start=j
				break


	# reduce the contig_to_sample to match the requested number of contigs

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
		assert(start<=interval_lengths[idx])
		assert(end<=interval_lengths[idx])
		assert(start>=0)


	# Sort 
	contig_to_sample.sort(key=lambda x: x[0:2])
	return contig_to_sample



contig_to_sample=sample_intervals(ll,n_contigs,BYSEQUENCE,AT_LEAST_ONE,length_avg,length_sd)

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
