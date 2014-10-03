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

SEQUENCE_CHARS=list("ABCDEFGHIJKLMNOPQRSTUVXYZabcdefghijklmnopqrstuvxyz")


ll=np.random.randint(MINSEQLENGTH,MAXSEQLENGTH,NSEQUENCES)
coverages=[]
for i in range(len(ll)):
	coverages.append([0]*ll[i])

tot_length=sum(ll)

ll_cs=np.hstack(([0],np.cumsum(ll)))




sample_points=sorted(random.sample(xrange(0,tot_length),n_contigs))
# set(sample_points).intersection(ll_cs)
# assign each sample_points to a sequence 


# opt: Generate all length at once ? 
# opt: mode at least one contig per sequence
# opt: mode minimal/maximal length boundaries 
# opt: account for any distribution of contig length 
# opt: by sequence 

contig_to_sample=[]

sequence_assigment=[-1]*len(sample_points)
for i_idx in xrange(len(sample_points)):
	# print i_idx,":",sample_points[i_idx],
	for j in xrange(0,len(ll_cs)):
		if ((sample_points[i_idx] >= ll_cs[j]) and (sample_points[i_idx] < ll_cs[j+1])):
			sequence_assigment[i_idx]=j
			# print "Found in [",ll_cs[j],ll_cs[j+1],"](",ll[j],"nt)",
			# Generate length, determine if sample_points is the start or the end of the contig 

			this_contig_length=min(int(math.ceil(random.normalvariate(length_avg,length_sd))),ll[j])

			# print "with length",this_contig_length,
			possible_samplings=[]

			if sample_points[i_idx]+this_contig_length <= ll_cs[j+1]: # can be the start 
				possible_samplings.append([j,sample_points[i_idx]-ll_cs[j],sample_points[i_idx]+this_contig_length-ll_cs[j]])

			if sample_points[i_idx]-this_contig_length >= ll_cs[j]: # can be the end
				possible_samplings.append([j,sample_points[i_idx]-this_contig_length-ll_cs[j],sample_points[i_idx]-ll_cs[j]])

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

				possible_samplings.append([j,left_end,right_end])

			# # correct for too short sequences and boundaries conditions
			# correct_sampling_units=[]
			# for a_samp in possible_samplings:

			# 	# Sequence is too short for the contig length, we truncate 
			# 	new_samp=list(a_samp)
			# 	if new_samp[1]<0:
			# 		new_samp[1]=0
			# 	# Sequence is too short for the contig length, we truncate 
			# 	if new_samp[2]>ll[j]:
			# 		new_samp[2]=ll[j]
			# 		# correct_sampling_units.append(new_samp)
			# 	if (a_samp[1]>=0) and (a_samp[2]<=ll[j]): #is valid, we keep as it is 
			# 		correct_sampling_units.append(a_samp)

			# 	correct_sampling_units.append(new_samp)

			# final_samp_unit=random.choice(correct_sampling_units)
			final_samp_unit=random.choice(possible_samplings)

			contig_to_sample.append(final_samp_unit)
			# print "=> will sample",sampling_unit

			break

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
		print SEQUENCE_CHARS[cont_idx]*ll[cont_idx],len([x for x in sequence_assigment if x==cont_idx]),"sample_points"
		last_cont_idx=cont_idx
	print " "*(start)+"."*(end-start)+" "*(ll[cont_idx]-end)
# Make statistics 
