# -*- coding: utf-8 -*-
from sklearn import neighbors
import pandas 
import pandas.rpy.common as com
import pandas.io.parsers
import scipy,scipy.sparse
import numpy as np

import prettytable

import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

import logging
import threading
import itertools

from Queue import *
NMAXTHREADS=6

if "logger" not in globals():
	logger = logging.getLogger('test_logger')
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s',"%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	logger.addHandler(ch)

logger = logging.getLogger('test_logger')
caret = importr("caret")

## Multi proc support
# Multi threading support, generic 

class Worker(threading.Thread):
	def __init__(self, function, in_queue, out_queue):
		self.function = function
		self.in_queue, self.out_queue = in_queue, out_queue
		super(Worker, self).__init__()

	def run(self):
		while True:
			try:
				if self.in_queue.empty(): 
					break
				data = self.in_queue.get()
				result = self.function(data)
				self.out_queue.put((data,result))
				self.in_queue.task_done()
				logger.info("Still %d to do",self.in_queue.qsize())
			except Exception as e:
				logger.critical('something happened!: Error on %s, %s',repr(data),repr(e))
				self.out_queue.put({})
				self.in_queue.task_done()
				break

def process(data, function, num_workers=1):
	in_queue = Queue()
	for item in data:
		in_queue.put(item)
	out_queue = Queue(maxsize=in_queue.qsize())
	print in_queue,in_queue.unfinished_tasks
	workers = [Worker(function, in_queue, out_queue) for i in xrange(num_workers)]
	for worker in workers: 
		worker.setDaemon(True)
		worker.start()
	in_queue.join()
	return out_queue


## Parse and build matrix 
input_kmers=pandas.io.parsers.read_csv("../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_111.csv",sep="\t")

# make a dense matrix 
input_kmers_counts=pandas.pivot_table(input_kmers,values="count",index=['sequence_description'],columns=["kmer"],fill_value=0)


kmer_colums=input_kmers_counts.columns

# Normalize by row counts 
row_sums=input_kmers_counts[kmer_colums].sum(axis=1)
normalized_counts=input_kmers_counts.apply(lambda l:l/row_sums,axis=0)
# for col in kmer_colums: 
# 	input_kmers_counts[col]=input_kmers_counts[col]/row_sums


# deduce class and species from seq description 

all_species=[]
all_classes=[]
for sd in input_kmers_counts.index:
	species,a_class,foo,contig_n,bar,bar=sd.split('_')
	all_species.append(species)
	all_classes.append(a_class)

input_kmers_counts["class"]=all_classes
input_kmers_counts["species"]=all_species
normalized_counts["class"]=all_classes
normalized_counts["species"]=all_species

# make split training  / testing 

n_rows=input_kmers_counts.shape[0]

indices=np.random.permutation(n_rows)
training_ratio=0.8
training_set=indices[0:int(n_rows*training_ratio)]

### HUGE BUG here making results totally invalid 
testing_set=indices[-int(n_rows*training_ratio):]



# fit a KNN on the normalized_counts
training_data=normalized_counts.loc[training_set]
testing_data=normalized_counts.loc[testing_set]

# clf = neighbors.KNeighborsClassifier(15, weights="uniform")
# clf.fit(training_data[kmer_colums], training_data["class"])
# #print "predicting"
# predicted_classes= clf.predict(testing_data[kmer_colums])
# # compute kappa stat 
# confusion_matrix(testing_data["class"],predicted_classes)
# # make a mapping 
# class_map=dict(zip(set(testing_data["class"]),range(0,4)))
# kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])

# ## Using rpy and caret for the stats 


# cm=caret.confusionMatrix(robjects.FactorVector(predicted_classes),robjects.FactorVector(testing_data["class"]))

# com.convert_robj(cm[2])
# com.convert_robj(cm[2])['Kappa']



# Test with a subset of k-mers 

# k_mer_subset=np.random.permutation(kmer_colums)[:16]

# training_subset=normalized_counts.loc[training_set][np.append(k_mer_subset,"class")]
# testing_subset=normalized_counts.loc[testing_set][np.append(k_mer_subset,"class")]

# clf = neighbors.KNeighborsClassifier(15, weights="uniform")
# clf.fit(training_subset[k_mer_subset], training_subset["class"])
# #print "predicting"
# predicted_classes= clf.predict(testing_subset[k_mer_subset])
# # compute kappa stat 
# confusion_matrix(testing_data["class"],predicted_classes)
# kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])

# make that into an objective function


def kNNClass(train_idx,test_idx,n_neighbors,k_mer_subset):
	logger.info('computing for %s'%(k_mer_subset))
	train_idx=train_idx
	test_idx=test_idx
	training_subset=normalized_counts.loc[train_idx][np.append(k_mer_subset,"class")]
	testing_subset=normalized_counts.loc[test_idx][np.append(k_mer_subset,"class")]
	clf = neighbors.KNeighborsClassifier(n_neighbors, weights="uniform")
	clf.fit(training_subset[k_mer_subset], training_subset["class"])
	#print "predicting"
	predicted_classes= clf.predict(testing_data[k_mer_subset])
	# compute kappa stat 
	confusion_matrix(testing_data["class"],predicted_classes)
	# make a mapping 
	class_map=dict(zip(set(testing_data["class"]),range(0,4)))
	kapp=kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])
	cm=caret.confusionMatrix(robjects.FactorVector(predicted_classes),robjects.FactorVector(	testing_data["class"]))
	logger.info("Finished for %s with kappa==%f"%(k_mer_subset,kapp))
	return kapp,cm


k_mer_subset=np.random.permutation(kmer_colums)[:32]
kNNClass(training_set,testing_set,15,k_mer_subset)



# we make it m-proc
all_kmer_subsets=[]
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:8] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:24] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:32] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:48] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:50] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:52] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:54] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:56] for i in range(60)])
# all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:58] for i in range(60)])
all_kmer_subsets.extend([np.random.permutation(kmer_colums)[:60] for i in range(60)])
all_kmer_subsets.append(kmer_colums)

inputs=itertools.product([training_set],[testing_set],[15],all_kmer_subsets)

results=process(inputs,lambda x:kNNClass(*x),4)


# we make it a pretty table 
k_mer_index = dict(zip(kmer_colums,range(0,len(kmer_colums))))
def map_to_str(l):
	indices = sorted([k_mer_index[x] for x in l])
	binary = [1 if i in indices else 0 for i in range(len(kmer_colums))]
	return ''.join(map(str,binary))


pp=prettytable.PrettyTable(['n kmers','kMers','Kappa'])
for q in results.queue:
	pp.add_row((len(q[0][3]),map_to_str(q[0][3]),q[1][0]))


print pp.get_string(sortby='Kappa')