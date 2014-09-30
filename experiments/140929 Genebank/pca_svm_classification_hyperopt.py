# -*- coding: utf-8 -*-
from sklearn import neighbors
import pandas 
import prettyplotlib as ppl

import matplotlib.pyplot as plt


from time import time
import numpy as np
import pandas.rpy.common as com
import pandas.io.parsers
import scipy,scipy.sparse
import numpy as np

import pylab as pl


import prettytable

import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import  train_test_split
from sklearn.grid_search import GridSearchCV,RandomizedSearchCV
from scipy.stats import randint as sp_randint

from sklearn.svm import SVC 

import logging
import threading
import itertools
import cPickle


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
input_kmers=pandas.io.parsers.read_csv("../../data/Genebank/GB_20mb_per_domain_k5_long.csv",sep="\t")

# make a dense matrix 
input_kmers_counts=pandas.pivot_table(input_kmers,values="count",index=['sequence_description'],columns=["kmer"],fill_value=0)


kmer_colums=input_kmers_counts.columns

# Normalize by row counts 
row_sums=input_kmers_counts[kmer_colums].sum(axis=1)
normalized_counts=input_kmers_counts.apply(lambda l:l/row_sums,axis=0)
# for col in kmer_colums: 
# 	input_kmers_counts[col]=input_kmers_counts[col]/row_sums


all_species=[]
all_classes=[]
for sd in input_kmers_counts.index:
	components=sd.split('_')
	species=components[0]
	a_class=components[1]
	contig_n=components[3]
	all_species.append(species)
	all_classes.append(a_class)

input_kmers_counts["class"]=all_classes
input_kmers_counts["species"]=all_species
normalized_counts["class"]=all_classes
normalized_counts["species"]=all_species

# Hyperopt search for RBF SVM using 60%training 2 folds 

X_train,X_test,Y_train,Y_test=train_test_split(normalized_counts[kmer_colums],normalized_counts["class"],test_size=0.4,random_state=42)
cv = StratifiedKFold(y=Y_train, n_folds=2)

pca=PCA()
rbfSVM=SVC()


pipe = Pipeline(steps=[('pca', pca), ('rbfSVM', rbfSVM)])


param_dist={
	"pca__n_components":sp_randint(10,700),
	"rbfSVM__C": scipy.stats.expon(scale=10),
	"rbfSVM__kernel": ["rbf"], 
	"rbfSVM__gamma": scipy.stats.expon(scale=0.01)
}
n_iter_search = 500
random_search = RandomizedSearchCV(pipe, param_distributions=param_dist, n_iter=n_iter_search,cv=cv,verbose=6,n_jobs=4)


random_search.fit(X_train,Y_train)
predicted_held_out=random_search.predict(X_test)
mmat=confusion_matrix(predicted_held_out,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predicted_held_out])


# We determine whether the variance of the number of components for the best CV 
all_scores=random_search.grid_scores_
all_scores.sort(key=lambda x:x.mean_validation_score)
with open("random_search_scores_1000iter_5mers.bdat","w") as f :
	cPickle.dump(all_scores,f)

# We generate a pandas data.frame with the results 
import pandas

all_score_dict=[]
for score in all_scores:
	this_score=score.parameters
	this_score['mean_cv']=score.mean_validation_score
	for i in range(len(score.cv_validation_scores)):
		this_score['cv_%s'%(i)]=score.cv_validation_scores[i]
	all_score_dict.append(this_score)

all_scores_df=pandas.DataFrame(all_score_dict)
all_scores_df.to_csv("random_search_scores_1000iter_5mers.csv")