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
	species,a_class,foo,contig_n,bar,bar=sd.split('_')
	all_species.append(species)
	all_classes.append(a_class)

input_kmers_counts["class"]=all_classes
input_kmers_counts["species"]=all_species
normalized_counts["class"]=all_classes
normalized_counts["species"]=all_species

# PCA 
normalized_counts[kmer_colums]=scale(normalized_counts[kmer_colums])
normalized_counts[kmer_colums].apply(scipy.mean,0)
normalized_counts[kmer_colums].apply(scipy.std,0)
non_zero=(normalized_counts[kmer_colums]!=0).apply(scipy.sum,0)
too_abundant_kmers=list(non_zero.order()[-10:].index)
kmer_colums_filt=list(set(kmer_colums).difference(too_abundant_kmers))


pca_trans=PCA(n_components=160)
pca_fitted=pca_trans.fit(normalized_counts[kmer_colums])
pca_coord=pca_fitted.transform(normalized_counts[kmer_colums])

# SVM classification and kappa estimates 

X_train,X_test,Y_train,Y_test=train_test_split(pca_coord,normalized_counts["class"],test_size=0.5,random_state=421)
clf=SVC(C=4.152687927300392,gamma=0.002448996894369464,kernel='rbf')
clf.fit(X_train,Y_train)
predictions=clf.predict(X_test)
mmat=confusion_matrix(predictions,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predictions])

