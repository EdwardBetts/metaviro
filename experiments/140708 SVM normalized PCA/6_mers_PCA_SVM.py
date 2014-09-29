# -*- coding: utf-8 -*-
from sklearn import neighbors
import pandas 

from time import time
import numpy as np
import pandas.rpy.common as com
import pandas.io.parsers
import scipy,scipy.sparse
import numpy as np
import pylab as pl
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


## Parse and build matrix 
input_kmers=pandas.io.parsers.read_csv("../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_111111.csv",sep="\t")

# make a dense matrix 
input_kmers_counts=pandas.pivot_table(input_kmers,values="count",index=['sequence_description'],columns=["kmer"],fill_value=0)


kmer_colums=input_kmers_counts.columns

# Normalize by row counts 
row_sums=input_kmers_counts[kmer_colums].sum(axis=1)
normalized_counts=input_kmers_counts[kmer_colums].apply(lambda l:l/row_sums,axis=0)
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

## univariate feature selection 
means=normalized_counts[kmer_colums].apply(scipy.mean,0).order()
std_vars=normalized_counts[kmer_colums].apply(scipy.std,0).order()

# We do a center and scale and PCA 
normalized_counts[kmer_colums]=scale(normalized_counts[kmer_colums])
with open("6_mers_normalized_matrices.bdat","w") as f:
	cPickle.dump((normalized_counts,kmer_colums),f)

with open("6_mers_normalized_matrices.bdat") as f:
	normalized_counts,kmer_colums = cPickle.load(f)

X_train,X_test,Y_train,Y_test=train_test_split(normalized_counts[kmer_colums],normalized_counts["class"],test_size=0.4,random_state=42)


cv = StratifiedKFold(y=Y_train, n_folds=2)
pca=PCA()
rbfSVM=SVC()
pipe = Pipeline(steps=[('pca', pca), ('rbfSVM', rbfSVM)])

param_dist={
	"pca__n_components":sp_randint(10,250),
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
class_map=dict(zip(set(normalized_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predicted_held_out])


# We determine whether the variance of the number of components for the best CV 
all_scores=random_search.grid_scores_
all_scores.sort(key=lambda x:x.mean_validation_score)
with open("random_search_scores_1000iter_6mers.bdat","w") as f :
	cPickle.dump(all_scores,f)

# We try the best parameters on an indepedent sample with 80% training
X_train,X_test,Y_train,Y_test=train_test_split(normalized_counts[kmer_colums],normalized_counts["class"],test_size=0.2,random_state=50)
pipe_best = Pipeline(steps=[('pca', pca), ('rbfSVM', rbfSVM)])
# pipe_best=pipe_best.set_params(pca__n_components=131,rbfSVM__kernel="rbf",rbfSVM__C=4,rbfSVM__gamma=0.004)
pipe_best=pipe_best.set_params(
	pca__n_components=random_search.best_params_['pca__n_components'],
	rbfSVM__kernel=random_search.best_params_['rbfSVM__kernel'],
	rbfSVM__C=random_search.best_params_['rbfSVM__C'],
	rbfSVM__gamma=random_search.best_params_['rbfSVM__gamma']
	)
pipe_best.fit(X_train,Y_train)
pipe_best_pred_HO=pipe_best.predict(X_test)
mmat=confusion_matrix(pipe_best_pred_HO,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in pipe_best_pred_HO])


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
all_scores_df.to_csv("random_search_scores_1000iter_6mers.csv")