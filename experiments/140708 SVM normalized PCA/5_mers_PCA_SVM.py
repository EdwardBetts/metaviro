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
input_kmers=pandas.io.parsers.read_csv("../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_11111.csv",sep="\t")

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


# We do a center and scale and PCA 


normalized_counts[kmer_colums]=scale(normalized_counts[kmer_colums])
normalized_counts[kmer_colums].apply(scipy.mean,0)
normalized_counts[kmer_colums].apply(scipy.std,0)
non_zero=(normalized_counts[kmer_colums]!=0).apply(scipy.sum,0)
too_abundant_kmers=list(non_zero.order()[-10:].index)
kmer_colums_filt=list(set(kmer_colums).difference(too_abundant_kmers))


pca_trans=PCA(n_components=100)
pca_fitted=pca_trans.fit(normalized_counts[kmer_colums])
pca_coord=pca_fitted.transform(normalized_counts[kmer_colums])


# We plot with some k-means 

kmeans = KMeans(init='k-means++', n_clusters=24, n_init=10)
kmeans.fit(pca_coord[:,0:2])

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 2     # point in the mesh [x_min, m_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = pca_coord[:, 0].min() - 0.1, pca_coord[:, 0].max() + 0.1
y_min, y_max = pca_coord[:, 1].min() - 0.1, pca_coord[:, 1].max() + 0.1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.figure(1)
pl.clf()
pl.imshow(Z, interpolation='nearest',
          extent=(xx.min(), xx.max(), yy.min(), yy.max()),
          cmap=pl.cm.Paired,
          aspect='auto', origin='lower')

# color mapping 
colors_m=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
colors=[colors_m[x] for x in input_kmers_counts['class']]
plt.scatter(pca_coord[:, 0], pca_coord[:, 1],c=colors)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_
pl.scatter(centroids[:, 0], centroids[:, 1],
           marker='x', s=169, linewidths=3,
           color='w', zorder=10)
pl.title('K-means clustering on the digits dataset (PCA-reduced data)\n'
         'Centroids are marked with white cross')
pl.xlim(x_min, x_max)
pl.ylim(y_min, y_max)
pl.xticks(())
pl.yticks(())
pl.show()


# SVM classification 


X_train,X_test,Y_train,Y_test=train_test_split(pca_coord,normalized_counts["class"],test_size=0.2,random_state=421)
clf=SVC()
clf.fit(X_train,Y_train)
predictions=clf.predict(X_test)
mmat=confusion_matrix(predictions,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predictions])






# Parameter tuning using grid search as in http://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html
X_train,X_test,Y_train,Y_test=train_test_split(pca_coord,normalized_counts["class"],test_size=0.2,random_state=42)
C_range = 10.0 ** np.arange(-2, 2)
gamma_range = 10.0 ** np.arange(-5, 0) # add 20
param_grid = dict(gamma=gamma_range, C=C_range)

cv = StratifiedKFold(y=Y_train, n_folds=4)
grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv,n_jobs=4)
grid.verbose=True
grid.fit(X_train, Y_train)

print "The best classifier is: ", grid.best_estimator_
# testing on the held out set 
predicted_held_out=grid.predict(X_test)
mmat=confusion_matrix(predicted_held_out,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predicted_held_out])




# We only keep K first dim; using the same model parameter as the best fit; inside a GridSearchCV 
X_train,X_test,Y_train,Y_test=train_test_split(normalized_counts[kmer_colums],normalized_counts["class"],test_size=0.4,random_state=42)
cv = StratifiedKFold(y=Y_train, n_folds=2)
pca=PCA()
rbfSVM=SVC()
n_components = [5,10,50,100,200,300,400,500,600,700]
C_range = [10.0,20,30]
gamma_range = [0.001,0.01]

pipe = Pipeline(steps=[('pca', pca), ('rbfSVM', rbfSVM)])

estimator=GridSearchCV(pipe,param_grid=dict(pca__n_components=n_components,rbfSVM__C=C_range,rbfSVM__gamma=gamma_range),cv=cv,verbose=5,n_jobs=4)
estimator.fit(X_train,Y_train)


predicted_held_out=estimator.predict(X_test)
mmat=confusion_matrix(predicted_held_out,Y_test)
print mmat
class_map=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
kappa([class_map[x] for x in Y_test],[class_map[x] for x in predicted_held_out])


# Same approach using a randomized param search (http://scikit-learn.org/stable/auto_examples/randomized_search.html#example-randomized-search-py)

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

# We try the best parameters on indepedent samples 
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
all_scores_df.to_csv("random_search_scores_1000iter_5mers.csv")