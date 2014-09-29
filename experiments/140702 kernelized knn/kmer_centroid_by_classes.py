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
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale


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


classes_mapping=dict(zip(set(input_kmers_counts["class"]),range(0,4)))
def map_classes(c):
	return classes_mapping(c)


# make split training  / testing 

n_rows=input_kmers_counts.shape[0]
indices=np.random.permutation(n_rows)
training_ratio=0.2
training_set=indices[0:int(n_rows*training_ratio)]
testing_set=indices[-int(n_rows*training_ratio):]

training_data=normalized_counts.loc[training_set]
testing_data=normalized_counts.loc[testing_set]



### K-means 

def bench_k_means(estimator, name, data,labels):
	t0 = time()
	estimator.fit(data)
	print('% 9s   %.2fs    %i   %.3f   %.3f   %.3f   %.3f   %.3f    %.3f'
		  % (name, (time() - t0), estimator.inertia_,
			metrics.homogeneity_score(labels, estimator.labels_),
			metrics.completeness_score(labels, estimator.labels_),
			metrics.v_measure_score(labels, estimator.labels_),
			metrics.adjusted_rand_score(labels, estimator.labels_),
			metrics.adjusted_mutual_info_score(labels,  estimator.labels_),
			metrics.silhouette_score(data, estimator.labels_, metric='euclidean')))

bench_k_means(KMeans(init='k-means++', n_clusters=4, n_init=10), name="k-means++", data=training_data[kmer_colums],labels=training_data["class"])

bench_k_means(KMeans(init='random', n_clusters=4, n_init=10), name="random", data=training_data[kmer_colums],labels=training_data["class"])

# in this case the seeding of the centers is deterministic, hence we run the
# kmeans algorithm only once with n_init=1
pca = PCA(n_components=4).fit(training_data[kmer_colums])
bench_k_means(KMeans(init=pca.components_, n_clusters=4, n_init=1), name="PCA-based", data=training_data[kmer_colums],labels=training_data["class"])


## plot of the boundaries 

reduced_data = PCA(n_components=2).fit_transform(input_kmers_counts[kmer_colums])
kmeans = KMeans(init='k-means++', n_clusters=24, n_init=10)
kmeans.fit(reduced_data)

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 2     # point in the mesh [x_min, m_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() - 0.1, reduced_data[:, 0].max() + 0.1
y_min, y_max = reduced_data[:, 1].min() - 0.1, reduced_data[:, 1].max() + 0.1
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
plt.scatter(reduced_data[:, 0], reduced_data[:, 1],c=colors)
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


#Cluster analysis using PCA 

reduced_data=PCA(n_components=32).fit(input_kmers_counts[kmer_colums])
pl.plot(reduced_data.explained_variance_ratio_)
pl.show(block=False)

plt.pcolor(reduced_data.components_)
pl.show()

# what are the loadings of the 1st component  ? 

loadings=pandas.DataFrame({"loadings":reduced_data.components_[0,:],"kmer":kmer_colums})
loadings.sort('loadings')

loadings=pandas.DataFrame({"loadings":reduced_data.components_[1,:],"kmer":kmer_colums})
loadings.sort('loadings')


# Scatter plot of GC rich vs AT rich 


gc_rich_kmers= [x for x in kmer_colums if set(x)==set(['G','C'])]
at_rich_kmers= [x for x in kmer_colums if set(x)==set(['A','T'])]

gc_rich=input_kmers_counts[gc_rich_kmers].apply(lambda x:sum(x),axis=1)
at_rich=input_kmers_counts[at_rich_kmers].apply(lambda x:sum(x),axis=1)
kmer_gc_only=pandas.concat((gc_rich,at_rich),axis=1)


ppl.scatter(kmer_gc_only[0], kmer_gc_only[1],c=colors)
pl.show()

# We make a data frame with the PCA coordinates and all annotations 


input_kmers_counts_output=pandas.pivot_table(input_kmers,values="count",index=['sequence_description','sequence_length','GC'],columns=["kmer"],fill_value=0)


reduced_data_coord=reduced_data.fit_transform(input_kmers_counts[kmer_colums])

output_table=pandas.DataFrame({
	"sequence_description":input_kmers_counts.index,
	"coord_1":reduced_data_coord[:,0],
	"coord_2":reduced_data_coord[:,1],
	# "gc_content":input_kmers_counts
	"class":input_kmers_counts['class']
	})
output_table=pandas.merge(output_table,input_kmers_counts_output,left_index=True,right_index=True)

output_table.to_csv('pca_coord_3mers.csv')


# Re-do PCA/kMeans on length normalized data 
## plot of the boundaries 

reduced_data = PCA(n_components=2).fit_transform(normalized_counts[kmer_colums])
kmeans = KMeans(init='k-means++', n_clusters=24, n_init=10)
kmeans.fit(reduced_data)

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 0.02     # point in the mesh [x_min, m_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() - 0.1, reduced_data[:, 0].max() + 0.1
y_min, y_max = reduced_data[:, 1].min() - 0.1, reduced_data[:, 1].max() + 0.1
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
plt.scatter(reduced_data[:, 0], reduced_data[:, 1],c=colors)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_
pl.scatter(centroids[:, 0], centroids[:, 1],
           marker='x', s=169, linewidths=3,
           color='w', zorder=10)
pl.title('K-means clustering on the  3 mers (PCA-reduced on length normalized)\n'
         'Centroids are marked with white cross')
pl.xlim(x_min, x_max)
pl.ylim(y_min, y_max)
pl.xticks(())
pl.yticks(())
pl.show()


input_kmers_counts_output=pandas.pivot_table(input_kmers,values="count",index=['sequence_description','sequence_length','GC'],columns=["kmer"],fill_value=0)

row_sums=input_kmers_counts_output[kmer_colums].sum(axis=1)
normalized_counts=input_kmers_counts_output[kmer_colums].apply(lambda l:l/row_sums,axis=0)

output_table=pandas.DataFrame({
	"sequence_description":input_kmers_counts.index,
	"coord_1":reduced_data[:,0],
	"coord_2":reduced_data[:,1],
	# "gc_content":input_kmers_counts
	"class":input_kmers_counts['class']
	})
output_table=pandas.merge(pandas.merge(output_table,input_kmers_counts_output,left_index=True,right_index=True),normalized_counts,left_index=True,right_index=True)

output_table.to_csv('pca_coord_length_normalized_3mers.csv')



## We normalize by GC too 

input_kmers_counts=pandas.pivot_table(input_kmers,values="count",index=['sequence_description'],columns=["kmer"],fill_value=0)
kmer_colums=input_kmers_counts.columns

all_species=[]
all_classes=[]
for sd in input_kmers_counts.index:
	species,a_class,foo,contig_n,bar,bar=sd.split('_')
	all_species.append(species)
	all_classes.append(a_class)



contigs_annot=input_kmers[['sequence_description',"sequence_length",'GC']].drop_duplicates()
contigs_annot=contigs_annot.set_index("sequence_description")
# make them the same order
contigs_annot =contigs_annot.ix[input_kmers_counts.index]

# normalized_counts=input_kmers_counts[kmer_colums].apply(lambda l:l/contigs_annot['GC'],axis=0)
normalized_counts=input_kmers_counts[kmer_colums].apply(lambda l:l/contigs_annot['sequence_length'],axis=0)
reduced_data = PCA(n_components=4).fit_transform(normalized_counts[kmer_colums])

output_table=pandas.DataFrame({
	"sequence_description":normalized_counts.index,
	"coord_1":reduced_data[:,0],
	"coord_2":reduced_data[:,1],
	"coord_3":reduced_data[:,2],
	"coord_4":reduced_data[:,3],
	# "gc_content":input_kmers_counts
	"class":all_classes,
	"species":all_species
	})

output_table=output_table.set_index('sequence_description')
output_table=pandas.merge(output_table,contigs_annot,left_index=True,right_index=True)
output_table=pandas.merge(output_table,normalized_counts,left_index=True,right_index=True)
output_table=pandas.merge(output_table,input_kmers_counts,left_index=True,right_index=True)

output_table.to_csv('pca_coord_length_gc_normalized_3mers.csv')