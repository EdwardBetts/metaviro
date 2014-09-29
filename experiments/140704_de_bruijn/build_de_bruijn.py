#! /usr/bin/env python
# encoding: utf-8
from sklearn import neighbors
import pandas 
import pandas.rpy.common as com
import pandas.io.parsers
import scipy,scipy.sparse
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

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

import numpy as np
import networkx
import random

import os,sys
import itertools
import pandas

import logging
if "logger" not in globals():
	logger = logging.getLogger('spaced_kmer')
	logger.setLevel(logging.DEBUG)

	# while len(logger.handlers()) > 0:
	#  	logger.pop()

	# create console handler and set level to debug
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s',"%Y-%m-%d %H:%M:%S")
	# formatter = logging.Formatter('%(asctime)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)

	# add ch to logger
	logger.addHandler(ch)




caret = importr("caret")

# read the input fasta, gather contigs annotations

from Bio import SeqUtils
from Bio import SeqIO

from Bio.Alphabet import generic_dna

input_fasta="../140613_more_patterns/long_kmers_mVar_n30.fa"

uid=0
# all_records=pandas.DataFrame(columns=('idx','species','class','sequence_description','sequence_length','GC','sequence'))
all_records=[]

all_sequences=[]
for record in SeqIO.parse(input_fasta, "fasta", generic_dna):
	dict_seq={}
	dict_seq["idx"]=uid
	dict_seq["sequence"]=str(record.seq)
	# dict_seq["class"]
	species,this_class,foo,sample_id=record.description.split('_')[:4]
	dict_seq['species']=species
	dict_seq["class"]=this_class
	dict_seq['sequence_length']=len(record.seq)
	dict_seq["GC"]=SeqUtils.GC(record.seq)
	dict_seq["sample_id"]=int(sample_id)
	dict_seq["sequence_description"]=record.description
	# fasta_keys=record.description.split('_')+[record.description,str(len(seq))]
	uid+=1
	# all_records.append(fasta_keys)
	# all_sequences.append(seq)
	all_records.append(dict_seq)
all_records=pandas.DataFrame(all_records)
all_records.set_index("idx")




# We choose random sequences and k merize them in a graph, then compute the degree
all_nodes_df=pandas.DataFrame()
random_sequences=all_records.ix[np.random.permutation(range(len(all_records)))[0:15000]]
kmer_len=5
all_graphs=[]
for i in random_sequences.index:
	node_attributes={}
	this_g=networkx.DiGraph()
	this_record=random_sequences.ix[i]
	logger.info("processing %s"%(this_record))
	this_seq=this_record["sequence"]
	if(len(this_seq)-2 < kmer_len):
		continue

	previous_kmer=None
	for i in range(0,len(this_seq)-kmer_len):
		this_kmer=str(this_record["idx"])+'_'+str(this_record["class"])+"_"+this_seq[i:i+kmer_len]
		if(previous_kmer):
			this_g.add_edge(previous_kmer,this_kmer)
		previous_kmer=this_kmer
		node_attributes[this_kmer]={'species':this_record["species"],'class':this_record["class"],"kmer":this_kmer.split("_")[2],"sequence_description":this_record["sequence_description"]}
	for n in this_g.node:
		this_g.node[n].update(node_attributes[n])

#	all_graphs.append(this_g)
	# networkx.write_gml(this_g,'test_%dmer_g.gml'%(kmer_len))

	node_df=pandas.DataFrame(node_attributes.values())
	node_df=pandas.concat([node_df,pandas.DataFrame(node_attributes.keys())],axis=1)
	# node_df.columns=["class","kmer","species","kmer_key"]

	degree=networkx.degree(this_g)
	sorted(degree.items(),key=lambda x:x[1])[-20:]
	degree_df=pandas.DataFrame(degree.items())

	# betweenness_centrality=networkx.betweenness_centrality(this_g)
	# sorted(betweenness_centrality.items(),key=lambda x:x[1])[-20:]
	# betweenness_centrality_df=pandas.DataFrame(betweenness_centrality.items())

	node_df=pandas.merge(node_df,degree_df)
	node_df.columns=["class","kmer","sequence_description","species","kmer_key","degree"]
	all_nodes_df=pandas.concat([all_nodes_df,node_df])


# networkx.write_gml(networkx.union_all(all_graphs),'test_%dmer_g.gml'%(kmer_len))

# we try a kNN based on the degree distribution 

all_nodes_df.to_csv('sparse_matrix_%dmers_10k_contigs.csv'%(kmer_len),index=False)


def kNNClass(train_idx,test_idx,n_neighbors,k_mer_subset):
	logger.info('computing for %s%'(k_mer_subset))
	train_idx=train_idx
	test_idx=test_idx
	training_subset=input_kmers_counts.loc[train_idx][np.append(k_mer_subset,"class")]
	testing_subset=input_kmers_counts.loc[test_idx][np.append(k_mer_subset,"class")]
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


input_kmers_counts=pandas.pivot_table(all_nodes_df,values="degree",index=['sequence_description'],columns=["kmer"],fill_value=0)
input_kmers_counts.to_csv('degree_matrix_%dmers_10k_contigs.csv'%(kmer_len))
count_colums=input_kmers_counts.columns



all_species=[]
all_classes=[]
for sd in input_kmers_counts.index:
	species,a_class,foo,contig_n,bar,bar=sd.split('_')
	all_species.append(species)
	all_classes.append(a_class)

input_kmers_counts["class"]=all_classes
input_kmers_counts["species"]=all_species

n_rows=len(input_kmers_counts)

indices=np.random.permutation(n_rows)
training_ratio=0.8
training_set=indices[0:int(n_rows*training_ratio)]
testing_set=indices[-int(n_rows*training_ratio):]

training_data=input_kmers_counts.loc[training_set]
testing_data=input_kmers_counts.loc[testing_set]


clf = neighbors.KNeighborsClassifier(15, weights="uniform")
clf.fit(training_data[count_colums], training_data["class"])
#print "predicting"
predicted_classes= clf.predict(testing_data[count_colums])
# compute kappa stat 
confusion_matrix(testing_data["class"],predicted_classes)
# make a mapping 
class_map=dict(zip(set(testing_data["class"]),range(0,4)))
kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])


# fit a KNN on the normalized_counts
# kNNClass(training_set,testing_set,15,count_colums)


# We focus on the ambiguous k-mers, approx 15k; basically all-kmer appear more than once

ambiguous_kmers=all_nodes_df[all_nodes_df["degree"]>2]
len(set(ambiguous_kmers['kmer']))
len(set(all_nodes_df['kmer']))

# We do a PCA on that 
amb_kmers_counts=pandas.pivot_table(ambiguous_kmers,values="degree",index=['sequence_description'],columns=["kmer"],fill_value=0)
kmer_colums=amb_kmers_counts.columns

# center scale and normalize
pca = PCA(n_components=16).fit(training_data[kmer_colums])