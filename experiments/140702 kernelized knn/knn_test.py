# -*- coding: utf-8 -*-
from sklearn import neighbors
import pandas 
import pandas.rpy.common as com
import pandas.io.parsers
import scipy,scipy.sparse
import numpy as np

import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix




## Parse and build matrix 
input_kmers=pandas.io.parsers.read_csv("../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_11111.csv",sep="\t")

# make a dense matrix 
input_kmers_counts=pandas.pivot_table(input_kmers,values="count",index=['sequence_description'],columns=["kmer"],fill_value=0)

kmer_colums=input_kmers_counts.columns
# deduce class and species from seq description 

all_species=[]
all_classes=[]
for sd in input_kmers_counts.index:
	species,a_class,foo,contig_n,bar,bar=sd.split('_')
	all_species.append(species)
	all_classes.append(a_class)

input_kmers_counts["class"]=all_classes
input_kmers_counts["species"]=all_species

# make split training  / testing 

n_rows=input_kmers_counts.shape[0]
indices=np.random.permutation(n_rows)
training_ratio=0.5
training_set=indices[0:int(n_rows*training_ratio)]
testing_set=indices[-int(n_rows*training_ratio):]
training_data=input_kmers_counts.loc[training_set]
testing_data=input_kmers_counts.loc[testing_set]


# fit a KNN 

clf = neighbors.KNeighborsClassifier(15, weights="uniform")
clf.fit(training_data[kmer_colums], training_data["class"])
#print "predicting"
predicted_classes= clf.predict(testing_data[kmer_colums])
# compute kappa stat 
confusion_matrix(testing_data["class"],predicted_classes)
# make a mapping 
class_map=dict(zip(set(testing_data["class"]),range(0,4)))
kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])

## Using rpy and caret for the stats 

caret = importr("caret")
cm=caret.confusionMatrix(robjects.FactorVector(predicted_classes),robjects.FactorVector(testing_data["class"]))

com.convert_robj(cm[2])
com.convert_robj(cm[2])['Kappa']



# make that into an objective function


def kNNClass(train_idx,test_idx,n_neighbors):
	training_data=input_kmers_counts.loc[train_idx]
	testing_data=input_kmers_counts.loc[test_idx]
	clf = neighbors.KNeighborsClassifier(n_neighbors, weights="uniform")
	clf.fit(training_data[kmer_colums], training_data["class"])
	#print "predicting"
	predicted_classes= clf.predict(testing_data[kmer_colums])
	# compute kappa stat 
	confusion_matrix(testing_data["class"],predicted_classes)
	# make a mapping 
	class_map=dict(zip(set(testing_data["class"]),range(0,4)))
	kapp=kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])
	cm=caret.confusionMatrix(robjects.FactorVector(predicted_classes),robjects.FactorVector(testing_data["class"]))
	return kapp,cm

