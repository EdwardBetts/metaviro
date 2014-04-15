# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""
# create matrices from mv file

from __future__ import with_statement
import time, re
from scipy.sparse import coo_matrix, csr_matrix, dok_matrix, lil_matrix
from sklearn.metrics.pairwise import euclidean_distances

import numpy as np
import scipy as sp
import random, multiprocessing, Queue, sys, time, subprocess, os.path
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets

from sklearn.metrics import confusion_matrix, f1_score, SCORERS
from skll.metrics import kappa


try:
   import cPickle as pickle
except:
   import pickle

filename = sys.argv[1]

#random.seed(123456)

n_neighbors = 15

print "~~~ %s ~~~" %filename

print "loading matrix"
with open(filename, "r") as fh:
    tag, length, targ, matrix, training_indices = pickle.load(fh)

print "sampling"




training = training_indices[0][0]
#training = random.sample(training_indices[0][0],20000)
testing = training_indices[0][1]
#testing = random.sample(training_indices[0][1],5000)

#print ("### sans norm")
#
#training_matrix = matrix[training,:]
#training_targets = [targ[i] for i in training]
#
#testing_matrix = matrix[testing,:]
#testing_targets = [targ[i] for i in testing]
#
#
#print "fitting"
#
#clf = neighbors.KNeighborsClassifier(n_neighbors, weights="distance")
#clf.fit(training_matrix, training_targets)
#print "predicting"
#classes=clf.predict(testing_matrix)
#
#print(confusion_matrix(classes,testing_targets))
            #
#print(kappa(classes,testing_targets))

#print ("### avec norm")

for i in range(len(length)):
    matrix[i,:]=matrix[i,:]/length[i]

training_matrix = matrix[training,:]
training_targets = [targ[i] for i in training]

testing_matrix = matrix[testing,:]
testing_targets = [targ[i] for i in testing]


#print "fitting"

clf = neighbors.KNeighborsClassifier(n_neighbors, weights="distance")
clf.fit(training_matrix, training_targets)
#print "predicting"
classes=clf.predict(testing_matrix)

print(confusion_matrix(classes,testing_targets))
            
print(kappa(classes,testing_targets))
