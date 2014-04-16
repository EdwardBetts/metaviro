# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:12:18 2014

@author: lschmitt
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""
# create matrices from mv file

#from __future__ import with_statement
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
from sklearn import svm

try:
   import cPickle as pickle
except:
   import pickle

filename = sys.argv[1] # first argument: data for a specific kmer length (.bdat file)
classifiers = sys.argv[2:] # other arguments: classifiers names as in classif_<name>.py with options
# example with options: rbf_svc,gamma=10,C=1

outfile = os.path.splitext(filename)[0] + ".csv"
if outfile.startswith('sample_matrix_'):
    outfile = outfile[len('sample_matrix_'):]
outfile = "-".join(classifiers) + "-" + outfile

print "~~~ %s ~~~" %filename

print "loading matrix"
with open(filename, "r") as fh:
    tag, length, targ, matrix, training_indices = pickle.load(fh)

with open("filtered_contigs.bdat", "r") as fh:
    filter_dict = pickle.load(fh)

with open(outfile, 'w') as outfh:
    for classifierstring in classifiers:
        print "Running ", classifierstring
        classifier = classifierstring.split(',')
        options = None
        if len(classifier) > 1:
            options = classifier[1:]
        exec "from classif_%s import classify as classify" % classifier[0]
        for batch_id in range(len(training_indices)):
            print "batch ", batch_id+1,'/',len(training_indices)
        
            training = training_indices[batch_id][0]
            testing = training_indices[batch_id][1]
            
            
            for i in range(len(length)):
                matrix[i,:]=matrix[i,:]/length[i]
            
            training_matrix = matrix[training,:]
            training_targets = [targ[i] for i in training]
            
            testing_matrix = matrix[testing,:]
            testing_targets = [targ[i] for i in testing]
                        
            testing_tags = [tag[i] for i in testing]
                        
            classes = classify(training_matrix, training_targets, testing_matrix, options)
            
            for i in range(len(testing_tags)):
                outfh.write( '\t'.join([str(batch_id), classifierstring, testing_tags[i], str(testing_targets[i]), str(classes[i]), str(filter_dict[testing_tags[i]])]) + '\n')

print "Done."
