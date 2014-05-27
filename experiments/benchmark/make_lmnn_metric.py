# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 17:49:16 2014

@author: lschmitt
"""
from metric_learn.lmnn import LMNN
import sys, os.path
import numpy as np
try:
   import cPickle as pickle
except:
   import pickle
   

filename = sys.argv[1] # first argument: data for a specific kmer length (.bdat file)
outfile = os.path.splitext(filename)[0] + ".metric"

print "~~~ %s ~~~" %filename

print "loading matrix"
with open(filename, "r") as fh:
    tag, length, targ, matrix, training_indices = pickle.load(fh)

# normalize data (only once ffs)
for i in range(len(length)):
    matrix[i,:]=matrix[i,:]/length[i]

# lmnnify the space (take that, Vador)
print "distorting space"
metric = LMNN(np.asarray(matrix), targ)
metric.fit()

print "serialize metric"
with open(outfile, 'w') as outfh:
    pickle.dump(metric, outfh)

print "done."
