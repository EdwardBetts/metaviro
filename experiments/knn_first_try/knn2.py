# create matrices from mv file

from __future__ import with_statement
import time, re
from scipy.sparse import coo_matrix, csr_matrix, dok_matrix, lil_matrix
from sklearn.metrics.pairwise import euclidean_distances

import numpy as np
import scipy as sp
import random
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets

from sklearn.metrics import confusion_matrix, f1_score, SCORERS
from skll.metrics import kappa


n_neighbors = 15
k=3
reigns = ["arch", "bact", "euk", "virus"]

data = []
rows = []
cols = []
targ = []
contig_to_row = {} # (reign, contig_id) -> row_id
rowid = 0
for i in range(len(reigns)):
    with open("%s_cont_numeric_k%d.mv"% (reigns[i], k), "r") as mvfh:
        contig=""
        for line in mvfh:
            line = line.split()
            contigid = int(line[0].split("_")[0])
            if (i, contigid) not in contig_to_row:
                contig_to_row[(i, contigid)] = rowid
                targ.append(i)
                rowid+=1
                
            rows.append(contig_to_row[(i, contigid)])
            cols.append(int(line[2], 4)) # base 4 to 10
            data.append(int(line[3]))
    
csr = csr_matrix((data, (rows, cols)), shape=(len(contig_to_row), 4**k))

training = random.sample(range(csr.shape[0]), int(0.8*csr.shape[0]))
testing = list(set(range(csr.shape[0])).difference(training))

clf = neighbors.KNeighborsClassifier(n_neighbors, weights="distance")
clf.fit(csr[training,:], [targ[i] for i in training])
classes=clf.predict(csr[testing,:])

print(confusion_matrix(classes,[targ[i] for i in testing]))

print(kappa(classes,[targ[i] for i in testing]))




