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

"""
try:
   import cPickle as pickle
except:
   import pickle
"""
import pickle

if __name__ == '__main__':

    processes = 5
    samplings = 5
    tasks_per_process = [samplings/processes]*processes
    for i in range(samplings%processes):
        tasks_per_process[i] += 1
    
    kvals = range(3,7)
    #kvals = [3]
    #kvals = [4]
    #reignset = [["arch", "bact", "euk", "virus"]]
    reignset = [["arch", "bact", "euk", "virusfiltered"]]
    #reignset = [["all"]]#, ["arch", "bact", "euk", "virusfiltered"]]
    
    git_commit = subprocess.check_output('cd ~/git && git log --format="%H" | head -n1', shell=True).strip()
    resampling_path = os.path.dirname(os.path.abspath("%s_cont_numeric_k%d.mv"% (reignset[0][0], kvals[0])))
    
    uid=0   
    for k in kvals:
        for reigns in reignset:
            print "loop: k=%d" %k
            data = []
            rows = []
            cols = []
            targ = []
            tag  = []
            length = []
            contig_to_row = {} # (reign, contig_id) -> row_id
            rowid = 0
            print "reading file..."
            for i in range(len(reigns)):
                with open("%s_cont_numeric_k%d.mv"% (reigns[i], k), "r") as mvfh:
                    contig=""
                    for line in mvfh:
                        line = line.split()
                        #contigid = int(line[0].split("_")[0])
                        contigid = line[0]
                        if (i, contigid) not in contig_to_row:
                            contig_to_row[(i, contigid)] = rowid
                            targ.append(i)
                            tag.append(line[0])
                            length.append(int(line[1]))
                            rowid+=1
                            
                        rows.append(contig_to_row[(i, contigid)])
                        cols.append(int(line[2], 4)) # base 4 to 10
                        data.append(float(line[3]))
            print "creating the matrix..."
            matrix = csr_matrix((data, (rows, cols)), shape=(len(contig_to_row), 4**k))
            
            with open("matrix_k%d_%s.bdat" % (k, "_".join(reigns)), "w") as fh:
                #pickle.dump([tag, length, targ, matrix.todense()],fh)
                pickle.dump([tag, length, targ, matrix],fh)
            uid +=1
print "Done."
