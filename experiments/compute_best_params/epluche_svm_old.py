# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:00:47 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

gamma = [90]
csvfile = "rbf_svc,g=90-k3_arch_bact_euk_virus_0.csv"

for g in gamma:
    batches = []
    for i in range(50):
        batches.append([])
    
    with open(csvfile, "r") as fh:
        for line in fh:
            line = line.split()
            if str(g) in line[1]:
                batches[int(line[0])].append(line)
    
    print "### gamma=%d ###"%g
    moy = 0
    
    for classifier in batches:
        classes = []
        predicted = []
        for contig in classifier:
    #        print contig[1], contig[3], contig[4]
            classes.append(int(contig[3]))
            predicted.append(int(contig[4]))
    
        print "kappa =", kappa(classes,predicted)
        moy += kappa(classes,predicted)
        print confusion_matrix(classes,predicted)
    print "moy = ", moy/50


