# -*- coding: utf-8 -*-
"""
Created on Fri May  2 14:07:47 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

nvals = [100,200,300,400,500]
csvfile = "rforest,n=100,crit='entropy'-rforest,n=200,crit='entropy'-rforest,n=300,crit='entropy'-rforest,n=400,crit='entropy'-rforest,n=500,crit='entropy'-k6_arch_bact_euk_virus_3.csv"

for n in nvals:
    print "### n=%s ###"%str(n)
    batches = []
    for i in range(50):
        batches.append([])
    c15 = 0
    with open(csvfile, "r") as fh:
        for line in fh:
            line = line.split()
            if str(n) in line[1]:
                batches[int(line[0])].append(line)
                if int(line[3]) == int(line[4]):
                    c15 += 1

    moy15 = 0
    for classifier in batches:
        classes = []
        predicted = []
        for contig in classifier:
    #        print contig[1], contig[3], contig[4]
            classes.append(int(contig[3]))
            predicted.append(int(contig[4]))
    
        print "kappa =", kappa(classes,predicted)
        moy15 += kappa(classes,predicted)
        
    #    print c15, "/", len(batches)
    #    print confusion_matrix(classes,predicted)
    
    print "moy = ", moy15/50