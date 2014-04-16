# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:00:47 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

g1 = []
for i in range(50):
    g1.append([])
g5 = []
for i in range(50):
    g5.append([])
g10 = []
for i in range(50):
    g10.append([])

with open("rbf_svc,gamma=1-rbf_svc,gamma=5-rbf_svc,gamma=10-k5_arch_bact_euk_virus_2.csv", "r") as fh:
    for line in fh:
        line = line.split()
        if "10" in line[1]:
            g10[int(line[0])].append(line)
        elif "5" in line[1]:
            g5[int(line[0])].append(line)
        else:
            g1[int(line[0])].append(line)


print "### gamma=1 ###"

for classifier in g1:
    classes = []
    predicted = []
    for contig in classifier:
#        print contig[1], contig[3], contig[4]
        classes.append(int(contig[3]))
        predicted.append(int(contig[4]))

    print "kappa =", kappa(classes,predicted)
#    print c15, "/", len(n15)
    print confusion_matrix(classes,predicted)

print "### gamma=5 ###"

for classifier in g5:
    classes = []
    predicted = []
    for contig in classifier:
#        print contig[1], contig[3], contig[4]
        classes.append(int(contig[3]))
        predicted.append(int(contig[4]))

    print "kappa =", kappa(classes,predicted)
    print confusion_matrix(classes,predicted)


print "### gamma=10 ###"

for classifier in g10:
    classes = []
    predicted = []
    for contig in classifier:
#        print contig[1], contig[3], contig[4]
        classes.append(int(contig[3]))
        predicted.append(int(contig[4]))

    print "kappa =", kappa(classes,predicted)
    print confusion_matrix(classes,predicted)


