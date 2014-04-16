# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 13:42:36 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix

n15 = []
for i in range(50):
    n15.append([])
n25 = []
for i in range(50):
    n25.append([])
c15 = 0
with open("knn,n=15-knn,n=25-k5_arch_bact_euk_virus_2.csv", "r") as fh:
    for line in fh:
        line = line.split()
        if "15" in line[1]:
            n15[int(line[0])].append(line)
            if int(line[3]) == int(line[4]):
                c15 += 1
        else:
            n25[int(line[0])].append(line)


print "### n=15 ###"

for classifier in n15:
    classes = []
    predicted = []
    for contig in classifier:
#        print contig[1], contig[3], contig[4]
        classes.append(int(contig[3]))
        predicted.append(int(contig[4]))

    print "kappa =", kappa(classes,predicted)
#    print c15, "/", len(n15)
#    print confusion_matrix(classes,predicted)

print "### n=25 ###"


for classifier in n25:
    classes = []
    predicted = []
    for contig in classifier:
#        print contig[1], contig[3], contig[4]
        classes.append(int(contig[3]))
        predicted.append(int(contig[4]))

    print "kappa =", kappa(classes,predicted)

print "kappa =", kappa(classes,predicted)