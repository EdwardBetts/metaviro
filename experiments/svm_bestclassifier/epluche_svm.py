# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:00:47 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix
import numpy as np
import os, sys
try:
   import cPickle as pickle
except:
   import pickle
   

g = 90
csvfile = "rbf_svc,g=90-k3_arch_bact_euk_virus_0.csv"
clffile = os.path.splitext(csvfile)[0] + ".clf"
#datafile = "sample_matrix_k3_inra_0.bdat"
datafile = "sample_matrix_k3_arch_bact_euk_virus_0.bdat"

kappas = []

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

    kappas.append(kappa(classes,predicted))

print "tous classifieurs:"
print "kappa moyen = %f" % np.mean(kappas)
print "écart-type = %f" % np.std(kappas)

best_batchid = kappas.index(max(kappas))
print "meilleur kappa = %f" % kappas[best_batchid]

print "loading best classifier..."
with open(clffile, "r") as clffh:
    bestclf = pickle.load(clffh)[best_batchid]
    
# Re-run the prediction with the best classifier

print "loading data..."
with open(datafile, "r") as fh:
    tag, length, targ, matrix, training_indices = pickle.load(fh)

# normalize data (only once ffs)
for i in range(len(length)):
    matrix[i,:]=matrix[i,:]/length[i]

print "running best classifier..."

predicted_classes = bestclf.predict(matrix)
kp = kappa(targ, predicted_classes)

correctly_assigned = 0.0
bases_assigned = 0.0
for i in range(len(targ)):
    if targ[i] == predicted_classes[i]:
        correctly_assigned +=1
        bases_assigned += length[i]

print "Other independant resampling:"
print "kappa = %f" % kp
print "% correct:", correctly_assigned / len(targ) * 100
print "% bases correct:", bases_assigned / np.sum(length) *100
print confusion_matrix(targ,predicted_classes)


#print "running best classifier on viruses only..."
#
#matrixv = matrix[np.array(targ)==3]
#targv = np.array(targ)[np.array(targ)==3]
#
#predicted_classes = bestclf.predict(matrixv)
#kp = kappa(targv, predicted_classes)
#
#correctly_assigned = 0.0
#for i in range(len(targv)):
#    if targv[i] == predicted_classes[i]:
#        correctly_assigned +=1
#
#print "Other independant resampling:"
#print "kappa = %f" % kp
#print "% correct:", correctly_assigned / len(targv) * 100





