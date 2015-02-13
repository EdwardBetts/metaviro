# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:00:47 2014

@author: lschmitt
"""

from skll.metrics import kappa
from sklearn.metrics import confusion_matrix
import numpy as np
import matplotlib.pyplot as plt
import os, sys
try:
   import cPickle as pickle
except:
   import pickle
   

csvfile = sys.argv[1]
#datafile = sys.argv[2]
fileprefix = os.path.splitext(csvfile)[0]

classes = []
predicted = []
with open(csvfile, "r") as fh:
    for line in fh:
		line = line.split()
		classes.append(int(line[2]))
		predicted.append(int(line[3]))

try:
	kp = kappa(classes, predicted)
	print "kappa = ", kp
except:
	kp = -1.0

categories = os.path.splitext(csvfile)[0].split("-k")[1].split("_")[1:]
print categories
confmat = np.zeros((len(categories),len(categories)))
for i in range(len(classes)):
	confmat[classes[i],predicted[i]] += 1
print confmat
print type(confmat)


conf_arr = confmat
norm_conf = []
for i in conf_arr:
	a = 0
	tmp_arr = []
	a = sum(i, 0)
	for j in i:
		tmp_arr.append(float(j)/float(a))
	norm_conf.append(tmp_arr)

fig = plt.figure()
plt.clf()
ax = fig.add_subplot(111)
ax.set_aspect(1)
res = ax.imshow(np.array(norm_conf), cmap=plt.cm.summer, interpolation='nearest', vmin=0, vmax=1)

width = len(conf_arr)
height = len(conf_arr[0])

for x in xrange(width):
	for y in xrange(height):
		ax.annotate(str(conf_arr[x][y]), xy=(y, x),
					horizontalalignment='center',
					verticalalignment='center')

cb = fig.colorbar(res)
plt.xticks(range(width), categories)
plt.yticks(range(height), categories)
plt.xlabel('Assigned categories')
plt.ylabel('Real categories')
plt.title("%s\nkappa=%f" % (fileprefix, kp))
#plt.show()
plt.savefig("%s.png" % fileprefix, format='png')




