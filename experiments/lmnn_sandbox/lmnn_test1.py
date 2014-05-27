# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 16:23:16 2014

@author: lschmitt
"""

import numpy as np
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets
from skll.metrics import kappa
import random
from metric_learn.lmnn import LMNN
n_neighbors = 15
weights = "distance"
iris = datasets.load_iris()

print type(iris.data)

X = iris.data[:, :2]  # we only take the first two features. We could
                      # avoid this ugly slicing by using a two-dim dataset
y = iris.target

training = random.sample(range(150), 100)
testing = list(set(range(150)).difference(training))

X_train = iris.data[training,:]
y_train = [y[i] for i in training]
X_test = iris.data[testing,:]
y_test = [y[i] for i in testing]

clf3 = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
clf3.fit(X_train, y_train)

print "kappa =", kappa(y_test,clf3.predict(X_test))

metric = LMNN(X_train, y_train)
metric.fit()
new_X_train = metric.transform()
new_X_test = metric.transform(X_test)
clf4 = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
clf4.fit(new_X_train, y_train)

print "kappa =", kappa(y_test,clf4.predict(new_X_test))


#
## Create color maps
#cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
#cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])
#
## meshgrid for both plots
#h = .02  # step size in the mesh
#x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
#y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
#xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
#
#
#metric = LMNN(X, y)
#metric.fit()
#new_X = metric.transform()
#
#
## color the mesh for regular knn
#clf = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
#clf.fit(X, y)
#Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
#
## Put the result into a color plot
#Z = Z.reshape(xx.shape)
## and fucking plot
#pl.subplot(1,2,1)
#pl.pcolormesh(xx, yy, Z, cmap=cmap_light)
#pl.scatter(X[:, 0], X[:, 1], c=y,cmap=cmap_bold)
#pl.xlim(xx.min(), xx.max())
#pl.ylim(yy.min(), yy.max())
#pl.title("kNN (k = %i, weights = '%s')" % (n_neighbors, weights))
#
## color the mesh for lmnn knn
#clf2 = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
#clf2.fit(new_X, y)
#Z2 = clf2.predict(metric.transform(np.c_[xx.ravel(), yy.ravel()]))
#
## Put the result into another color plot
#Z2 = Z2.reshape(xx.shape)
## and fucking plot
#pl.subplot(1,2,2)
#pl.pcolormesh(xx, yy, Z2, cmap=cmap_light)
#pl.scatter(X[:, 0], X[:, 1], c=y,cmap=cmap_bold)
#pl.xlim(xx.min(), xx.max())
#pl.ylim(yy.min(), yy.max())
#pl.title("LMNN + kNN (k = %i, weights = '%s')" % (n_neighbors, weights))
#
#pl.show()


