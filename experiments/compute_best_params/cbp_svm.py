print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
#from sklearn.datasets import load_iris
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV

import sys, pickle

##############################################################################
# Load and prepare data set
#
# dataset for grid search
#iris = load_iris()
#X = iris.data
#Y = iris.target

testfile = sys.argv[1] # first argument: training data for a specific kmer length (.bdat file)
print "loading testing matrix"
with open(testfile, "r") as fh:
    data = pickle.load(fh)
    testdata = {}
    testdata["tag"] = data[0]
    testdata["length"] = data[1]
    testdata["targ"] = data[2]
    testdata["matrix"] = data[3].todense()
X = testdata["matrix"]
Y = testdata["targ"]

# dataset for decision function visualization
#X_2d = X[:, :2]
#X_2d = X_2d[Y > 0]
#Y_2d = Y[Y > 0]
#Y_2d -= 1

# It is usually a good idea to scale the data for SVM training.
# We are cheating a bit in this example in scaling all of the data,
# instead of fitting the transformation on the training set and
# just applying it on the test set.

scaler = StandardScaler()

X = scaler.fit_transform(X)
#X_2d = scaler.fit_transform(X_2d)

##############################################################################
# Train classifier
#
# For an initial search, a logarithmic grid with basis
# 10 is often helpful. Using a basis of 2, a finer
# tuning can be achieved but at a much higher cost.

C_range = 10.0 ** np.arange(-2, 9)
gamma_range = 10.0 ** np.arange(-5, 4)
param_grid = dict(gamma=gamma_range, C=C_range)
cv = StratifiedKFold(y=Y, n_folds=4)
grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv, n_jobs=20)
grid.fit(X, Y)

print("The best classifier is: ", grid.best_estimator_)

# Now we need to fit a classifier for all parameters in the 2d version
# (we use a smaller set of parameters here because it takes a while to train)
#C_2d_range = [1, 1e2, 1e4]
#gamma_2d_range = [1e-1, 1, 1e1]
#classifiers = []
#for C in C_2d_range:
    #for gamma in gamma_2d_range:
        #clf = SVC(C=C, gamma=gamma)
        #clf.fit(X_2d, Y_2d)
        #classifiers.append((C, gamma, clf))
