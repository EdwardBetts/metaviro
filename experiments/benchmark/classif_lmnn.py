# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 11:35:05 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, data_filename):
    from sklearn import neighbors
    import numpy as np
    import os
    try:
       import cPickle as pickle
    except:
       import pickle
    
    n = 15
    
    if options:
        for option in options:
            exec option
    
    # Load the pre-calculated metric
    metric_dump = os.path.splitext(data_filename)[0] + ".metric"
    with open(metric_dump, "r") as fh:
        metric = pickle.load(fh)
    
    # Transform both training and testing coordinates with the metric
    new_X_train = metric.transform(np.asarray(training_matrix))
    new_X_test = metric.transform(np.asarray(testing_matrix))
    
    # Classic kNN
    clf = neighbors.KNeighborsClassifier(n, weights="distance")
    clf.fit(new_X_train, training_targets)
    return clf.predict(new_X_test)