# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, sata_filename):
    from sklearn import neighbors
    
    n = 15
    
    if options:
        for option in options:
            exec option
    
    clf = neighbors.KNeighborsClassifier(n, weights="distance")
    clf.fit(training_matrix, training_targets)
    #print "predicting"
    return clf.predict(testing_matrix)
