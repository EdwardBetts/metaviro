# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options):
    from sklearn import svm
    
    gamma = 0.7
    C = 1.0
    
    if options:
        for option in options:
            exec option
    
    #print "fitting"
    rbf_svc = svm.SVC(kernel='rbf', gamma=gamma, C=C).fit(training_matrix, training_targets)
    
    #print "predicting"
    return rbf_svc.predict(testing_matrix)
