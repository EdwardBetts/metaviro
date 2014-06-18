# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, idbatch):
    from sklearn import svm
    
    g = 0.7
    C = 1.0
    
    if options:
        for option in options:
            exec option
    
    #print "fitting"
    rbf_svc = svm.SVC(kernel='rbf', gamma=g, C=C).fit(training_matrix, training_targets)
    
    #print "predicting"
    return rbf_svc.predict(testing_matrix), rbf_svc
