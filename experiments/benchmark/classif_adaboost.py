# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 16:31:08 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, data_filename):
    from sklearn.ensemble import AdaBoostClassifier
    
    n = 100    
    
    if options:
        for option in options:
            exec option
    
    # Spawn classifier and predict
    clf = AdaBoostClassifier(n_estimators=n).fit(training_matrix, training_targets)
    return clf.predict(testing_matrix)