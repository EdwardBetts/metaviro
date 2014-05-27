# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 16:21:27 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, data_filename):
    from sklearn.ensemble import RandomForestClassifier
    
    n=10
    crit='gini'
    
    if options:
        for option in options:
            exec option
    
    # Spawn classifier and predict
    clf = RandomForestClassifier(n_estimators=n, criterion=crit).fit(training_matrix, training_targets)
    return clf.predict(testing_matrix)