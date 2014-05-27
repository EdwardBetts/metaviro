# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 16:34:03 2014

@author: lschmitt
"""

def classify(training_matrix, training_targets, testing_matrix, options, data_filename):
    from sklearn.naive_bayes import MultinomialNB
    
    # NO OPTIONS
    
    # Spawn classifier and predict
    clf = MultinomialNB().fit(training_matrix, training_targets)
    return clf.predict(testing_matrix)