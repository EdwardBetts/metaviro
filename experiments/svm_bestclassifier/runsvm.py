# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:12:18 2014

@author: lschmitt
"""

import sys, os.path
try:
   import cPickle as pickle
except:
   import pickle

filename = sys.argv[1] # first argument: data for a specific kmer length (.bdat file)
classifiers = sys.argv[2:] # other arguments: classifiers names as in classif_<name>.py with options
# example with options: rbf_svc,gamma=10,C=1

outfile = os.path.splitext(filename)[0] + ".csv"

if outfile.startswith('sample_matrix_'):
    outfile = outfile[len('sample_matrix_'):]
outfile = "-".join(classifiers) + "-" + outfile

clffile = os.path.splitext(outfile)[0] + ".clf"


print "~~~ %s ~~~" %filename

print "loading matrix"
with open(filename, "r") as fh:
    tag, length, targ, matrix, training_indices = pickle.load(fh)

# normalize data (only once ffs)
for i in range(len(length)):
    matrix[i,:]=matrix[i,:]/length[i]

#with open("filtered_contigs.bdat", "r") as fh:
#    filter_dict = pickle.load(fh)
clfs = []
with open(outfile, 'w') as outfh:
    for classifierstring in classifiers:
        print "Running ", classifierstring
        classifier = classifierstring.split(',')
        options = None
        if len(classifier) > 1:
            options = classifier[1:]
        exec "from classif_%s import classify as classify" % classifier[0]
        for batch_id in range(len(training_indices)):
            sys.stdout.write("batch %d/%d\r"%(batch_id+1,len(training_indices)))
            sys.stdout.flush()
        
            training = training_indices[batch_id][0]
            testing = training_indices[batch_id][1]
            
            training_matrix = matrix[training,:]
            training_targets = [targ[i] for i in training]
            
            testing_matrix = matrix[testing,:]
            testing_targets = [targ[i] for i in testing]
                        
            testing_tags = [tag[i] for i in testing]
                        
            classes, classifier = classify(training_matrix, training_targets, testing_matrix, options, batch_id)
            clfs.append(classifier)
            for i in range(len(testing_tags)):
#                outfh.write( '\t'.join([str(batch_id), classifierstring, testing_tags[i], str(testing_targets[i]), str(classes[i]), str(filter_dict[testing_tags[i]])]) + '\n')
                outfh.write( '\t'.join([str(batch_id), classifierstring, testing_tags[i], str(testing_targets[i]), str(classes[i])]) + '\n')
            
        print

with open(clffile, 'w') as fh:
    pickle.dump(clfs, fh)

print "Done."