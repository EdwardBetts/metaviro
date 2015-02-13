# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:12:18 2014

@author: lschmitt
"""

import sys, os.path
#try:
   #import cPickle as pickle
#except:
   #import pickle
import pickle, datetime
from collections import Counter

trainfile = sys.argv[1] # first argument: training data for a specific kmer length (.bdat file)
testfile = sys.argv[2] # second argument: testing data for a specific kmer length (.bdat file)
clfstring = sys.argv[3] # third argument: classifier name as in classif_<name>.py with options
# example with options: rbf_svc,gamma=10,C=1
# example with options: svm_pca,gamma=10,C=1,n=131
# (removed the possibility to run more than one to avoid long as F file names)
try:
	cores = int(sys.argv[4])
except:
	cores = 10
print "cores used: ", cores

dirpath = os.path.split(testfile)[0]
trainfilename = os.path.split(trainfile)[1]
outfile = os.path.splitext(trainfilename)[0]
if outfile.startswith('matrix_'):
    outfile = outfile[len('matrix_'):]
outfile = clfstring + "-" + outfile + ".csv"
outfile = os.path.join(dirpath, outfile)

#clffile = os.path.splitext(trainfile)[0] + ".clf"
clffile = os.path.splitext(trainfile)[0]
if clffile.startswith('matrix_'):
    clffile = clffile[len('matrix_'):]
clffile = clfstring + "-" + clffile + ".clf"
clffile_exists = os.path.isfile(clffile)

print "~~~ %s VS %s ~~~" %(trainfile, testfile)

if clffile_exists:
	print "loading pre-existing classifier model"
	with open(clffile, "r") as fh:
		pca, rbf_svc = pickle.load(fh)
else:
	print "loading training matrix"
	with open(trainfile, "r") as fh:
		data = pickle.load(fh)
		traindata = {}
		traindata["tag"] = data[0]
		traindata["length"] = data[1]
		traindata["targ"] = data[2]
		traindata["matrix"] = data[3].todense()

print "loading testing matrix"
with open(testfile, "r") as fh:
    data = pickle.load(fh)
    testdata = {}
    testdata["tag"] = data[0]
    testdata["length"] = data[1]
    testdata["targ"] = data[2]
    testdata["matrix"] = data[3].todense()

#print Counter(targ)

# normalize data (only once ffs)
print "Normalizing data"
if not clffile_exists:
	for i in range(len(traindata["length"])):
		traindata["matrix"][i,:]=traindata["matrix"][i,:]/traindata["length"][i]

for i in range(len(testdata["length"])):
    testdata["matrix"][i,:]=testdata["matrix"][i,:]/testdata["length"][i]

#with open("filtered_contigs.bdat", "r") as fh:
#    filter_dict = pickle.load(fh)
#clfs = []
with open(outfile, 'w') as outfh:
	classifier = clfstring.split(',')
	print "[%s] Running " % datetime.datetime.now(), clfstring
	options = None
	if len(classifier) > 1:
		options = classifier[1:]
	if clffile_exists:
		exec "from classif_%s import classify_only as classify" % classifier[0]
		classes = classify(pca, rbf_svc, testdata["matrix"], cores)
	else:
		exec "from classif_%s import learn_and_classify as classify" % classifier[0]
		classes = classify(traindata["matrix"], traindata["targ"], testdata["matrix"], options, clffile, cores)
	print "[%s] Saving results" % datetime.datetime.now()
	for i in range(len(testdata["tag"])):
		#outfh.write( '\t'.join([str(batch_id), clfstring, testing_tags[i], str(testing_targets[i]), str(classes[i]), str(filter_dict[testing_tags[i]])]) + '\n')
		outfh.write( '\t'.join([clfstring, testdata["tag"][i], str(testdata["targ"][i]), str(classes[i])]) + '\n')
			
print "[%s] Done." % datetime.datetime.now()
