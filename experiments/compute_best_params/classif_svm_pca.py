# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:58:48 2014

@author: lschmitt
"""
import multiprocessing, sys

class predicter(multiprocessing.Process):
    """
    Hello, I'm a process class, and I'm in charge of predicting a subset 
    of contigs because I can. Yeah.
    """
    #TODO: Make a decent description
    
    def __init__(self, result_queue, model, submat_id, submatrix):
        # base class initialization
        multiprocessing.Process.__init__(self)
        # job management stuff
        self.result_queue = result_queue
        self.kill_received = False
        self.model = model
        self.submat_id = submat_id
        self.submatrix = submatrix

    def run(self):
		"""
		"""            
		classes_subset = self.model.predict(self.submatrix)

		self.result_queue.put((self.submat_id, classes_subset))
		sys.stdout.write("Done with worker %d\n" % (self.submat_id))

def split_predict(testing_matrix, rbf_svc, cores):
	from numpy import array_split, concatenate
	# Let's prepare the queues and subprocesses
	results_queue = multiprocessing.Queue()

	# Split matrix
	sections = array_split(testing_matrix, cores)
	results = [None]*len(sections)

	# Spawn workers
	print "launching subprocesses..."
	for i in range(cores):
		worker = predicter(results_queue, rbf_svc, i, sections[i])
		worker.start()

	# Collect the results off the queue
	print "Collecting results..."
	for i in range(cores): 
		submat_id,classes_subset  = results_queue.get()
		if len(classes_subset) > 0:
			results[submat_id] = classes_subset

	# Merge results
	return concatenate(results)

def learn_and_classify(training_matrix, training_targets, testing_matrix, options, clffile, cores=1):
    from sklearn import svm
    from sklearn.decomposition import PCA
    import pickle, datetime, multiprocessing, Queue

    # Default values for the options
    g = 0.7
    C = 1.0
    n = 100
    
    # If options provided, replace the default values
    if options:
        for option in options:
            exec option
    
    # Dimension reduction
    pca = PCA(n_components=n)
    print "[%s] fit & transform the training matrix" % datetime.datetime.now()
    pca.fit_transform(training_matrix, training_targets)
    print "[%s] transform the testing matrix" % datetime.datetime.now()
    pca.transform(testing_matrix)
    
    # SVM fitting
    print "[%s] learning" % datetime.datetime.now()
    rbf_svc = svm.SVC(kernel='rbf', gamma=g, C=C).fit(training_matrix, training_targets)
    
    # Saving model
    print "[%s] saving model" % datetime.datetime.now()
    with open(clffile, 'w') as fh:
		pickle.dump((pca, rbf_svc), fh)
    
    #print "predicting"
    print "[%s] classifying" % datetime.datetime.now()
    
    return split_predict(testing_matrix, rbf_svc, cores)

def classify_only(pca, rbf_svc, testing_matrix, cores=1):
	from sklearn import svm
	from sklearn.decomposition import PCA
	import datetime, multiprocessing, Queue
	
	# reducing testing matrix
	pca.transform(testing_matrix)
	#print "predicting"
	return split_predict(testing_matrix, rbf_svc, cores)
