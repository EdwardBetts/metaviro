# create matrices from mv file

from __future__ import with_statement
import time, re
from scipy.sparse import coo_matrix, csr_matrix, dok_matrix, lil_matrix
from sklearn.metrics.pairwise import euclidean_distances

import numpy as np
import scipy as sp
import random, multiprocessing, Queue, sys, time, subprocess, os.path
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets

from sklearn.metrics import confusion_matrix, f1_score, SCORERS
from skll.metrics import kappa



class KNNrunner(multiprocessing.Process):
    """
    Hello, I'm a process class, and I'm in charge of processing contigs 
    because I can. Can't think of a better description for now.
    """
    #TODO: Make a decent description
    
    def __init__(self, matrix, result_queue, ntasks, trainingratio, info_log):
        # base class initialization
        multiprocessing.Process.__init__(self)
        # job management stuff
        self.matrix = matrix
        self.result_queue = result_queue
        self.kill_received = False
        self.ntasks = ntasks
        self.trainingratio = trainingratio
        self.info_log = info_log
        
        # reference info
#        self.ref_vector = ref_vector

    def run(self):
        """
        run Forrest
        """
        
        for loopcount in range(self.ntasks):
            seed = time.time()
            resultsline = []
            # do stuff
            training = random.sample(range(self.matrix.shape[0]), int(self.trainingratio*self.matrix.shape[0]))
            testing = list(set(range(self.matrix.shape[0])).difference(training))
            
            print  self.matrix[training,:].shape, len([targ[i] for i in training]), self.matrix[testing,:].shape
            
            clf = neighbors.KNeighborsClassifier(n_neighbors, weights="uniform")
#            clf.fit(self.matrix[training,:].todense(), [targ[i] for i in training])
            clf.fit(self.matrix[training,:], [targ[i] for i in training])
#            classes=clf.predict(self.matrix[testing,:].todense())
            classes=clf.predict(self.matrix[testing,:])
            
#            print(confusion_matrix(classes,[targ[i] for i in testing]))
            
#            print(kappa(classes,[targ[i] for i in testing]))

            resultsline = []
            resultsline = resultsline+info_log
            resultsline.append(seed)
            resultsline.append(kappa(classes,[targ[i] for i in testing]))
            conf = []
            for row in confusion_matrix(classes,[targ[i] for i in testing]):
                conf = conf+list(row)
            resultsline = resultsline+conf
            
            self.result_queue.put([str(item) for item in resultsline])
            
        # store the results in the results queue once all the contigs are processed
        # please run please run please run please run please run please run 
#        
        sys.stdout.write("Done with worker for %d tasks: %d loops done\n" % (self.ntasks, loopcount+1))


if __name__ == '__main__':

    processes = 5
    samplings = 5
    tasks_per_process = [samplings/processes]*processes
    for i in range(samplings%processes):
        tasks_per_process[i] += 1
    
    n_neighbors = 10
    trainingratio = 0.8
#    kvals = range(3,8)
    kvals = [7]
    reigns = ["arch", "bact", "euk", "virus"]
    git_commit = subprocess.check_output('cd ~/git && git log --format="%H" | head -n1', shell=True).strip()
    resampling_path = os.path.dirname(os.path.abspath("balanced_%d_mers.tsv"% kvals[0]))
    
    for k in kvals:
        print "loop: k=%d" %k 
        data = []
        rows = []
        cols = []
        targ = []
        
        print "reading file..."
        with open("balanced_%d_mers.tsv"%k, "r") as mvfh:
            mvfh.readline()
            mvfh.readline()
            rowid = 0
            contig=""
            for line in mvfh:
                line = line.split()
                for reign in reigns:
                    if reign in line[1]:
                        targ.append(reigns.index(reign))
                        break
                for i in range(len(line[4:])):
                    rows.append(rowid)
                    cols.append(i) # base 4 to 10
                    data.append(int(line[4+i]))
                rowid +=1
        print "creating the matrix..."
        csr = csr_matrix((data, (rows, cols)), shape=(rowid, 4**k))
        
        
        
        
        # Let's prepare the queues and subprocesses
        results_queue = multiprocessing.Queue()
        all_results = []
        
        info_log = [time.time(), sys.argv[0], git_commit, k, trainingratio, "knn_%s"%n_neighbors, ";".join(reigns), resampling_path]
    #    info_log = [time.time(), sys.argv[0], git_commit, k, trainingratio, n_neighbors, "knn_%s"%n_neighbors, repr(reigns), resampling_path]
        
        # Spawn workers
        print "launching subprocesses..."
        for ntasks in tasks_per_process:
            worker = KNNrunner(csr, results_queue, ntasks, trainingratio, info_log)
            worker.start()
        
        print "Collecting results..."
        for i in range(samplings): 
            all_results.append(results_queue.get())
            print len(all_results), '/', samplings
        
        with open("results_for_k%s_%s.csv"%(k,reigns[3]), "w") as resultsfh:
            for result in all_results:
                resultsfh.write( "\t".join(result) + "\n")

