# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:33:26 2014

@author: lschmitt
"""


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
