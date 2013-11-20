# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:39:53 2013

Relative Oligonucleotide Frequency from pseudocounts vectors for the 
Jensen-Shannon divergence

@author: lschmitt
"""


import argparse, csv, multiprocessing, Queue, sys, numpy, os
from collections import defaultdict
from scipy.stats import pearsonr

import warnings

#def fxn():
#    warnings.warn("runtime", RuntimeWarning)
#
#with warnings.catch_warnings(record=True) as w:
#    # Cause all warnings to always be triggered.
#    warnings.simplefilter("error")
#    # Trigger a warning.
#    fxn()
#    # Verify some things
#    assert len(w) == 1
#    assert issubclass(w[-1].category, RuntimeWarning)
#    assert "runtime" in str(w[-1].message)

class basioConverter(multiprocessing.Process):
    """
    Hello, I'm a process class, and I'm in charge of processing contigs 
    because I can. Can't think of a better description for now.
    """
    #TODO: Make a decent description
    
    def __init__(self, work_queue, result_queue, number):
        # base class initialization
        multiprocessing.Process.__init__(self)
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.number = number
        # reference info
#        self.ref_vector = ref_vector

    def run(self):
        """
        """
        loop_count = 0
        while not self.kill_received:
            # get a contig's results to process
            try:
                contig_info = self.work_queue.get(True, 5)
            except Queue.Empty:
                sys.stdout.write("Queue empty for worker %d\n" % (self.number))
                break
            if contig_info == 0:
                break
            
            # Just because. Could be useful.
            loop_count = loop_count +1
            
            # getting ref_vector:
            ref_vector = contig_info[2]
            
            # creating query vector 
#            query_vector = [1.0]*len(ref_vector)
            query_vector = [0.0]*len(ref_vector)
            for line in contig_info[1]:
                query_vector[int(line[2])] = float(line[3]) #+1

            
            # make frequency vectors
            freq_vectors = [[],[]]
            query_total = numpy.sum(query_vector)
            ref_total = numpy.sum(ref_vector)
            for i in range(len(ref_vector)):
                if query_vector[i] != 0.0:
                    freq_vectors[0].append(query_vector[i]/query_total)
                    freq_vectors[1].append(ref_vector[i]/ref_total)
#            print freq_vectors
#            with open('%s_%s.vectors' % (contig_info[3], contig_info[0]), 'w') as vfh:
#                vfh.write(' '.join(map(str,query_vector)) + "\n")
#                vfh.write(' '.join(map(str,ref_vector)) + "\n")
#                vfh.write(' '.join(map(str,freq_vectors[0])) + "\n")
#                vfh.write(' '.join(map(str,freq_vectors[1])) + "\n")
#            
            kldiv = 0
            for i in range(len(freq_vectors[0])):
                try:
                    with warnings.catch_warnings(record=True) as w:
                        # Cause all warnings to always be triggered as exceptions.
                        warnings.simplefilter("error")
                        kldiv = kldiv + ( freq_vectors[0][i] * numpy.log2( freq_vectors[0][i] / freq_vectors[1][i] ) )
                except:
                    sys.stderr.write( str(freq_vectors[0][i]) + ' ' + str(freq_vectors[1][i]) + "\n")
#                    sys.stderr.write( str(query_vector[i]) + ' ' + str(i) + "\n")
            
            kldiv2 = 0
            for i in range(len(freq_vectors[0])):
                try:
                    with warnings.catch_warnings(record=True) as w:
                        # Cause all warnings to always be triggered as exceptions.
                        warnings.simplefilter("error")
                        kldiv2 = kldiv2 + ( freq_vectors[1][i] * numpy.log2( freq_vectors[1][i] / freq_vectors[0][i] ) )
                except:
                    sys.stderr.write( str(freq_vectors[0][i]) + ' ' + str(freq_vectors[1][i]) + "\n")
#                    sys.stderr.write( str(query_vector[i]) + ' ' + str(i) + "\n")
            
            self.result_queue.put((contig_info[0], numpy.mean([kldiv, kldiv2])))
            
        # store the results in the results queue once all the contigs are processed
        # please run please run please run please run please run please run 
#        
        sys.stdout.write("Done with worker %d, %d loops\n" % (self.number, loop_count))

if __name__ == '__main__':
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-t", "--cpu", help="Max number of subprocesses", default = "10", type=int)
    argparser.add_argument("query", help=".mvdist query file")
    argparser.add_argument("ref", help=".mvref reference files (accepts multiple files)", nargs='+')
    args = argparser.parse_args()
    
    #TODO: verify if the files exist (or maybe wait for the shit to hit the fan anyway)
    
    # Getting the prefixes (for convenience)
    query_prefix = os.path.splitext(args.query)[0]
    ref_prefixes = [os.path.splitext(reffile)[0] for reffile in args.ref]
    
    
    # Let's prepare the queues and subprocesses
    work_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue()
    all_results = defaultdict(list)
    
    # One loop per reference
    for ref_prefix in ref_prefixes:
        print ref_prefix
        
        ref_segments_converter = None
        ref_vector = None
        
        print "reading reference segments"
        with open('%s.mvseg' % ref_prefix, 'rb') as segfh:
            segfile = csv.reader(segfh, delimiter=' ')
            ref_segments_converter = {rows[1]:int(rows[0]) for rows in segfile}
            ref_id_max = len(ref_segments_converter) - 1
        print "reading reference dist"
        with open('%s.mvref' % ref_prefix, 'rb') as reffh:
            # We only store the counts since we rely on the ids, not on the
            # segments' sequence. Therefore we'll use the arrays indices.
            ref_vector = [-1]*len(ref_segments_converter)
            for line in csv.reader(reffh, delimiter=' '):
                ref_vector[int(line[2])] = float(line[3])
        
        query_contigs = defaultdict(list)
        query_segments = None        
        
        print "reading query segments"
        with open('%s.mvseg' % query_prefix, 'rb') as segfh:
            segfile = csv.reader(segfh, delimiter=' ')
            query_segments = list(zip(*segfile)[1])
        print "reading query dist"
        with open('%s.mvdist' % query_prefix, 'rb') as reffh:
            for line in csv.reader(reffh, delimiter=' '):
                # Let's convert the ids to the reference equivalent
                #TODO: Adapt this to the PATATAS
                line[2] = int(line[2])
                try:
                    if line[2] != ref_segments_converter[query_segments[line[2]]]:
                        line[2] = ref_segments_converter[query_segments[line[2]]]
                except KeyError: # The segments is not referenced in the hash
                    ref_segments_converter[query_segments[line[2]]] = len(ref_segments_converter)
                    ref_vector.append(0.1)
                    if line[2] != ref_segments_converter[query_segments[line[2]]]:
                        line[2] = ref_segments_converter[query_segments[line[2]]]
                query_contigs[line[0]].append(line)
        
        print "loading queue"
        for contig in query_contigs:
            work_queue.put((contig, query_contigs[contig], ref_vector))
        for i in range(args.cpu):
            work_queue.put(0) # To forcibly stop the workers when they're done
        
        # Spawn workers
        print "launching subprocesses..."
        for i in range(args.cpu):
            worker = basioConverter(work_queue, results_queue, i+1)
            worker.start()
        
        print "Collecting results..."
        for i in range(len(query_contigs)): 
            distance = results_queue.get()
            if distance:
    #            print distance
                all_results[distance[0]].append(distance[1])
    
    print "saving to file..."
    maxcounts = [0]*len(ref_prefixes)
    with open('%s.js' % query_prefix, 'w') as outfh:
        outfh.write('contig_id\t')
        sys.stdout.write('contig_id\t')
        outfh.write('\t'.join(ref_prefixes))
        sys.stdout.write('\t'.join(ref_prefixes))
        outfh.write('\tmin\n')
        sys.stdout.write('\tmin\n')
        for contig in all_results:
            outfh.write('%s\t' % contig)
            outfh.write('\t'.join([ str(val) for val in all_results[contig] ]))
            corrmax = all_results[contig].index(min(all_results[contig]))
            outfh.write('\t%d\n' % corrmax)
            maxcounts[corrmax] += 1
        outfh.write('cumulative_min:\t')
        sys.stdout.write('cumulative_min:\t')
        outfh.write('\t'.join([ str(val) for val in maxcounts]))
        sys.stdout.write('\t'.join([ str(val) for val in maxcounts]))
        outfh.write('\t%s\n' % ref_prefixes[maxcounts.index(max(maxcounts))])
        sys.stdout.write('\t%s\n' % ref_prefixes[maxcounts.index(max(maxcounts))])
    
    print "done"
    
