# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:55:15 2013

K-NN

@author: lschmitt
"""

import argparse, csv, multiprocessing, Queue, sys, numpy, os
from collections import defaultdict
#~ from scipy.stats import pearsonr
from scipy.spatial.distance import euclidean

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
            
            # getting ref_dict:
            ref_dict = contig_info[2]
            query_dict = contig_info[1]
            
            query_vector = []
            ref_vector = []
            
            for seg in query_dict.keys():
				try:
					ref_vector.append(float(ref_dict[seg][0]))
					query_vector.append(float(query_dict[seg][0]))
				except KeyError:
					pass
            #~ print "ref", len(ref_vector), "query", len(query_vector)
            
            
            # make frequency vectors
            freq_vectors = [[],[]]
            query_total = numpy.sum(query_vector)
            ref_total = numpy.sum(ref_vector)
            for i in range(len(ref_vector)):
                freq_vectors[0].append(query_vector[i]/query_total)
                freq_vectors[1].append(ref_vector[i]/ref_total)
            
            
            self.result_queue.put((contig_info[0], euclidean(freq_vectors[0],freq_vectors[1])))
            
        # store the results in the results queue once all the contigs are processed
        # please run please run please run please run please run please run 
#        
        sys.stdout.write("Done with worker %d, %d loops\n" % (self.number, loop_count))

if __name__ == '__main__':
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-t", "--cpu", help="Max number of subprocesses", default = "10", type=int)
    argparser.add_argument("-k", "--kval", help="k-mer length", required=True, type=int)
    argparser.add_argument("-o", "--outdir", help="Output dir", default = None)
    argparser.add_argument("query", help=".mvdist query file")
    argparser.add_argument("ref", help=".mvref reference files (accepts multiple files)", nargs='+')
    args = argparser.parse_args()
    
    #TODO: verify if the files exist (or maybe wait for the shit to hit the fan anyway)
    
    # Getting the prefixes (for convenience)
    query_prefix = os.path.splitext(args.query)[0]
    ref_prefixes = [os.path.splitext(reffile)[0] for reffile in args.ref]
    if args.outdir:
		if not os.path.isdir(args.outdir):
			os.makedirs(args.outdir)
		outfile_prefix = os.path.join(args.outdir, os.path.splitext(os.path.basename(os.path.abspath(args.query)))[0])
    else:
		outfile_prefix = query_prefix
    
    # Let's prepare the queues and subprocesses
    work_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue()
    all_results = defaultdict(list)
    
    
    data = []
    rows = []
    cols = []
    targ = []
    contig_to_row = {} # (reign, contig_id) -> row_id
    rowid = 0
	for i in range(len(reigns)):
		with open("%s_cont_numeric_k%d.mv"% (reigns[i], k), "r") as mvfh:
			contig=""
			for line in mvfh:
				line = line.split()
				contigid = int(line[0].split("_")[0])
				if (i, contigid) not in contig_to_row:
					contig_to_row[(i, contigid)] = rowid
					targ.append(i)
					rowid+=1
					
				rows.append(contig_to_row[(i, contigid)])
				cols.append(int(line[2], 4)) # base 4 to 10
				data.append(int(line[3]))
				
	csr = csr_matrix((data, (rows, cols)), shape=(len(contig_to_row), 4**k))
    
    # One loop per reference
    for ref_prefix in ref_prefixes:    
        
        ref_dict = {}
            
        print "reading reference dist"
        with open('%s.mv' % ref_prefix, 'rb') as reffh:
            # We only store the counts since we rely on the ids, not on the
            # segments' sequence. Therefore we'll use the arrays indices.
            for line in csv.reader(reffh, delimiter=' '):
                ref_dict[line[2]] = [x + y for x, y in zip([int(line[3]),int(line[1])], ref_dict.setdefault(line[2], [0,0]))]
        
        ref_segments = ref_dict.keys()
        
        
        query_contigs = defaultdict(dict)
        
        print "reading query dist"
        with open('%s.mv' % query_prefix, 'rb') as reffh:
            for line in csv.reader(reffh, delimiter=' '):
                query_contigs[line[0]][line[2]] = [line[3],line[1]]
        
        
        print "loading queue"
        for contig in query_contigs:
            work_queue.put((contig, query_contigs[contig], ref_dict))
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
    with open('%s.eucl' % outfile_prefix, 'w') as outfh:
        outfh.write('contig_id\t')
        outfh.write('\t'.join(ref_prefixes))
        outfh.write('\tmax\n')
        for contig in all_results:
            outfh.write('%s\t' % contig)
            outfh.write('\t'.join([ str(val) for val in all_results[contig] ]))
            distmin = all_results[contig].index(min(all_results[contig]))
            outfh.write('\t%d\n' % distmin)
            maxcounts[distmin] += 1
        outfh.write('cumulative_max:\t')
        outfh.write('\t'.join([ str(val) for val in maxcounts]))
        outfh.write('\t%s\n' % ref_prefixes[maxcounts.index(max(maxcounts))])
    
    print "done"
    
