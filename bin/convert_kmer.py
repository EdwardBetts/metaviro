# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 18:40:19 2013

@author: Louise-Amélie Schmitt
"""

import multiprocessing, Queue, sys, argparse, os, subprocess, re, pprint
import oak
from io import StringIO


class basioConverter(multiprocessing.Process):
    """
    Hello, I'm a process class, and I'm in charge of processing a subset of 
    segments because I can. Can't think of a better description for now.
    """
    #TODO: Make a decent description
    
    def __init__(self, work_queue, result_queue, fastafile, number):
        # base class initialization
        multiprocessing.Process.__init__(self)
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.number = number
        # You must return here with a shrubbery... or else you will never pass through this wood...   alive.
        self.shrub = oak.plant()
        self.fastafile = fastafile

    def run(self):
        """
        """
        loop_count = 0
        while not self.kill_received:
            # get a contig's results to process
            try:
                sgmblock = self.work_queue.get(True, 10)
            except Queue.Empty:
                sys.stdout.write("Queue empty for worker %d\n" % (self.number))
                break
        
            # Just because. Could be useful.
            loop_count = loop_count +1
            
            # Get the sequence of the contig with samtools and the header
            cmd = "module load samtools && samtools faidx %s \"%s\"" %(fastafile, sgmblock["header"])
#            sys.stderr.write(cmd + "\n")
            record = subprocess.check_output(cmd, shell=True).split()
            contig = ''.join(record[1:])
            if not re.match('^[ATGC]*$', str(contig).upper()):
#                sys.stderr.write("\ndafuq did I just read?\n")
#                sys.stderr.write("skipping contig: %s\n\n" % contig)
                continue
            # add the new sequence to the worker's shrub
            oak.start_new_sequence(self.shrub, len(contig), sgmblock["header"])
            # and the segments one by one
            for kmerline in sgmblock["kmers"]:
                kmerline = kmerline.split()
                oak.add_segment(self.shrub, kmerline[0], int(kmerline[1]) )
            # and fill the blanks if any (not necessarily relevant for kmers though)
            oak.complete_tree(self.shrub)
            
        # store the results in the results queue once all the contigs are processed
        # please run please run please run please run please run please run 
        self.result_queue.put(self.shrub)
        sys.stdout.write("Done with worker %d, %d loops\n" % (self.number, loop_count))


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-r", "--isreference", help="Create a reference vector for a defined family", action="store_true")
    argparser.add_argument("-d", "--dumpmode", help="how to store the data on disk", choices=["mv","pickle", "both"], default="mv")
    argparser.add_argument("-t", "--cpu", help="Max number of subprocesses", default = "1", type=int)
    argparser.add_argument("-f", "--fastafile", help="Multifasta file path (default: deduced from frags file name)", default=None)
    argparser.add_argument("fragsfile", help=".frags file path", default="./")
    args = argparser.parse_args()
    
    fastafile = args.fastafile
    if not args.fastafile:
        fastafile = re.sub(r'_[a-zA-Z]+[0-9\.]+\.[a-z]+$', '.fasta', args.fragsfile)
        pass 
    
    # Also, the fasta file should already have a fasta index so we'll skip that conditional step for now
    
    # Let's prepare the queues and subprocesses
    work_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue()
    queue_length = 0
    
    # Read the .sgm file linearly and push the contents in the input queue
    print "Loading queue..."
    fragsfile = open(args.fragsfile, "r")
    sgmblock = None
    for line in fragsfile.readlines():
        if line.startswith(">"):
            if sgmblock:
                work_queue.put(sgmblock)
                queue_length = queue_length +1
                pass
            sgmblock = {}
            sgmblock["header"] = line[1:].strip()
        elif len(line) != 0:
            if "kmers" in sgmblock.keys():
                sgmblock["kmers"].append(line)
            else:
                sgmblock["kmers"] = [line]
    work_queue.put(sgmblock)
    queue_length = queue_length +1
    
    # Spawn workers
    print "launching subprocesses..."
    for i in range(args.cpu):
        worker = basioConverter(work_queue, results_queue, fastafile, i+1)
        worker.start()
    
    # Collect the results off the queue and merge them
    print "Collecting results..."
    distree = oak.plant()
    for i in range(args.cpu): 
        graft = results_queue.get()
        if graft:
            sys.stdout.write("graft %d ...\n" % (i + 1))
            oak.graft(distree, graft)
            sys.stdout.write("graft %d ok\n" % (i + 1))
    
    # Save the results:
    print "Saving results..."
    if args.dumpmode == "pickle" or args.dumpmode == "both":
        pass #TODO: implement that sh** too
    if args.dumpmode == "mv" or args.dumpmode == "both":
#        oak.save_to_file(distree, os.path.splitext(args.fragsfile)[0], args.isreference)
        oak.save_to_file(distree, os.path.splitext(args.fragsfile)[0])
#        file1 = open("%s.mvseg" % os.path.splitext(args.fragsfile)[0], "w")
#        seglist = distree[2].keys()
#        for i in range(len(seglist)):
#            file1.write("%d %s\n" % (i, seglist[i]))
#        file1.close()
#        file2 = open("%s.mvdist" % os.path.splitext(args.fragsfile)[0], "w")
#        for i in range(len(distree[0])):
#            for j in range(len(seglist)):
#                if distree[2][seglist[j]][i] != 0:
#                    file2.write("%s %d %d %d\n" % (distree[0][i], distree[1][i], j, distree[2][seglist[j]][i]))
#        file2.close()
    
    print "Done :)"
    
    
    
    
