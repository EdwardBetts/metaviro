# create matrices from mv file

from __future__ import with_statement
import time, re
from scipy.sparse import coo_matrix, csr_matrix, dok_matrix, lil_matrix
from sklearn.metrics.pairwise import euclidean_distances


class Timer(object):
    def __enter__(self):
        self.__start = time.time()

    def __exit__(self, type, value, traceback):
        # Error handling here
        self.__finish = time.time()

    def duration_in_seconds(self):
        return self.__finish - self.__start

timer = Timer()

#ncontigs = 0
#with open("virus_cont_numeric.fasta", "r") as fastafile:
#    for line in fastafile:
#        if re.match('>', line):
#            ncontigs+=1
#
#print "-- COO --"
#for k in range(3, 8):
#    data = []
#    rows = []
#    cols = []
#    
#    with timer:
#        with open("virus_cont_numeric_k%d.mv"%k, "r") as mvfh:
#            for line in mvfh:
#                line = line.split()
#                rows.append(int(line[0].split("_")[0]))
#                cols.append(int(line[2], 4)) # base 4 to 10
#                data.append(int(line[3]))
#    readtime = timer.duration_in_seconds()
#            
#    with timer:
#        coo = coo_matrix((data, (rows, cols)))
#    
#    
#
#    print k, readtime, timer.duration_in_seconds()

reigns = ["arch", "bact", "euk", "virus"]

print "-- CSR --"
for k in range(7, 8):
    data = []
    rows = []
    cols = []
    targ = []
    contig_to_row = {} # (reign, contig_id) -> row_id
    rowid = 0
    with timer:
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
                        
                    
    readtime = timer.duration_in_seconds()
            
    with timer:
        csr = csr_matrix((data, (rows, cols)), shape=(len(contig_to_row), 4**k))

    print k, readtime, timer.duration_in_seconds()
#
#print "-- DOK --"
#for k in range(3, 8):
    #data = []
    #rows = []
    #cols = []
    #
    #with timer:
        #with open("virus_cont_numeric_k%d.mv"%k, "r") as mvfh:
            #for line in mvfh:
                #line = line.split()
                #rows.append(int(line[0].split("_")[0]))
                #cols.append(int(line[2], 4)) # base 4 to 10
                #data.append(int(line[3]))
    #readtime = timer.duration_in_seconds()
            #
    #with timer:
        #dok = dok_matrix((ncontigs, 4**k))
        #for row, col, value in zip(rows, cols, data):
            #dok[row, col] = value
#
    #print k, readtime, timer.duration_in_seconds()
#
#print "-- LIL --"
#for k in range(3, 8):
    #data = []
    #rows = []
    #cols = []
    #
    #with timer:
        #with open("virus_cont_numeric_k%d.mv"%k, "r") as mvfh:
            #for line in mvfh:
                #line = line.split()
                #rows.append(int(line[0].split("_")[0]))
                #cols.append(int(line[2], 4)) # base 4 to 10
                #data.append(int(line[3]))
    #readtime = timer.duration_in_seconds()
            #
    #with timer:
        #lil = lil_matrix((ncontigs, 4**k))
        #for row, col, value in zip(rows, cols, data):
            #lil[row, col] = value
#
    #print k, readtime, timer.duration_in_seconds()
