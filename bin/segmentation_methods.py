# -*- coding: utf-8 -*-

import os, tempfile, oak, subprocess, glob, threading, sys
from Bio import SeqIO
#from cStringIO import StringIO
# or, in Py3/Py2.6+:
from io import StringIO
from itertools import ifilter


class segmentationRunner(threading.Thread):
    """
    Hello, I'm a thread class, and I'm in charge of one particular fasta
    sequence per instance/thread. I'll prepare the data, launch the segmentation
    subprocess, safely add the results to the shared tree, and clean up
    afterwards 'cause I'm a good guy. Have a nice day!
    """

    def __init__(self, args, record, distree, dontCreateAMonster_Lock):
        threading.Thread.__init__(self)
        self.args = args
        self.record = record
        self.distree = distree
        self.dontCreateAMonster_Lock = dontCreateAMonster_Lock
        self.fastafile = ""
        self.shrub = oak.plant()
        self.bin = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "bin"))

    def log(self, msg):
        sys.stdout.write("[%s] %s\n" % (self.record.id, msg))

    def run(self):
        """
        Triggers the appropriate function as specified in the args,
        manages the boring stuff of the tree and creates a temp file for each
        sequence.
        """
        # BEFORE adding all the segment_methods, we add the sequence length
        # (since it's used to normalize stuff)
        oak.start_new_sequence(self.shrub, len(self.record.seq))
        # Create a temporary file containing the current FASTA sequence
        # (since most of the guys can't handle multifasta)
        tmpfd, self.fastafile = tempfile.mkstemp(suffix=".fa")
    #    print tmpfilepath
        SeqIO.write([self.record], self.fastafile, "fasta")
        # Do stuff with the temporary files
        exec "self.%s()" % self.args.segment_method
        # Delete the tempfile
        os.close(tmpfd) # To prevent file descriptors exhaustion
        os.remove(self.fastafile)
        self.log("Locking access to shared data")
        self.dontCreateAMonster_Lock.acquire()
        try:
            oak.graft(self.distree, self.shrub)
        finally:
            self.log("releasing lock")
            self.dontCreateAMonster_Lock.release()


    def kmers(self):
        """
        kmers counter (Thomas' version)
        """
        # Get k or the range of k values from options
        k_values = []
        k_from_args = self.args.k.strip('"').split("-")
        if len(k_from_args) == 2:
            k_values = range( int(k_from_args[0]), int(k_from_args[1]) + 1 )
        elif len(k_from_args) == 1:
            k_values = [int(k) for k in k_from_args]
        else:
            raise Exception("Inappropriate k value or range")
    #    print "k values: ", k_values

        for k in k_values:
    #        print "k = ", k
    #        print "file = ", self.fastafile
            cmd = "kmer %d %s" % (k, self.fastafile)
            exec_command(cmd)
            dump = open("%s.frags" % self.fastafile.split(".")[0])
            for kmer in dump.readlines():
                kmer = kmer.split()
                if len(kmer) > 1:
                    oak.add_segment(self.shrub, kmer[0], int(kmer[1]))
    #                if str(k) in kcounter and not kmer[0] in kcounter[str(k)]:
    ##                    print "adding ", kmer[1]
    #                    kcounter[str(k)].append(kmer[0])
    #                else:
    ##                    print "creating ", str(k), " - ", kmer[1]
    #                    kcounter[str(k)] = [kmer[0]]
            dump.close()
        try:
            os.remove("%s.frags" % self.fastafile.split(".")[0])
        except:
            self.log("missing files to delete")


    def basio(self):
        """
        Basio - Bayesian segmentation
        """
        #TODO: check if args.bip is an actual float!!!!!
        fasta_splitpath = os.path.splitext(self.fastafile)
        # First: Run basio
        if len(self.record.seq) > 15000: # for long sequences (probably not viral)
            self.log("splitting big file and converting...")
            cmd = "%s/basio_splitb --length=10000 --overlap=1000 %s" % (self.bin, self.fastafile)
            exec_command(cmd)
            self.log("segmenting...")
            for segbin in glob.glob( "%s*_0.sgm" % fasta_splitpath[0] ):
                cmd = "%s/basio_basio --bip=%s %s.sgm" % ( self.bin, self.args.bip, segbin )
                exec_command(cmd)
            self.log("merging...")
            cmd = "%s/basio_merge %s" % (self.bin, fasta_splitpath[0])
            exec_command(cmd)
            for segbin in glob.glob( "%s_0*.sgm" % fasta_splitpath[0] ):
                os.remove(segbin)
            self.log("start optimisation")
        else:
            self.log("converting file")
            cmd = "%s/basio_splitb %s" % (self.bin, self.fastafile)
            exec_command(cmd)
            self.log("start segmentation")
        cmd = "%s/basio_basio --bip=%s %s.sgm" % ( self.bin, self.args.bip, fasta_splitpath[0] )
        exec_command(cmd)
        cmd = "%s/basio_control -v %s.sgm" % (self.bin, fasta_splitpath[0])
        exec_command(cmd)
        # Second: retrieve the segments' sequences and feed the tree
        self.log("filling the tree...")
        contig = StringIO(unicode(self.record.seq) )
        previousline = None
        with open("%s.sgm" % fasta_splitpath[0], 'r') as segfile:
            for line in ifilter(lambda line: line.startswith('BR'), segfile):
                line = line.split()
                if previousline:
                    segment = contig.read( int(line[1]) - int(previousline[1]) )
                    oak.add_segment(self.shrub, segment)
                previousline = line
        segfile.close()
        # Erase all the evidence
        os.remove( "%s.sgm" % fasta_splitpath[0] )
        self.log("sequence done.")





def exec_command( cmd ):
    """
    Launches the execution of a command in a shell.
    """
    proc = subprocess.Popen( args=cmd,
                 shell=True )
    returncode = proc.wait()
    if returncode != 0:
        raise Exception( "Bad exit status for command: %s" % cmd )
#
#def run(args, record, distree):
#    """
#    Triggers the appropriate function as specified in the args,
#    manages the boring stuff of the tree and creates a temp file for each
#    sequence.
#    """
#    # BEFORE adding all the segment_methods, we add the sequence length
#    # (since it's used to normalize stuff)
#    oak.start_new_sequence(distree, len(record.seq))
#    # Create a temporary file containing the current FASTA sequence
#    # (since most of the guys can't handle multifasta)
#    tmpfd, tmpfilepath = tempfile.mkstemp(suffix=".fa")
##    print tmpfilepath
#    SeqIO.write([record], tmpfilepath, "fasta")
#    # Do stuff with the temporary files
#    exec "%s(tmpfilepath, record, distree, args)" %args.segment_method
#    # Delete the tempfile
#    os.close(tmpfd) # To prevent file descriptors exhaustion
#    os.remove(tmpfilepath)
#
#def kmers(fastafile, record, tree, args):
#    """
#    kmers counter (Thomas' version)
#    """
#    # Get k or the range of k values from options
#    k_values = []
#    k_from_args = args.k.strip('"').split("-")
#    if len(k_from_args) == 2:
#        k_values = range( int(k_from_args[0]), int(k_from_args[1]) + 1 )
#    elif len(k_from_args) == 1:
#        k_values = [int(k) for k in k_from_args]
#    else:
#        raise Exception("Inappropriate k value or range")
##    print "k values: ", k_values
#
#    for k in k_values:
##        print "k = ", k
##        print "file = ", fastafile
#        cmd = "kmer %d %s" % (k, fastafile)
#        exec_command(cmd)
#        dump = open("%s.frags" % fastafile.split(".")[0])
#        for kmer in dump.readlines():
#            kmer = kmer.split()
#            if len(kmer) > 1:
#                oak.add_segment(tree, kmer[0], int(kmer[1]))
##                if str(k) in kcounter and not kmer[0] in kcounter[str(k)]:
###                    print "adding ", kmer[1]
##                    kcounter[str(k)].append(kmer[0])
##                else:
###                    print "creating ", str(k), " - ", kmer[1]
##                    kcounter[str(k)] = [kmer[0]]
#        oak.complete_tree(tree)
#        dump.close()
#    try:
#        os.remove("%s.frags" % fastafile.split(".")[0])
#    except:
#        print "missing files to delete"
#
#
#def basio(fastafile, record, tree, args):
#    """
#    Basio - Bayesian segmentation
#    """
#    #TODO: check if args.bip is an actual float!!!!!
#    fasta_splitpath = os.path.splitext(fastafile)
#    # First: Run basio
#    if len(record.seq) > 15000: # for long sequences (probably not viral)
#        print "splitting big file and converting..."
#        cmd = "basio_splitb --length=10000 --overlap=1000 %s" % fastafile
#        exec_command(cmd)
#        print "segmenting..."
#        for segbin in glob.glob( "%s*_0.sgm" % fasta_splitpath[0] ):
#            cmd = "basio_segment --bip=%s %s.sgm" % ( args.bip, segbin )
#            exec_command(cmd)
#        print "merging..."
#        cmd = "basio_merge %s" % fasta_splitpath[0]
#        exec_command(cmd)
#        for segbin in glob.glob( "%s_0*.sgm" % fasta_splitpath[0] ):
#            os.remove(segbin)
#        print "start optimisation"
#    else:
#        print "converting file"
#        cmd = "basio_splitb %s" % fastafile
#        exec_command(cmd)
#        print "start segmentation"
#    cmd = "basio_segment --bip=%s %s.sgm" % ( args.bip, fasta_splitpath[0] )
#    exec_command(cmd)
#    cmd = "basio_control -v %s.sgm" % fasta_splitpath[0]
#    exec_command(cmd)
#    # Second: retrieve the segments' sequences and feed the tree
#    print "filling the tree..."
#    contig = StringIO(unicode(record.seq) )
#    previousline = None
#    with open("%s.sgm" % fasta_splitpath[0], 'r') as segfile:
#        for line in ifilter(lambda line: line.startswith('BR'), segfile):
#            line = line.split()
#            if previousline:
#                segment = contig.read( int(line[1]) - int(previousline[1]) )
#                oak.add_segment(tree, segment)
#            previousline = line
#    segfile.close()
#    # Erase all the evidence
#    os.remove( "%s.sgm" % fasta_splitpath[0] )
#    print "sequence done."
#
