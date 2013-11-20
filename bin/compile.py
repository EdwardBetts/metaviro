# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:15:17 2013

@author: Louise-Amélie Schmitt

For each method, we need 2 fonctions: one that checks if the binaries exist
in the bin directory, and one for the compilation.

You can freely add as many as you want, as long as you are coherent in the
names of your methods, i.e. make sure the new directory in the src dir is 
named with the same thing as you use here. (in the list below and the names 
of the new fonctions)
"""

import sys, os.path
from lib.segmentation_methods import exec_command

### List of segmentation methods! ###
methods = ["basio", "kmers"]

class NastyCompilationError(Exception):
    def __init__(self, e):
        self.msg = "Compilation failed. (\"%s\")" % (str(e))
    def __str__(self):
        return self.msg

# Just because. Laziness. Takes the codename of the method as an argument.
def get_abs_paths(name):
    src_abspath = os.path.abspath(os.path.join(os.path.dirname(__file__), "src/%s" % name))
    bin_abspath = os.path.abspath(os.path.join(os.path.dirname(__file__), "bin"))
    return src_abspath, bin_abspath

# One function to rule them all, One function to find them,
# One function to bring them all and in the library bind them
# In the Land of MetaviroII where the libraries lie
def check_binaries():
    for method in methods:
        exec """if not %s_binaries_exist():
    print \"Missing files for %s\"
    compile_%s()""" % (method, method, method)

############
## K-mers ##
############

def kmers_binaries_exist():
    srcdir, bindir = get_abs_paths("kmers")
    if os.path.isfile("%s/kmer" % bindir):
        return True
    return False

def compile_kmers():
    srcdir, bindir = get_abs_paths("kmers")
    cmd = "module load gcc && g++ -W -Wall -std=c++11 -O2 -o %s/kmer %s/kmer.cpp" % (bindir, srcdir)
    try:
        print "compiling kmer binaries..."
        exec_command(cmd)
        print "done compiling"
    except Exception as e:
        raise NastyCompilationError(e)

###########
## Basio ##
###########

def basio_binaries_exist():
    srcdir, bindir = get_abs_paths("kmers")
    binaries=["basio", "control", "filter", "merge", "splitb"]
    for binary in binaries:
        if not os.path.isfile("%s/basio_%s" % (bindir, binary)):
            return False
    return True

def compile_basio():
    binaries=["basio", "control", "filter", "merge", "splitb"]
    srcdir, bindir = get_abs_paths("basio")
    print "compiling kmer binaries..."
    for binary in binaries:
        cmd = "gcc -lm -Wall -o %s/basio_%s %s/%s.c" % (bindir, binary, srcdir, binary)
        try:
            print "processing %s ..." % binary
            exec_command(cmd)
        except Exception as e:
            raise NastyCompilationError(e)
    print "done compiling"


#####################
## Run as a script ##
#####################

if __name__ == '__main__':
    check_binaries()