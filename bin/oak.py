# -*- coding: utf-8 -*-

import re, random, sys, os, numpy



## Save the forest, plant a tree
#def plant():
#    print 'Seed planted! Let\'s water it now.'
#    return [ [], {"A":[],"T":[],"G":[],"C":[]} ]

def plant():
    return [[], [], {}]

# Adds a new length the sequence length list
# TO BE RUN FIRST
def start_new_sequence(tree, length, seqid):
    for key in tree[2].keys():
        while len(tree[1]) > len(tree[2][key]):
            tree[2][key].append(0.0)
    tree[0].append(seqid)
    tree[1].append(length)

## Add 1 to a specific segment counter (and creates required nodes if needed)
#def add_segment(tree, segment, step=1):
#    segment = segment.upper()
#    if not re.match('^[ATGC]*$', segment):
#        print "segment %s not valid! ignored." % segment
#        return False
#    current_node = tree
#    for i in range(len(segment)):
#        if current_node[1][segment[i]] == []:
#            current_node[1][segment[i]] = [ [], {"A":[],"T":[],"G":[],"C":[]} ]
#        current_node = current_node[1][segment[i]]
#    while len(current_node[0]) < len(tree[1]):
#        current_node[0].append(0.0)
#    current_node[0][len(current_node[0]) - 1]= current_node[0][len(current_node[0]) - 1] + int(step)
#    return True

def add_segment(tree, segment, step=1):
    segment = segment.upper()
    if not re.match('^[ATGC]*$', segment):
        print "segment %s not valid! ignored." % segment
        return False
    if segment not in tree[2]:
        tree[2][segment] = []
    while len(tree[2][segment]) < len(tree[1]):
        tree[2][segment].append(0)
    tree[2][segment][len(tree[2][segment]) - 1] = tree[2][segment][len(tree[2][segment]) - 1] + int(step)

def complete_tree(tree):
    length = len(tree[1])
    for key in tree[2].keys():
        while len(tree[2][key]) < length:
            tree[2][key].append(0)

def old_graft(host, graft):
    for i in range(len(graft[0])):
        host[0].append(graft[0][i])
        host[1].append(graft[1][i])
        for segment in graft[2].keys():
            add_segment(host, segment, graft[2][segment][i])
    complete_tree(host)

def graft(host, graft):
    lenhost = len(host[0])
    lengraft = len(graft[0])
    hostsegs = set(host[2].keys())
    graftsegs =set(graft[2].keys())
    host[0] = host[0] + graft[0]
    host[1] = host[1] + graft[1]
    for commonseg in hostsegs & graftsegs:
        host[2][commonseg] = host[2][commonseg] + graft[2][commonseg]
    for hostonly in hostsegs - graftsegs:
        host[2][hostonly] = host[2][hostonly] + [0]*lengraft
    for graftonly in graftsegs - hostsegs:
        host[2][graftonly] = [0]*lenhost + graft[2][graftonly]
        

def save_to_file_old(tree, filenameprefix):
    seglist = tree[2].keys()
    mvfile = open("%s.mv" % filenameprefix, "w")
    mvfile.write("seqid seqlen %s\n" % ' '.join(seglist) )
    for i in range(len(tree[0])):
        sys.stdout.write("%d\n" % i)
        mvfile.write( "%s %s" % (tree[0][i], tree[1][i]) )
        for seg in seglist:
            mvfile.write( " %s" % tree[2][seg][i] )
        mvfile.write("\n")
    mvfile.close()
    
def save_to_file_old2(tree, filenameprefix):
    seglist = tree[2].keys()
    mvfile = open("%s.mv" % filenameprefix, "w")
    mvfile.write("seqid %s\n" % ' '.join(tree[0]) )
    mvfile.write("seqlen %s\n" % ' '.join(str(e) for e in tree[1]) )
    for seg in seglist:
#        sys.stdout.write("%s\n" % seg)
        mvfile.write( "%s %s\n" % ( seg, ' '.join(str(e) for e in tree[2][seg]) ) )
    mvfile.close()

def save_to_file_old3(tree, filenameprefix, isref):
    # Saving the segments in the .mvseg file
#    file1 = open("%s.mvseg" % os.path.splitext(filenameprefix)[0], "w")
    file1 = open("%s.mvseg" % filenameprefix, "w")
    seglist = tree[2].keys()
    for i in range(len(seglist)):
        file1.write("%d %s\n" % (i, seglist[i]))
    file1.close()
    # Saving the counts in the .mvdist file
    
    if isref:
#        file2 = open("%s.mvref" % os.path.splitext(filenameprefix)[0], "w")
        file2 = open("%s.mvref" % filenameprefix, "w")
        for j in range(len(seglist)):
            file2.write("%s %d %d %d\n" % (filenameprefix, numpy.sum(tree[1]), j, numpy.sum(tree[2][seglist[j]]) ))
    else:
#        file2 = open("%s.mvdist" % os.path.splitext(filenameprefix)[0], "w")
        file2 = open("%s.mvdist" % filenameprefix, "w")
        for i in range(len(tree[0])):
            for j in range(len(seglist)):
                if tree[2][seglist[j]][i] != 0:
                    file2.write("%s %d %d %d\n" % (tree[0][i], tree[1][i], j, tree[2][seglist[j]][i]))
    file2.close()
    
def save_to_file(tree, filenameprefix):
    # Saving the counts in the .mvdist file
    file2 = open("%s.mv" % filenameprefix, "w")
    for i in range(len(tree[0])):
        for seg in tree[2].keys():
            if tree[2][seg][i] != 0:
                file2.write("%s %d %s %d\n" % (tree[0][i], tree[1][i], seg, tree[2][seg][i]))
    file2.close()

def read_from_file(filename):
    mvfile = open(filename, "r")
    headers = mvfile.readline().split()
    if headers[0] == "seqid" and headers[1] == "seqlen":
        tree = plant()
        for line in mvfile.readlines():
            line = line.split()
            tree[0].append(line[0])
            tree[1].append(line[1])
            for i in range(2, len(headers)):
                tree[2][headers[i]].append(line[i])
        mvfile.close()
        return tree
    else:
        print "Wrong file?"
        mvfile.close()
        return None





