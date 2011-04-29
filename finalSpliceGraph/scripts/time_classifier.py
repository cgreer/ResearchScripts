#!/usr/local/bin/python

"""
time_classifier.py:

    time a certain task using Classifier.py

"""

import os, sys, time
import SpliceGraph, Classifier

repeats = 1

# input splice graph files (decreasing file size)
# obtained by: `ls -ltr /r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/ | grep '\.sg$' | sort -n -r -k 5 | less`
## sgFiles = [ \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:9222362:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:153121268:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:93009603:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:45656702:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:150456552:-.sg', \
##     '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:170568784:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:224577096:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:115012666:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:157982178:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:224634717:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:148040270:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:241353868:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:151005728:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:157126135:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:162332108:+.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:148549113:-.sg', \
##     #'/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:24722925:+.sg', \
##     ]

inDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr18'
sgFiles = [ "%s/%s" % (inDir,f) for f in os.listdir(inDir) if f.endswith('.sg') ]

# do preliminary work
sys.stdout.write("averaging %i repeats of each task\n" % repeats)
sgs = {}
for sgFile in sgFiles:
    try:
        sg = SpliceGraph.SpliceGraph(filename=sgFile)
    except:
        pass
    else:
        sgs[sg.name] = sg
classifier = Classifier.Classifier()
fh = open('tmp_1', 'w')
started = time.time()

# do task 1
for rep in xrange(repeats):
    for name in sgs.iterkeys():
#        CEs = classifier.classify_CE_obsolete(sgs[name])
#        for string in classifier.prettyString(CEs, 'CE', name):
#            fh.write(string)
#        CIs = classifier.classify_CI_obsolete(sgs[name])
#        for string in classifier.prettyString(CIs, 'CI', name):
#            fh.write(string)
        MXEs = classifier.classify_MXEobsolete(sgs[name])
        for string in classifier.prettyString(MXEs, 'MXEobsolete', name):
            fh.write(string)
        
# measure time interval 1
interval1 = float(time.time() - started)
sys.stdout.write("task 1 average: %.3f seconds of real time\n" % (interval1/repeats))

# do intermediate work
fh.close()
fh = open('tmp_2', 'w')
events = {}
started = time.time()

# do task 2
for rep in xrange(repeats):
    for name in sgs.iterkeys():
#        CEs = classifier.classify_CE(sgs[name])
#        for string in classifier.prettyString(CEs, 'CE', name):
#            fh.write(string)
#        CIs = classifier.classify_CI(sgs[name])
#        for string in classifier.prettyString(CIs, 'CI', name):
#            fh.write(string)
        MXEs = classifier.classify_MXE(sgs[name])
        for string in classifier.prettyString(MXEs, 'MXE', name):
            fh.write(string)

# measure time interval 2
interval2 = float(time.time() - started)
sys.stdout.write("task 2 average: %.3f seconds of real time\n" % (interval2/repeats))
sys.stdout.write("                %.2f%% of task 1\n" % (interval2*100/interval1))

# compare the results
fh.close()
sys.stdout.write("compare results in files tmp_1 (task 1) and tmp_2 (task 2)\n")
