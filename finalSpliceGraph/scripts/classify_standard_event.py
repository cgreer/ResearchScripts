#!/usr/bin/env python

"""
classify_standard_event.py:

    classify standard splicing event for all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph, Classifier, LocusFilter
from optparse import OptionParser, make_option

import psyco
#psyco.log()
#psyco.profile()
psyco.full()

# default directory
defaultAssembly = 'hg17'
defaultDir = '/r100/burge/shared/splice_graphs/%s/FlatFiles' % defaultAssembly
defaultFilter = None

# digest command line
opts = OptionParser(option_list=                                                                             \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                             \
                                 help='input directory with SpliceGraph* dirs [default: %s]' % defaultDir),  \
                     make_option('-a','--assembly', dest='assembly', default=defaultAssembly,                \
                                 help='assembly (used for LocusFilter) [default: %s]' % defaultAssembly),    \
                     make_option('-e','--event', dest='event', default='MXE',                                \
                                 help='standard event type to classify [default: MXE]'),                     \
                     make_option('-p','--inprefix', dest='inprefix', default='SpliceGraphs_',                \
                                 help='iterate over dirs with this prefix [default: SpliceGraphs_]'),        \
                     make_option('-f','--filter', dest='filter', default=defaultFilter,                      \
                                 help='skip splice graphs by LocusFilter regions from this file [default: %s]' % defaultFilter),        \
                     make_option('-o','--outprefix', dest='outprefix', default='SpliceGraphs_',              \
                                 help='store events in dirs with this prefix [default: StandardEvents_]'),   \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    lf = None
    if options.filter != None:
        lf = LocusFilter.LocusFilter(filename=options.filter, verbose=0)
        sys.stderr.write("going to use LocusFilter:\n\t%s" % lf)
    c = Classifier.Classifier()
    classifyMethod = getattr(c, "classify_%s" % options.event)
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(options.inDir, prefix=options.inprefix)
    print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)

    # open output file
    dirname  = iterChrom.topdir + '/' + options.outprefix + (os.path.basename(iterChrom.current_dir()).split('_',1)[-1])
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        os.chmod(dirname, 0770)
    filename = options.event + '_events'
    if os.path.exists('%s/%s' % (dirname, filename)):
        os.remove('%s/%s' % (dirname, filename))
    outFile  = open('%s/%s' % (dirname, filename), 'w')
    for string in c.prettyString(None, options.event):
        outFile.write(string) #write header to output file

    # iterate over splice graphs
    sg = iterChrom.next_sg()
    current_dir = iterChrom.current_dir()
    nb = 0
    sys.stderr.write("classifying splice graphs in %s\n" % current_dir)
    while ( sg ):
        nb += 1
        sys.stderr.write("analyzing splice graph %d of %d (%s)\r" % (nb, iterChrom.number(), sg.name))

        if lf == None or ( isinstance(lf, LocusFilter.LocusFilter) and not lf.overlaps_sg( sg ) ):
            events = classifyMethod(sg)
            for string in c.prettyString(events, options.event, sg.name):
                outFile.write(string)
        else:
            sys.stderr.write("skipped %s: overlapps %s locus                       \n" % (sg.name, lf.lastOverlap()))

        # get the next splice graph
        sg = iterChrom.next_sg()

        # begin new output file if a new splice graph directory was started
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            # iterChrom has gone to the next chromosome --> begin new outFile
            outFile.close()
            os.chmod('%s/%s' % (dirname, filename), 0660)
            dirname  = iterChrom.topdir + '/' + options.outprefix + (os.path.basename(iterChrom.current_dir()).split('_',1)[-1])
            if not os.path.exists(dirname):
                os.mkdir(dirname)
                os.chmod(dirname, 0770)
            if os.path.exists('%s/%s' % (dirname, filename)):
                os.remove('%s/%s' % (dirname, filename))
            outFile  = open('%s/%s' % (dirname, filename), 'w')
            for string in c.prettyString(None, options.event):
                outFile.write(string) #write header to output file
            current_dir = iterChrom.current_dir()
            nb = 0
            sys.stderr.write("classifying splice graphs in %s\n" % current_dir)

    outFile.close()
    os.chmod('%s/%s' % (dirname, filename), 0660)
