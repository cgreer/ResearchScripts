#!/usr/bin/python

"""
visualize_splice_graphs.py:

    visualize all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph, Classifier
from optparse import OptionParser, make_option

# default directory
#defaultDir = '/r100/burge/shared/splice_graphs/mm6/FlatFiles'
defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles'
defaultPrefix = 'SpliceGraphs_'

# digest command line
opts = OptionParser(option_list=                                                                             \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                             \
                                 help='input directory with splice graph sub-dirs [default: %s]' % defaultDir),  \
                     make_option('-p','--inprefix', dest='inprefix', default=defaultPrefix,                  \
                                 help='iterate over dirs with this prefix [default: %s]' % defaultPrefix),   \
                     make_option('-c','--chr', dest='chr', default=None,                                     \
                                 help='only do one chromosome [default: do all chromosomes]'),   \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(options.inDir, prefix=options.inprefix)
    print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)

    # iterate over splice graphs
    current_dir = iterChrom.current_dir()
    if options.chr != None:
        while current_dir != "%s/%s%s" % (options.inDir, options.inprefix, options.chr) and current_dir != None:
            current_dir = iterChrom.next_dir()
    sg = iterChrom.next_sg()
    nb = 0
    sys.stderr.write("visualizing splice graphs in %s\n" % current_dir)
    while ( sg ):
        nb += 1
        sys.stderr.write("drawing splice graph %d of %d\r" % (nb, iterChrom.number()))
        
        try:
            sg.visualize("%s/%s.ps" % (iterChrom.current_dir(), sg.name))
        except:
            sys.stderr.write("WARNING: Could not visualize gene %s: %s\n" % (sg.name, sys.exc_info()[1]) )

        # get the next splice graph
        sg = iterChrom.next_sg()

        # a new splice graph directory was started?
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            if options.chr != None:
                # only one chromosome to do --> finished
                sys.exit(0)
            current_dir = iterChrom.current_dir()
            nb = 0
            sys.stderr.write("classifying splice graphs in %s\n" % current_dir)
