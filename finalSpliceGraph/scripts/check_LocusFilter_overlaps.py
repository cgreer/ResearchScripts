#!/usr/bin/env python

"""
check_LocusFilter_overlaps.py:

    For all *.sg files for all chromosomes in a given directory, check overlaps with LocusFilter regions

"""

import os, sys, SpliceGraph, LocusFilter
from optparse import OptionParser, make_option

import psyco
psyco.full()


# default directory
defaultAssembly = 'hg17'
defaultDir      = '/r100/burge/shared/splice_graphs/%s/FlatFiles' % defaultAssembly
#defaultPrefix   = 'SpliceGraphsFiltered_chr14'
defaultPrefix   = 'SpliceGraphsFiltered_'

# digest command line
opts = OptionParser(option_list=                                                                             \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                             \
                                 help='input directory with SpliceGraph* dirs [default: %s]' % defaultDir),  \
                     make_option('-p','--prefix', dest='prefix', default=defaultPrefix,                      \
                                 help='iterate over dirs with this prefix [default: %s]' % defaultPrefix),   \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(options.inDir, prefix=options.prefix)
    sys.stderr.write("\tthere are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir))

    lf = LocusFilter.LocusFilter()

    # iterate over splice graphs
    sg = iterChrom.next_sg()
    current_dir = iterChrom.current_dir()
    sys.stderr.write("\n\tsearching *.sg in %s\n" % current_dir)
    while ( sg ):
        if lf.overlaps_sg(sg):
            print '%s : %s/%s' % (lf.lastOverlap(), iterChrom.current_dir(), iterChrom.diriter.files[iterChrom.diriter.current_index()])
        sg = iterChrom.next_sg()
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            current_dir = iterChrom.current_dir()
            sys.stderr.write("\n\tsearching *.sg in %s\n" % current_dir)
