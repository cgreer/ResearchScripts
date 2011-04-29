#!/usr/local/bin/python

"""
find_elements.py:

    find elements with certain criteria in all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph, Classifier
from optparse import OptionParser, make_option

import psyco
psyco.full()

# default directory
#defaultDir = '/r100/burge/shared/splice_graphs/mm6/FlatFiles'
defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles'
defaultMinCoverage = '2'
defaultCanonical = '0'
defaultOutfile = 'found_elements.txt'

# digest command line
opts = OptionParser(option_list=                                                                                  \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                                  \
                                 help='input directory with SpliceGraph_* dirs [default: %s]' % defaultDir),      \
                     make_option('-m','--mincoverage', dest='mincoverage', default=defaultMinCoverage,            \
                                 help='min. coverage of elements to find [default: %s]' % defaultMinCoverage),    \
                     make_option('-c','--canonical', dest='canonical', default=defaultCanonical,                  \
                                 help='find canonical elements (1:yes, 0:no)? [default: %s]' % defaultCanonical), \
                     make_option('-o','--outfile', dest='outfile', default=defaultOutfile,                        \
                                 help='output filename [default: %s]' % defaultOutfile), \
                     ])
(options, args) = opts.parse_args()
options.canonical = bool(int(options.canonical))
options.mincoverage = int(options.mincoverage)


if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(options.inDir, prefix='SpliceGraphsFiltered_')
    print "#there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)

    # open output file
    outFile = open(options.outfile, 'w')
    outFile.write("#output of find_elements.py\n")
    outFile.write("#mincoverage: %i\n" % options.mincoverage)
    outFile.write("#canoncial: %s\n" % (options.canonical and "yes" or "no"))
    outFile.write("#outfile: %s\n" % options.outfile)
    outFile.write("#gnId\tel\tcoverage\tcanoncial\n")

    # iterate over splice graphs
    sg = iterChrom.next_sg()
    current_dir = iterChrom.current_dir()
    nb = 0
    sys.stderr.write("classifying splice graphs in %s\n" % current_dir)
    while ( sg ):
        nb += 1
        sys.stderr.write("analyzing splice graph %d of %d\r" % (nb, iterChrom.number()))


        # find elements
        for el in sg.allElements():

            if sg.getType(el) in ['5ss','3ss']:
                elSS = SpliceGraph.SS(el, sg.getType(el))
                if sg.coverage(el) >= options.mincoverage and \
                   ( (options.canonical is True  and     elSS.isCanonical()) or \
                     (options.canonical is False and not elSS.isCanonical()) ):
                    outFile.write("%s\t%s\t%i\t%s\n" % (sg.name, el, sg.coverage(el), (elSS.isCanonical() and '1' or '0')))


        # get the next splice graph
        sg = iterChrom.next_sg()

        # begin new output directory if a new splice graph directory was started
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            # iterChrom has gone to the next chromosome --> begin new directory
            current_dir = iterChrom.current_dir()
            nb = 0
            sys.stderr.write("classifying splice graphs in %s\n" % current_dir)

    outFile.close()
