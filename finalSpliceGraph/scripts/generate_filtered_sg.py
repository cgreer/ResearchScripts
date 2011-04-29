#!/usr/local/bin/python

"""
generate_filtered_sg.py:

    generate filtered versions for all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph, Classifier
from optparse import OptionParser, make_option

# default directory
#defaultDir = '/r100/burge/shared/splice_graphs/mm6/FlatFiles'
defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles'
defaultMethod = '1'
defaultInPrefix = 'SpliceGraphs_'
defaultOutPrefix = 'SpliceGraphsCanonical_'

# digest command line
opts = OptionParser(option_list=                                                                             \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                             \
                                 help='input directory with SpliceGraph_* dirs [default: %s]' % defaultDir), \
                     make_option('-m','--method', dest='method', default=defaultMethod,                      \
                                 help='filtering method (1:remove_non_canonical, 2:remove_by_coverage_and_ss, 3:tweak_non_canonical) [default: 1]'),  \
                     make_option('-i','--inprefix', dest='inprefix', default=defaultInPrefix,                \
                                 help='input dir prefix [default: %s]' % defaultInPrefix),                   \
                     make_option('-o','--outprefix', dest='outprefix', default=defaultOutPrefix,             \
                                 help='output dir prefix [default: %s]' % defaultOutPrefix),                 \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(options.inDir, prefix=options.inprefix)
    print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)

    # open output file
    dirname  = iterChrom.topdir + '/' + options.outprefix + (os.path.basename(iterChrom.current_dir()).split('_',1)[-1])
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        os.chmod(dirname, 0770)

    # iterate over splice graphs
    sg = iterChrom.next_sg()
    current_dir = iterChrom.current_dir()
    nb = 0
    sys.stderr.write("filtering splice graphs in %s\n" % current_dir)
    while ( sg ):
        nb += 1
        sys.stderr.write("analyzing splice graph %d of %d\r" % (nb, iterChrom.number()))

        ### FILTERING IS HERE ###
        if options.method == '1':
            newsg = sg.remove_non_canonical()
        elif options.method == '2':
            newsg = sg.remove_by_coverage_and_ss(coverage=2)
        elif options.method == '3':
            newsg = sg.tweak_non_canonical(nbBases=3)
        else:
            sys.stderr.write("ERROR: unknown filtering mode %s. Aborting" % options.method)
            sys.exit(0)

        try:
            newsg.store("%s/%s.sg" % (dirname, sg.name))
        except:
            sys.stderr.write("WARNING: Could not store filtered gene %s: %s\n" % (sg.name, sys.exc_info()[1]) )

        # get the next splice graph
        sg = iterChrom.next_sg()

        # begin new output directory if a new splice graph directory was started
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            # iterChrom has gone to the next chromosome --> begin new directory
            dirname  = iterChrom.topdir + '/' + options.outprefix + (os.path.basename(iterChrom.current_dir()).split('_',1)[-1])
            if not os.path.exists(dirname):
                os.mkdir(dirname)
                os.chmod(dirname, 0770)
            current_dir = iterChrom.current_dir()
            nb = 0
            sys.stderr.write("filtering splice graphs in %s\n" % current_dir)
