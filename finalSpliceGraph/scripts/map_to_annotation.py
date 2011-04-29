#!/usr/bin/python

"""
map_to_annotation.py:

    map splice graphs to given annotation for all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph, Mapper
from optparse import OptionParser, make_option

import psyco
#psyco.log()
#psyco.profile()
psyco.full()

# default directory
#defaultDir = '/r100/burge/shared/splice_graphs/mm6/FlatFiles'
defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles'


# @@@@ remember to add  the new annotation type in Mapper.py in Tools directory.

# digest command line
opts = OptionParser(option_list=                                                                              \
                    [make_option('-d','--dir', dest='inDir', default=defaultDir,                              \
                                 help='input directory with SpliceGraph* dirs [default: %s]' % defaultDir),   \
                     make_option('-t','--type', dest='type', default='u133v2',                            \
                                 help='annotation type [default: u133v2]'),                               \
                     make_option('-m','--method', dest='method', default='overlap',                           \
                                 help='method for mapping [default: overlap]'),                               \
                     make_option('-p','--inprefix', dest='inprefix', default='SpliceGraphsFiltered_',         \
                                 help='iterate over dirs with this prefix [default: SpliceGraphsFiltered_]'), \
                     make_option('-s','--species', dest='species', default='hsa',                             \
                                 help='species [default: hsa]'),                                              \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    mapper  = Mapper.Mapper(configFile='/r100/burge/xiao/Tools/splice_graphs/configuration.txt',verbose=1, species=options.species)
    #mapper  = Mapper.Mapper(verbose=1)

    inDirs = [d for d in os.listdir(options.inDir) if d.startswith(options.inprefix)]

    for inDir in inDirs:
        chr = inDir.split(options.inprefix)[-1]
        mapper.mapAllInDir(indir   = "%s/%s" % (options.inDir, inDir),
                           outfile = "%s/Mappings_%s_%s.txt" % (options.inDir, options.type, chr),
                           type    = options.type,
                           method  = options.method)
