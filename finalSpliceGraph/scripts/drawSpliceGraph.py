#!/usr/bin/python

"""
drawSpliceGraph:

    draw splice graph corresponding to a gene

"""

help_string = """
    drawSpliceGraph-%s
    Michael Stadler, stadler@mit.edu

    USAGE:

      -g str     : genome (e.g. hg17 [default], mm6)
      -i str     : gene id (e.g. chr21:220441176:+)
      -f file    : *.sg file (overrides -g and -i)
      -e str     : emphasize element type (e.g. SE)
      -o file    : output filename

      \n""" % "$Revision: 1.4 $"

import os
import sys
import SpliceGraph
import Classifier
from optparse import OptionParser

# set the correct path to directories of .sg files
directories = {"hg17" : '/r100/burge/shared/splice_graphs/hg17/FlatFiles',
               #"hg16" : '/r100/burge/shared/splice_graphs/hg16/FlatFiles',
               #"mm5"  : '/r100/burge/shared/splice_graphs/mm5/FlatFiles',
               "mm6"  : '/r100/burge/shared/splice_graphs/mm6/FlatFiles'}

# set functions for possible emphasize elements
emphasize = {'SE'   : 'classify_SE',   \
             'A3SS' : 'classify_A3SS', \
             'A5SS' : 'classify_A5SS', \
             'CE'   : 'classify_CE',   \
             'CI'   : 'classify_CI',   \
             }

# digest command line
opts = OptionParser()
opts.add_option('-g','--genome', dest='genome')
opts.add_option('-i','--id', dest='id') 
opts.add_option('-f','--infile', dest='infile') 
opts.add_option('-e','--emph',  dest='emph')    
opts.add_option('-o','--outfile',  dest='outfile')    
(options, args) = opts.parse_args()

# check arguments
if options.genome is None:
    options.genome = 'hg17'

if not directories.has_key(options.genome):
    print >> sys.stderr, "ERROR: unknown genome %s, should be one of %s" % \
          (options.genome, ",".join(directories.keys()))
    sys.exit(0)

if options.emph is not None and not emphasize.has_key(options.emph):
    print >> sys.stderr, "ERROR: unknown Classifier element %s, should be one of %s" % \
          (options.emph, ",".join(emphasize.keys()))
    sys.exit(0)

if (options.id is None and options.infile is None) or options.outfile is None:
    sys.stdout.write(help_string)
else:
    inname = ""
    if options.infile is not None:
        inname = options.infile
    else:
        (chr, coord, strand) = options.id.split(':')
        inname = "%s/SpliceGraphs_%s/%s.sg" % (directories[options.genome], chr, options.id)
    try:
        sg = SpliceGraph.SpliceGraph(filename=inname)
    except:
        print >> sys.stderr, 'ERROR creating SpiceGraph object: %s' % sys.exc_info()[1]
    else:
        emphElements = None

        if options.emph is not None:
            c = Classifier.Classifier()
            emphMethod = getattr(c, emphasize[options.emph])
            emphElements = emphMethod(sg)

            if options.emph == 'A5SS':
                tmpDict = {}
                for key in emphElements.keys():
                    if not tmpDict.has_key(emphElements[key]["Upstream3ss"]):
                        tmpDict[emphElements[key]["Upstream3ss"]] = {}
                    for altSS in key.split(','):
                        tmpDict[emphElements[key]["Upstream3ss"]][altSS] = 1
                emphElements = tmpDict
            elif options.emph == 'A3SS':
                tmpDict = {}
                for key in emphElements.keys():
                    for altSS in key.split(','):
                        if not tmpDict.has_key(altSS):
                            tmpDict[altSS] = {}
                        tmpDict[altSS][emphElements[key]["Downstream5ss"]] = 1
                emphElements = tmpDict

            for i in emphElements.keys():
                for j in emphElements[i].keys():
                    emphElements[i][j] = options.emph #+ str(emphElements[i][j])

        sg.visualize(filename=options.outfile, emph=emphElements)
