#!/usr/bin/python

"""
classify_special_event.py:

    classify special splicing event for all *.sg files for all chromosomes in a given directory

"""

import os
import sys
import SpliceGraph

import psyco
#psyco.log()
#psyco.profile()
psyco.full()

# classifyer methods
def find_relaxed_A5SS(sg):
    "look for internal exons with common 3'SS and alternative 5'SS --> list of lists"
    events = []

    if isinstance(sg, SpliceGraph.SpliceGraph):
            
        # look for an internal exon start (common upstream 3'SS)
        for up3ss in sg.allElements():

            if sg.is3ss(up3ss):
                # are there several alternative 5ss?
                alt5sites = [alt5ss for alt5ss in sg.downstreamConnectedElements(up3ss) if sg.is5ss(alt5ss)]
                if len(alt5sites) > 1:

                    # store the events: up3ss -> alt5sites (=common[down3ss]) -> down3ss
                    #gnId   nbAlt5SS        alt5SS  coverage        upstream3SS     downstream3SS
                    events.append([sg.name,
                                   str(len(alt5sites)),
                                   ','.join(alt5sites),
                                   ','.join([str(sg.coverage(up3ss,alt5ss)) for alt5ss in alt5sites]),
                                   up3ss,
                                   '-'
                                   ])
     
    else:
        sys.stderr.write("ERROR in find_relaxed_A5SS: invalid splice graph")

    return events

def find_relaxed_A3SS(sg):
    "look for internal exons with common 5'SS and alternative 3'SS --> list of lists"
    events = []

    if isinstance(sg, SpliceGraph.SpliceGraph):
            
        # look for an internal exon start (common downstream 5'SS)
        for down5ss in sg.allElements():

            if sg.is5ss(down5ss):
                # are there several alternative 3ss?
                alt3sites = [alt3ss for alt3ss in sg.upstreamConnectedElements(down5ss) if sg.is3ss(alt3ss)]
                if len(alt3sites) > 1:

                    # store the events: up3ss -> alt5sites (=common[down3ss]) -> down3ss
                    #gnId   nbAlt3SS        alt3SS  coverage        upstream5SS     downstream5SS
                    events.append([sg.name,
                                   str(len(alt3sites)),
                                   ','.join(alt3sites),
                                   ','.join([str(sg.coverage(alt3ss,down5ss)) for alt3ss in alt3sites]),
                                   '-',
                                   down5ss
                                   ])
     
    else:
        sys.stderr.write("ERROR in find_relaxed_A3SS: invalid splice graph")

    return events


############################################################
# global variables
############################################################
inDir    = '/r100/burge/shared/splice_graphs/mm6/FlatFiles'
#inDir    = '/r100/burge/shared/splice_graphs/hg17/FlatFiles'
#outfile  = 'mm6.relaxed_A5SS.txt'
outfile  = 'mm6.relaxed_A3SS.txt'
inprefix = 'SpliceGraphsFiltered_'
#classifyMethod = find_relaxed_A5SS
classifyMethod = find_relaxed_A3SS
############################################################

if not os.path.isdir(inDir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes(inDir, prefix=inprefix)
    print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)

    # open output file
    outFH = open(outfile, 'w')
    #gnId   nbAlt5SS        alt5SS  coverage        upstream3SS     downstream3SS
    #outFH.write('\t'.join(['gnId','nbAlt5SS','alt5SS','coverage','upstream3SS','downstream3SS']) + '\n')
    #gnId   nbAlt3SS        alt3SS  coverage        upstream5SS     downstream5SS
    outFH.write('\t'.join(['gnId','nbAlt3SS','alt3SS','coverage','upstream5SS','downstream5SS']) + '\n')

    # iterate over splice graphs
    sg = iterChrom.next_sg()
    current_dir = iterChrom.current_dir()
    nb = 0
    sys.stderr.write("classifying splice graphs in %s\n" % current_dir)
    while ( sg ):
        nb += 1
        sys.stderr.write("analyzing splice graph %d of %d (%s)\r" % (nb, iterChrom.number(), sg.name))

        events = classifyMethod(sg)
        if len(events) > 0:
            for event in events:
                outFH.write('\t'.join(event) + '\n')

        # get the next splice graph
        sg = iterChrom.next_sg()

        # begin new output file if a new splice graph directory was started
        if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
            current_dir = iterChrom.current_dir()
            nb = 0
            sys.stderr.write("classifying splice graphs in %s\n" % current_dir)

    outFH.close()
