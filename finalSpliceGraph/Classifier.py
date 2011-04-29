#---------------------------------------------------------------------------------------
#   $Id: Classifier.py,v 1.34 2006/03/01 22:10:44 stadler Exp $ 
#---------------------------------------------------------------------------------------
"""Framework for classifying constitutive/alternative exons and introns based on splice graphs

Synopsis:
    import Classifyer

    c = Classifier.Classifier()

    list1 = c.classify_SE('file1.sg')

    import SpliceGraph
    sg = SpliceGraph.SpliceGraph(filename='file2.sg')

    list2 = c.classify_SE(sg)

"""

import sys, os, configuration, Errors, SpliceGraph

class Classifier:
    "Classifier: Object for classifying events based on splice graphs"

    def __init__(self, conf=None, verbose=False):
        "Classifier: constructor"

        self.verbose=verbose
        
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        else:
            self.conf = configuration.Configuration(verbose=self.verbose)


    def classify_SE(self, sg):
        "identify SE (skipped exons) from splice graph"

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an internal exon start (3ss element)
            for exStart in sg.allElements():

                if sg.is3ss(exStart):

                    # store start(s) of upstream intron(s)
                    upConnected = sg.upstreamConnectedElements(exStart)

                    # look for a corresponding exon end (5ss element)
                    for exEnd in sg.downstreamConnectedElements(exStart):

                        # make sure exEnd is a 5ss element
                        if sg.is5ss(exEnd):

                            nbIncl = sg.areDirectlyConnected(exStart, exEnd)
                            
                            # store end(s) of downstream intron(s)
                            downConnected = sg.downstreamConnectedElements(exEnd)

                            # check for the skipping event(s): connection between upConnected and downConnected
                            nbSkip = 0
                            for up in upConnected:
                                for down in downConnected:
                                    nbSkip += abs(sg.areDirectlyConnected(up, down))

                            if nbSkip > 0:
                                # store the exon
                                if not events.has_key(exStart):
                                    events[exStart] = {}
                                events[exStart][exEnd] = (nbSkip, nbIncl)
                    
        else:
            self.__inform("ERROR in classify_SE: invalid splice graph")
            raise Errors.ArgumentError("Classify_SE", "argument must either be filename or instance of SpliceGraph")

        return events

    def classify_SEP(self, sg):
        "identify SE (skipped exons) from splice graph, ADDED UPSTREAM AND DOWNSTREAM EXONS TO OUTPUT FILE"

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an internal exon start (3ss element)
            for exStart in sg.allElements():

                if sg.is3ss(exStart):

                    # store start(s) of upstream intron(s)
                    upConnected = sg.upstreamConnectedElements(exStart)
                    
                    #!!!!CG HACK
                    upExList = []
                    for upExStart in upConnected:
                    	for el in sg.upstreamConnectedElements(upExStart):
                    		upExList.append('%s_%s' % (upExStart, el))

                    # look for a corresponding exon end (5ss element)
                    for exEnd in sg.downstreamConnectedElements(exStart):

                        # make sure exEnd is a 5ss element
                        if sg.is5ss(exEnd):

                            nbIncl = sg.areDirectlyConnected(exStart, exEnd)
                            
                            # store end(s) of downstream intron(s)
                            downConnected = sg.downstreamConnectedElements(exEnd)
                            
                            #!!!!CG HACK
                            downExList = []
                            for downExStart in downConnected:
                            	for el in sg.downstreamConnectedElements(downExStart):
                            		downExList.append('%s_%s' % (downExStart, el))                            			

                            # check for the skipping event(s): connection between upConnected and downConnected
                            nbSkip = 0
                            for up in upConnected:
                                for down in downConnected:
                                    nbSkip += abs(sg.areDirectlyConnected(up, down))

                            if nbSkip > 0:
                                # store the exon
                                if not events.has_key(exStart):
                                    events[exStart] = {}
                                events[exStart][exEnd] = (nbSkip, nbIncl, ','.join(upExList), ','.join(downExList))
                    
        else:
            self.__inform("ERROR in classify_SE: invalid splice graph")
            raise Errors.ArgumentError("Classify_SE", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_RI(self, sg):
        "identify RI (retained introns) from splice graph"

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an internal intron start (5ss element)
            for inStart in sg.allElements():

                if sg.is5ss(inStart):

                    # store start(s) of upstream exon(s)
                    upConnected = sg.upstreamConnectedElements(inStart)

                    # look for a corresponding intron end (3ss element)
                    for inEnd in sg.downstreamConnectedElements(inStart):

                        # make sure exEnd is a 3ss element
                        if sg.is3ss(inEnd):

                            nbSpliced = sg.areDirectlyConnected(inStart, inEnd)
                            
                            # store end(s) of downstream exon(s)
                            downConnected = sg.downstreamConnectedElements(inEnd)

                            # check for the retaining event(s): connection between upConnected and downConnected
                            nbRetained = 0
                            for up in upConnected:
                                for down in downConnected:
                                    nbRetained += abs(sg.areDirectlyConnected(up, down))

                            if nbRetained > 0:
                                # store the event
                                if not events.has_key(inStart):
                                    events[inStart] = {}
                                events[inStart][inEnd] = (nbRetained, nbSpliced)
                    
        else:
            self.__inform("ERROR in classify_RI: invalid splice graph")
            raise Errors.ArgumentError("Classify_RI", "argument must either be filename or instance of SpliceGraph")

        return events
    
    def classify_RIP(self, sg):
        "identify RI (retained introns) from splice graph, ADDED UPSTREAM AND DOWNSTREAM EXONS TO OUTPUT FILE"

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an internal intron start (5ss element)
            for inStart in sg.allElements():

                if sg.is5ss(inStart):

                    # store start(s) of upstream exon(s)
                    upConnected = sg.upstreamConnectedElements(inStart)
                    
                    #!!!!CG HACK
                    upExList = []
                    for upExStart in upConnected:
                    	for el in sg.upstreamConnectedElements(upExStart):
                    		upExList.append('%s_%s' % (upExStart, el))

                    # look for a corresponding intron end (3ss element)
                    for inEnd in sg.downstreamConnectedElements(inStart):

                        # make sure exEnd is a 3ss element
                        if sg.is3ss(inEnd):

                            nbSpliced = sg.areDirectlyConnected(inStart, inEnd)
                            
                            # store end(s) of downstream exon(s)
                            downConnected = sg.downstreamConnectedElements(inEnd)
                            
                            #!!!!CG HACK
                            downExList = []
                            for downExStart in downConnected:
                            	for el in sg.downstreamConnectedElements(downExStart):
                            		downExList.append('%s_%s' % (downExStart, el))  

                            # check for the retaining event(s): connection between upConnected and downConnected
                            nbRetained = 0
                            for up in upConnected:
                                for down in downConnected:
                                    nbRetained += abs(sg.areDirectlyConnected(up, down))

                            if nbRetained > 0:
                                # store the event
                                if not events.has_key(inStart):
                                    events[inStart] = {}
                                events[inStart][inEnd] = (nbRetained, nbSpliced, ','.join(upExList), ','.join(downExList))
                    
        else:
            self.__inform("ERROR in classify_RI: invalid splice graph")
            raise Errors.ArgumentError("Classify_RI", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_A5SSobsolete(self, sg):
        """identify alternative 5'splice site usage from splice graph, where the two
        alternative exons feed both into the same downstream exon"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # look for an internal exon start (3ss element)
            for exStart in sg.allElements():
                
                if sg.is3ss(exStart):

                    # dictionary to store 3'ss -> downstream 5'ss links, temporarily
                    exEnd_and_downConnected = {}
                    exEnd_and_downConnected["exEnds"] = []
                    exEnd_and_downConnected["downstreamExStarts"] = []

                    # look for corresponding exon end (5ss elements)
                    for exEnd in sg.downstreamConnectedElements(exStart):
                        
                        # make sure exEnd is a 5ss element
                        if sg.is5ss(exEnd):
                            downConnected = sg.downstreamConnectedElements(exEnd)
                            # store 5'ss and the downstream 3'ss
                            exEnd_and_downConnected["exEnds"].append(exEnd)
                            for downstreamExStart in downConnected:
                                exEnd_and_downConnected["downstreamExStarts"].append(downstreamExStart)
                            

                    # see if there are two, or more, downstream 5'ss
                    if len( exEnd_and_downConnected["exEnds"] ) > 1:
                        exEnds = exEnd_and_downConnected["exEnds"]
                        nrExEnds = len(exEnds)
                        downstreamExStarts = exEnd_and_downConnected["downstreamExStarts"]
                                               
                        for i in xrange( nrExEnds ):
                            for j in xrange(i+1 , nrExEnds, 1):
                                # see if their path's converge
                                for downstreamExStart in downstreamExStarts:
                                    if sg.areDirectlyConnected(exEnds[i],downstreamExStart) and sg.areDirectlyConnected(exEnds[j],downstreamExStart):
                                        # save alternative 5'ss event
                                        new_key = '%s,%s' % (exEnds[i],exEnds[j])
                                        events[new_key] = {}
                                        events[new_key]["Upstream3ss"] = exStart
                                        events[new_key]["Downstream3ss"] = downstreamExStart
                                        events[new_key]["Coverage"] = (abs(sg.areDirectlyConnected(exEnds[i],downstreamExStart)),
                                                                       abs( sg.areDirectlyConnected(exEnds[j],downstreamExStart)))  
     
        else:
            self.__inform("ERROR in classify_A5SS: invalid splice graph")
            raise Errors.ArgumentError("Classify_A5SS", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_A5SS(self, sg):
        """identify alternative 5'splice site usage from splice graph, defined as:
        two or more 5'SS with common upstream and downstream 3'SS

        This method is a complete rewrite that avoids creating of redundant events in cases of more than two 5'SS

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # look for an internal exon start (common upstream 3'SS)
            for up3ss in sg.allElements():

                if sg.is3ss(up3ss):
                    # are there several alternative 5ss?
                    alt5sites = [alt5ss for alt5ss in sg.downstreamConnectedElements(up3ss) if sg.is5ss(alt5ss)]
                    if len(alt5sites) > 1:

                        # collect alt5sites that have a common downstream 3'SS
                        common = {} # format: common{down3ss} = [alt5site1, alt5site2, ...]
                        for alt5site in alt5sites:
                            for down3ss in [down3ss for down3ss in sg.downstreamConnectedElements(alt5site) if sg.is3ss(down3ss)]:
                                if not common.has_key(down3ss):
                                    common[down3ss] = []
                                common[down3ss].append(alt5site)
                         
                        # store the events: up3ss -> alt5sites (=common[down3ss]) -> down3ss
                        for down3ss in common.iterkeys():
                            if len(common[down3ss]) > 1:
                                new_key = ','.join(common[down3ss])
                                events[new_key] = {}
                                events[new_key]["Upstream3ss"] = up3ss
                                events[new_key]["Downstream3ss"] = down3ss
                                events[new_key]["Coverage"] = [ abs(sg.areDirectlyConnected(alt5site,down3ss)) \
                                                                for alt5site in common[down3ss] ]  
     
        else:
            self.__inform("ERROR in classify_A5SS: invalid splice graph")
            raise Errors.ArgumentError("Classify_A5SS", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_A3SSobsolete(self, sg):
        """identify alternative 3'splice site usage from splice graph, where the two alternative exons feed both into the same downstream exon"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # look for an internal exon start (3ss element)
            for exEnd in sg.allElements():
                
                if sg.is5ss(exEnd):

                    # dictionary to store 3'ss upstream of the 5'ss, temporarily
                    exEnd_and_downConnected = {}
                    exEnd_and_downConnected["exStarts"] = []
                    exEnd_and_downConnected["upstreamExEnds"] = []

                    # look for corresponding exon end (5ss elements)
                    for exStart in sg.upstreamConnectedElements(exEnd):
                        
                        # make sure exEnd is a 5ss element
                        if sg.is3ss(exStart):
                            upConnected = sg.upstreamConnectedElements(exStart)
                            # store 5'ss and the downstream 3'ss
                            exEnd_and_downConnected["exStarts"].append(exStart)
                            for upstreamExEnd in upConnected:
                                exEnd_and_downConnected["upstreamExEnds"].append(upstreamExEnd)
                            

                    # see if there are two, or more, upstream 3'ss
                    if len( exEnd_and_downConnected["exStarts"] ) > 1:
                        exStarts = exEnd_and_downConnected["exStarts"]
                        nrExStarts = len(exStarts)
                        upstreamExEnds = exEnd_and_downConnected["upstreamExEnds"]

                        for i in xrange( nrExStarts ):
                            for j in xrange(i+1 , nrExStarts, 1):
                                # see if their path's converge
                                for upstreamExEnd in upstreamExEnds:
                                                                     
                                    if sg.areDirectlyConnected(upstreamExEnd,exStarts[i]) and sg.areDirectlyConnected(upstreamExEnd,exStarts[j]):
                                        # save alternative 5'ss event
  
                                        new_key = '%s,%s' % (exStarts[i],exStarts[j])
                                        events[new_key] = {}
                                        events[new_key]["Upstream5ss"] = upstreamExEnd
                                        events[new_key]["Downstream5ss"] = exEnd
                                        events[new_key]["Coverage"] = (abs(sg.areDirectlyConnected(upstreamExEnd,exStarts[i])),
                                                                       abs( sg.areDirectlyConnected(upstreamExEnd,exStarts[j])))  
     
        else:
            self.__inform("ERROR in classify_A3SS: invalid splice graph")
            raise Errors.ArgumentError("Classify_A3SS", "argument must either be filename or instance of SpliceGraph")

        return events
        

    def classify_A3SS(self, sg):
        """identify alternative 3'splice site usage from splice graph, defined as:
        two or more 3'SS with common upstream and downstream 5'SS

        This method is a complete rewrite that avoids creating of redundant events in cases of more than two 3'SS

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # look for an internal exon end (common upstream 5'SS)
            for up5ss in sg.allElements():

                if sg.is5ss(up5ss):
                    # are there several alternative 3ss?
                    alt3sites = [alt3ss for alt3ss in sg.downstreamConnectedElements(up5ss) if sg.is3ss(alt3ss)]
                    if len(alt3sites) > 1:

                        # collect alt3sites that have a common downstream 5'SS
                        common = {} # format: common{down5ss} = [alt3site1, alt3site2, ...]
                        for alt3site in alt3sites:
                            for down5ss in [down5ss for down5ss in sg.downstreamConnectedElements(alt3site) if sg.is5ss(down5ss)]:
                                if not common.has_key(down5ss):
                                    common[down5ss] = []
                                common[down5ss].append(alt3site)
                         
                        # store the events: up5ss -> alt3sites (=common[down5ss]) -> down5ss
                        for down5ss in common.iterkeys():
                            if len(common[down5ss]) > 1:
                                new_key = ','.join(common[down5ss])
                                events[new_key] = {}
                                events[new_key]["Upstream5ss"] = up5ss
                                events[new_key]["Downstream5ss"] = down5ss
                                events[new_key]["Coverage"] = [ abs(sg.areDirectlyConnected(alt3site,down5ss)) \
                                                                for alt3site in common[down5ss] ]  
     
        else:
            self.__inform("ERROR in classify_A3SS: invalid splice graph")
            raise Errors.ArgumentError("Classify_A3SS", "argument must either be filename or instance of SpliceGraph")

        return events
    
    def classify_A3SSP(self, sg):
        """identify alternative 3'splice site usage from splice graph, defined as:
        two or more 3'SS with common upstream and downstream 5'SS

        This method is a complete rewrite that avoids creating of redundant events in cases of more than two 3'SS

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # look for an internal exon end (common upstream 5'SS)
            for up5ss in sg.allElements():

                if sg.is5ss(up5ss):

                    upExList = []
                    for el in sg.upstreamConnectedElements(up5ss):
                        upExList.append('%s_%s' % (up5ss, el))

                    # are there several alternative 3ss?
                    alt3sites = [alt3ss for alt3ss in sg.downstreamConnectedElements(up5ss) if sg.is3ss(alt3ss)]
                    if len(alt3sites) > 1:

                        # collect alt3sites that have a common downstream 5'SS
                        common = {} # format: common{down5ss} = [alt3site1, alt3site2, ...]
                        for alt3site in alt3sites:
                            for down5ss in [down5ss for down5ss in sg.downstreamConnectedElements(alt3site) if sg.is5ss(down5ss)]:
                                if not common.has_key(down5ss):
                                    common[down5ss] = []
                                common[down5ss].append(alt3site)
                         
                        # store the events: up5ss -> alt3sites (=common[down5ss]) -> down5ss
                        for down5ss in common.iterkeys():
                            if len(common[down5ss]) > 1:
                                new_key = ','.join(common[down5ss])
                                events[new_key] = {}
                                events[new_key]["Upstream5ss"] = up5ss
                                events[new_key]["upExons"] = ','.join(upExList)
                                events[new_key]["Downstream5ss"] = down5ss
                                events[new_key]["Coverage"] = [ abs(sg.areDirectlyConnected(alt3site,down5ss)) \
                                                                for alt3site in common[down5ss] ]  
     
        else:
            self.__inform("ERROR in classify_A3SS: invalid splice graph")
            raise Errors.ArgumentError("Classify_A3SS", "argument must either be filename or instance of SpliceGraph")

        return events

    def classify_CE_obsolete(self, sg, intronStrict=False):
        """Identify Constitutive exons, defined as: an exon edge that if removed destroys any path from upstream to downstream elements
        intronStrict : if true, there have to be single flanking introns for a CE

        OBSOLETE: do not use any more
                  this method is not efficient for large splice graphs and does not give correct results in cases where
                  there is an overlap between the putative CE and the last edge in a parallel but fragmentary path
                  whose final node is inside of CE

        --> dict
        """
        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            # do the classification
            
            # look for an internal exon start (3ss element)
            for exStart in sg.allElements():
                
                if sg.is3ss(exStart):
                    for exEnd in sg.downstreamConnectedElements(exStart):
                        
                        # make sure exEnd is a 5ss element
                        if sg.is5ss(exEnd):

                            # is it the only edge between exStart,exEnd?
                            if len( sg.upstreamConnectedElements(exEnd) )     == 1 and \
                               len( sg.downstreamConnectedElements(exStart) ) == 1:

                                if (intronStrict is True and \
                                    len( sg.upstreamConnectedElements(exStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(exEnd) ) == 1 ) or \
                                    intronStrict is False:

                                    # if intronStrict is True:
                                    # this is an exon for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    CE = True

                                    if sg.ListsAreConnectedNotVia(sg.upstreamElements(exStart),
                                                                  sg.downstreamElements(exEnd),
                                                                  exStart, exEnd):
                                        CE = False

                                    if CE is True:
                                        events[exStart] = {}
                                        events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)
        else:
            self.__inform("ERROR in classify_CE: invalid splice graph")
            raise Errors.ArgumentError("Classify_CE", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_CE(self, sg, intronStrict=False):
        """Identify Constitutive exons, defined as:
        exon edges that do not overlap any other exon/intron edges in the current splice graph
        intronStrict : if true, there have to be single flanking introns for a CE
        --> dict
        """
        events = {} # format: events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            if sg.nbExons() > 0:
                # construct coverage vector
                totStart, totEnd = sg.genomicRange()
                coverage = []
                for i in xrange(totEnd - totStart + 1):
                    coverage.append(0)

                elements = sg.allElements()

                for startEl,endEl in sg.allExons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1
                for startEl,endEl in sg.allIntrons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    startCoord += 1 #do not overcount first/last base
                    endCoord   -= 1
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1

                # look for an internal exon start (3ss element)
                for exStart in sg.allElements():
                
                    if sg.is3ss(exStart):
                        for exEnd in sg.downstreamConnectedElements(exStart):
                        
                            # make sure there is only one connection exStart,exEnd and exEnd is a 5ss element
                            if len( sg.upstreamConnectedElements(exEnd) )     == 1 and \
                               len( sg.downstreamConnectedElements(exStart) ) == 1 and \
                               sg.is5ss(exEnd):

                                if (intronStrict is True and \
                                    len( sg.upstreamConnectedElements(exStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(exEnd) ) == 1 ) or \
                                    intronStrict is False:

                                    # if intronStrict is True:
                                    # this is an exon for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    startCoord = int(exStart.split(':')[1])
                                    endCoord   = int(exEnd.split(':')[1])
                                    if endCoord < startCoord:
                                        startCoord, endCoord = endCoord, startCoord
                                    pos = startCoord - totStart
                                    CE = True
                                    while pos <= endCoord - totStart:
                                        if coverage[pos] != 1:
                                            CE = False
                                            break
                                        pos += 1
                                    if CE:
                                        events[exStart] = {}
                                        events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)
        else:
            self.__inform("ERROR in classify_CE: invalid splice graph")
            raise Errors.ArgumentError("Classify_CE", "argument must either be filename or instance of SpliceGraph")

        return events
        
    def classify_CEP(self, sg, intronStrict=False):
        """Identify Constitutive exons, defined as:
        exon edges that do not overlap any other exon/intron edges in the current splice graph
        intronStrict : if true, there have to be single flanking introns for a CE
        --> dict
        
        !!! give list of directly connected upstream and downstream elements as well
        """
        events = {} # format: events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            if sg.nbExons() > 0:
                # construct coverage vector
                totStart, totEnd = sg.genomicRange()
                coverage = []
                for i in xrange(totEnd - totStart + 1):
                    coverage.append(0)

                elements = sg.allElements()

                for startEl,endEl in sg.allExons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1
                for startEl,endEl in sg.allIntrons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    startCoord += 1 #do not overcount first/last base
                    endCoord   -= 1
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1

                # look for an internal exon start (3ss element)
                for exStart in sg.allElements():
                
                    if sg.is3ss(exStart):
                        for exEnd in sg.downstreamConnectedElements(exStart):
                        
                            # make sure there is only one connection exStart,exEnd and exEnd is a 5ss element
                            if len( sg.upstreamConnectedElements(exEnd) )     == 1 and \
                               len( sg.downstreamConnectedElements(exStart) ) == 1 and \
                               sg.is5ss(exEnd):

                                if (intronStrict is True and \
                                    len( sg.upstreamConnectedElements(exStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(exEnd) ) == 1 ) or \
                                    intronStrict is False:

                                    # if intronStrict is True:
                                    # this is an exon for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    startCoord = int(exStart.split(':')[1])
                                    endCoord   = int(exEnd.split(':')[1])
                                    if endCoord < startCoord:
                                        startCoord, endCoord = endCoord, startCoord
                                    pos = startCoord - totStart
                                    CE = True
                                    while pos <= endCoord - totStart:
                                        if coverage[pos] != 1:
                                            CE = False
                                            break
                                        pos += 1
                                    if CE:
                                        events[exStart] = {}
                                    	events[exStart][exEnd] = [] #this will contain the three things
                                        events[exStart][exEnd].append(sg.areDirectlyConnected(exStart,exEnd))
                                        events[exStart][exEnd].append(','.join(sg.upstreamConnectedElements(exStart)))
                                        events[exStart][exEnd].append(','.join(sg.downstreamConnectedElements(exEnd)))
        else:
            self.__inform("ERROR in classify_CE: invalid splice graph")
            raise Errors.ArgumentError("Classify_CE", "argument must either be filename or instance of SpliceGraph")

        return events



    def classify_CE_graphBased(self, sg):
        """Identify Constitutive exons, defined as: single exon edges that connect two articulation points (cut vertices) --> dict
        """
        events = {} # format: events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            # get a list of bridge edges (when removed, will disconnect the graph. bridge edges connect articulation points)
            bridge_tuples = sg.bridgeEdges()
            bridges = {}
            for b in bridge_tuples:
                if not bridges.has_key(b[0]):
                    bridges[b[0]] = {}
                bridges[b[0]][b[1]] = 1
            
            # look for an internal exon start (3ss element)
            for exStart in sg.allElements():
                
                if sg.is3ss(exStart):
                    for exEnd in sg.downstreamConnectedElements(exStart):
                        
                        # make sure exEnd is the only element and a 5ss element
                        if len( sg.downstreamConnectedElements(exStart) ) == 1 and sg.is5ss(exEnd):

                            # this is an exon
                            if bridges.has_key(exStart) and bridges[exStart].has_key(exEnd):
                                events[exStart] = {}
                                events[exStart][exEnd] = sg.areDirectlyConnected(exStart,exEnd)
        else:
            self.__inform("ERROR in classify_CE: invalid splice graph")
            raise Errors.ArgumentError("Classify_CE", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_CI_obsolete(self, sg, exonStrict=False):
        """Identify Constitutive introns, defined as:
           exonStrict : if true, there have to be single flanking exons for a CE

        OBSOLETE: do not use any more
                  this method is not efficient for large splice graphs and does not give correct results in cases where
                  there is an overlap between the putative CI and the last edge in a parallel but fragmentary path
                  whose final node is inside of CI

        --> dict
        """
        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            # do the classification
            
            # look for an internal exon start (5ss element)
            for inStart in sg.allElements():
                
                if sg.is5ss(inStart):
                    for inEnd in sg.downstreamConnectedElements(inStart):
                        
                        # make sure exEnd is a 3ss element
                        if sg.is3ss(inEnd):

                            # this is an intron
                            if len( sg.downstreamConnectedElements(inStart) ) == 1:
                                
                                if (exonStrict is True and \
                                    len( sg.upstreamConnectedElements(inStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(inEnd) ) == 1 ) or \
                                    exonStrict is False:

                                    # if exonStrict is True:
                                    # this is an intron for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    CI = True

                                    if sg.ListsAreConnectedNotVia(sg.upstreamElements(inStart),
                                                                  sg.downstreamElements(inEnd),
                                                                  inStart, inEnd):
                                        CI = False

                                    if CI is True:
                                        events[inStart] = {}
                                        events[inStart][inEnd] = sg.areDirectlyConnected(inStart, inEnd)
        else:
            self.__inform("ERROR in classify_CI: invalid splice graph")
            raise Errors.ArgumentError("Classify_CI", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_CI(self, sg, exonStrict=False):
        """Identify Constitutive introns, defined as:
        intron edges that do not overlap any other exon/intron edges in the current splice graph
        exonStrict : if true, there have to be single flanking exons for a CI
        --> dict
        """
        events = {} # format: events[inStart][inEnd] = sg.areDirectlyConnected(inStart,inEnd)

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            if sg.nbIntrons() > 0:
                # construct coverage vector
                totStart, totEnd = sg.genomicRange()
                coverage = []
                for i in xrange(totEnd - totStart + 1):
                    coverage.append(0)

                elements = sg.allElements()

                for startEl,endEl in sg.allExons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1
                for startEl,endEl in sg.allIntrons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    startCoord += 1 #do not overcount first/last base
                    endCoord   -= 1
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1

                # look for an internal intron start (5ss element)
                for inStart in sg.allElements():

                    if sg.is5ss(inStart):
                        for inEnd in sg.downstreamConnectedElements(inStart):

                            # make sure there is only one connection inStart,inEnd and inEnd is a 3ss element
                            if len( sg.upstreamConnectedElements(inEnd) )     == 1 and \
                               len( sg.downstreamConnectedElements(inStart) ) == 1 and \
                               sg.is3ss(inEnd):

                                if (exonStrict is True and \
                                    len( sg.upstreamConnectedElements(inStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(inEnd) ) == 1 ) or \
                                    exonStrict is False:

                                    # if exonStrict is True:
                                    # this is an intron for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    startCoord = int(inStart.split(':')[1])
                                    endCoord   = int(inEnd.split(':')[1])
                                    if endCoord < startCoord:
                                        startCoord, endCoord = endCoord, startCoord
                                    pos = startCoord - totStart
                                    CI = True
                                    while pos <= endCoord - totStart:
                                        if coverage[pos] != 1:
                                            CI = False
                                            break
                                        pos += 1
                                    if CI:
                                        events[inStart] = {}
                                        events[inStart][inEnd] = sg.areDirectlyConnected(inStart,inEnd)
        else:
            self.__inform("ERROR in classify_CI: invalid splice graph")
            raise Errors.ArgumentError("Classify_CI", "argument must either be filename or instance of SpliceGraph")

        return events
        
    def classify_CIP(self, sg, exonStrict=False):
        """Identify Constitutive introns, defined as:
        intron edges that do not overlap any other exon/intron edges in the current splice graph
        exonStrict : if true, there have to be single flanking exons for a CI
        --> dict
        """
        events = {} # format: events[inStart][inEnd] = sg.areDirectlyConnected(inStart,inEnd)

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            if sg.nbIntrons() > 0:
                # construct coverage vector
                totStart, totEnd = sg.genomicRange()
                coverage = []
                for i in xrange(totEnd - totStart + 1):
                    coverage.append(0)

                elements = sg.allElements()

                for startEl,endEl in sg.allExons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1
                for startEl,endEl in sg.allIntrons():
                    startCoord = int(startEl.split(':')[1])
                    endCoord   = int(endEl.split(':')[1])
                    if endCoord < startCoord:
                        startCoord, endCoord = endCoord, startCoord
                    startCoord += 1 #do not overcount first/last base
                    endCoord   -= 1
                    pos        = startCoord - totStart
                    while pos <= endCoord - totStart:
                        coverage[pos] += 1
                        pos += 1

                # look for an internal intron start (5ss element)
                for inStart in sg.allElements():

                    if sg.is5ss(inStart):
                        for inEnd in sg.downstreamConnectedElements(inStart):

                            # make sure there is only one connection inStart,inEnd and inEnd is a 3ss element
                            if len( sg.upstreamConnectedElements(inEnd) )     == 1 and \
                               len( sg.downstreamConnectedElements(inStart) ) == 1 and \
                               sg.is3ss(inEnd):

                                if (exonStrict is True and \
                                    len( sg.upstreamConnectedElements(inStart) ) == 1 and \
                                    len( sg.downstreamConnectedElements(inEnd) ) == 1 ) or \
                                    exonStrict is False:

                                    # if exonStrict is True:
                                    # this is an intron for which no alternative splicing occurs at its
                                    # 5' nor 3' splice site

                                    startCoord = int(inStart.split(':')[1])
                                    endCoord   = int(inEnd.split(':')[1])
                                    if endCoord < startCoord:
                                        startCoord, endCoord = endCoord, startCoord
                                    pos = startCoord - totStart
                                    CI = True
                                    while pos <= endCoord - totStart:
                                        if coverage[pos] != 1:
                                            CI = False
                                            break
                                        pos += 1
                                    if CI:
                                        events[inStart] = {}
                                        events[inStart][inEnd] = []
                                        events[inStart][inEnd].append(sg.areDirectlyConnected(inStart,inEnd))
                                        events[inStart][inEnd].append(','.join(sg.upstreamConnectedElements(inStart)))
                                        events[inStart][inEnd].append(','.join(sg.downstreamConnectedElements(inEnd)))
        else:
            self.__inform("ERROR in classify_CI: invalid splice graph")
            raise Errors.ArgumentError("Classify_CI", "argument must either be filename or instance of SpliceGraph")

        return events


    def classify_MXEobsolete(self, sg):
        "Identifies mutually exclusive exons in a splice graph"

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            exons = sg.allExons()

            for exon1 in exons:
                for exon2 in exons:

                    if exon1 == exon2:
                        continue

                    (ex1Start, ex1End) = exon1
                    (ex2Start, ex2End) = exon2

                    if ex1Start == ex2Start or ex1End == ex2End:
                        continue

                    # both exons should splice to the same upstream and downstream ss
                    # and to no other splice sites (alternative splice site usage is not allowed)
                    upConnected1 = sg.upstreamConnectedElements(ex1Start)
                    upConnected2 = sg.upstreamConnectedElements(ex2Start)

                    if len(upConnected1) != 1 or len(upConnected2) != 1:
                        continue
                    if upConnected1 != upConnected2:
                        continue

                    upEx5ss = upConnected1[0]

                    downConnected1 = sg.downstreamConnectedElements(ex1End)
                    downConnected2 = sg.downstreamConnectedElements(ex2End)

                    if len(downConnected1) != 1 or len(downConnected2) != 1:
                        continue
                    if downConnected1 != downConnected2:
                        continue

                    downEx3ss = downConnected1[0]

                    # the upstream 5ss and downstream 3ss should not have other splices than the two exons under investigation
                    # except for a possible direct splice between them

                    maxNb = 2
                    if sg.areDirectlyConnected(upEx5ss, downEx3ss):
                        maxNb += 1

                    if len(sg.downstreamConnectedElements(upEx5ss)) > maxNb or len(sg.upstreamConnectedElements(downEx3ss)) > maxNb:
                        continue

                    # there should not be any links between the exons, i.e. they should be mutually exclusive
                    if sg.areDirectlyConnected(ex1Start,ex2End) or sg.areDirectlyConnected(ex1End,ex2Start) or \
                       sg.areDirectlyConnected(ex2End,ex1Start) or sg.areDirectlyConnected(ex2Start,ex1End):
                        continue

                    # otherwise a MXE pair found
                    if events.has_key(upEx5ss):
                        # already reported
                        continue
                    
                    events[upEx5ss] = {}
                    events[upEx5ss]["exon1"] = (ex1Start, ex1End, sg.coverage(ex1Start, ex1End))
                    events[upEx5ss]["exon2"] = (ex2Start, ex2End, sg.coverage(ex2Start, ex2End))
                    events[upEx5ss]["bothSkipped"] = sg.coverage(upEx5ss, downEx3ss)
                    events[upEx5ss]["downstream3ss"] = downEx3ss
        return events

    def build_connection_matrix(self, sg):
        "build an matrix of connections between splice sites -> list of lists"
        # y-axis : 5'ss and TER
        # x-axis : 3'ss and TSS

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            nb_xElements = 0
            nb_yElements = 0
            mtx = []
            elementsMap = {}

            elements = sg.allElements()
            nrElements = len(elements)


            # get matrix size
            for element in elements:
                if sg.isTER(element) or sg.is5ss(element):
                    elementsMap[element] = nb_yElements
                    nb_yElements += 1
                    mtx.append([])
                else:
                    elementsMap[element] = nb_xElements
                    nb_xElements += 1

            # build matrix
            for y in xrange(nb_yElements):
                mtx[y] = [0*x for x in range(nb_xElements)]

            # fill in connections
            #for connection in sg.allConnections():
            for el1 in xrange(nrElements):
                for el2 in xrange(nrElements):
                    if sg.areDirectlyConnected(elements[el1],elements[el2]): 
                        if sg.isTER(elements[el1]) or sg.is5ss(elements[el1]):
                            mtx[elementsMap[elements[el1]]][elementsMap[elements[el2]]] -= 1
                        else:
                            mtx[elementsMap[elements[el2]]][elementsMap[elements[el1]]] += 1
                            
            for list in mtx:
                print list
            

                            
    def classify_MXE(self, sg):
        """identify mutually exclusive exons from splice graph, defined as:
        two or more non-overapping exons with common common upstream and downstream 5'SS

        This method is a complete rewrite that allows the presence of more than two MXE and enforces them to be non-overlapping

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an exon end (common upstream 5ss: upEx5ss)
            for upEx5ss in sg.allElements():

                if sg.is5ss(upEx5ss):
                    # are there several alternative 3ss?
                    MXEstarts = sg.downstreamConnectedElements(upEx5ss)
                    if len(MXEstarts) > 1:

                        # collect candidate MXE (common downstream 3ss: downEx3ss, non-overlapping)
                        candidates = {} # format: candidates[downEx3ss] = [ (start1, end1), (start2, end2), ... ]
                        for start in MXEstarts:
                            MXEends = sg.downstreamConnectedElements(start)
                            for end in MXEends:
                                for downEx3ss in sg.downstreamConnectedElements(end):
                                    if not candidates.has_key(downEx3ss):
                                        candidates[downEx3ss] = [ (start,end) ]
                                    else:
                                        # make sure they are non-overlapping and not co-occurring in same transcript
                                        keepThis = True
                                        for (otherStart, otherEnd) in candidates[downEx3ss]:
                                            if ( sg.strand() == '+' and \
                                                 int(end.split(':')[1])   >= int(otherStart.split(':')[1]) and \
                                                 int(start.split(':')[1]) <= int(otherEnd.split(':')[1])           ) or \
                                               ( sg.strand() == '-' and \
                                                 int(start.split(':')[1]) >= int(otherEnd.split(':')[1]) and \
                                                 int(end.split(':')[1])   <= int(otherStart.split(':')[1])         ):
                                                # overlaps existing MXE candidate
                                                keepThis = False
                                                break
                                            if sg.areDirectlyConnected(end, otherStart) or sg.areDirectlyConnected(otherEnd, start):
                                                # is present in the same transcript as existing MXE candidate
                                                keepThis = False
                                                break
                                        if keepThis:
                                            candidates[downEx3ss].append( (start,end) )

                        # store the events: upEx5ss -> MXE (=candidates[downEx3ss]) -> downEx3ss
                        for downEx3ss in candidates.iterkeys():
                            if len(candidates[downEx3ss]) > 1:
                                new_key = ','.join(['%s;%s' % (start,end) for start,end in candidates[downEx3ss] ])
                                events[new_key] = {}
                                events[new_key]["Upstream5ss"] = upEx5ss
                                events[new_key]["Downstream3ss"] = downEx3ss
                                events[new_key]["Coverage"] = [ sg.areDirectlyConnected(start,end) \
                                                                for start,end in candidates[downEx3ss] ]

        return events

    def classify_MSEobsolete(self, sg):
        """identify multiple skipped exons from splice graph, defined as:
        two or more non-overapping exons which are skipped together in at least on transcript

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            # look for an exon end (common upstream 5ss: upEx5ss)
            for upEx5ss in sg.allElements():

                if sg.is5ss(upEx5ss):
                    # are there several alternative 3ss?
                    start3ss = sg.downstreamConnectedElements(upEx5ss)
                    if len(start3ss) > 1:

                        self.__inform("upEx5ss: %s \n" % upEx5ss,5)

                        # collect candidate MSE (common downstream 3ss: downEx3ss, non-overlapping)
                        candidates = {} # format: candidates[downEx3ss] = [ (start1, end1), (start2, end2), ... ]
                        # counts the exons travelled through
                        exonsTravelled = 0
                       
                        for start in start3ss:
                            exStarts = [start] # a stack of 3ss encountered in downstream exons
                            exonHistory = {start: 0} # tracks each exons history
                            exonPath = {}
                        
                            while len(exStarts) > 0:
                                self.__inform(" ".join(s for s in exStarts)+"\n",5)
                                
                                newStart = exStarts.pop(0)
                                ends = sg.downstreamConnectedElements(newStart)
                                
                                for end in ends:
                                    for downEx3ss in sg.downstreamConnectedElements(end):
                                        if not downEx3ss in start3ss:
                                            self.__inform("\t\tadding %s\n"%downEx3ss,5)
                                            exStarts.append(downEx3ss) # add to stack
                                        
                                        # needs to collect whole path not only the nb of exons
                                        # also at present it may erase some of the other history of 2 paths converge again
                                        if exonPath.has_key(newStart):
                                            exonPath[downEx3ss] = exonPath[newStart] + [(newStart,end)]
                                        else:
                                            exonPath[downEx3ss] =  [(newStart,end)]
                                            
                                        exonHistory[downEx3ss] = exonHistory[newStart] + 1 
                                        if downEx3ss in start3ss and exonHistory[downEx3ss] > 1:
                                            if candidates.has_key(downEx3ss):
                                                # candidates[downEx3ss].append(exonPath[downEx3ss]) # (newStart,end, downEx3ss, exonHistory[downEx3ss]))
                                                pass
                                            else:
                                                candidates[downEx3ss] = exonPath[downEx3ss] #[(newStart,end, downEx3ss, exonHistory[downEx3ss])]
                                                
                                self.__inform("\t"+" ".join(k for k in exonHistory.keys()) +"\n",5)
                                    
                            #if len(candidates.keys()) > 0:
                            #    for key in candidates.keys():
                            #        print key, candidates[key]

                            for downEx3ss in candidates.iterkeys():
                                new_key = "%s;%s" % (upEx5ss,downEx3ss)
                                #print new_key, candidates[downEx3ss]
                                if not events.has_key(new_key):
                                    events[new_key] = {}
                                    events[new_key]["path"] = candidates[downEx3ss]
                                    events[new_key]["skCoverage"] = sg.areDirectlyConnected(upEx5ss,downEx3ss)
                                    events[new_key]["nbSkExons"] = len(candidates[downEx3ss])
                                
        return events

    def classify_MSE(self, sg):
        """identify multiple skipped exons from splice graph, defined as:
        two or more non-overapping exons which are skipped together in at least on transcript

        --> dict"""

        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            mult5ss = [] # list of 5ss with multiple downstream 3ss splices
            mult3ss = [] # list of 3ss with multiple upstream 5ss splices
            

            # look for an exon end (common upstream 5ss: upEx5ss)
            for el in sg.allElements():
                if sg.is5ss(el) and sg.downstreamConnectedElements(el) > 1:
                    mult5ss.append(el)
                elif sg.is3ss(el) and sg.upstreamConnectedElements(el) > 1:
                    mult3ss.append(el)

            for m5ss in mult5ss:
            	#print '5ss: %s' % m5ss
                for m3ss in mult3ss:
                    #print '5ss: %s 3ss: %s' % (m5ss,m3ss)

                    if sg.areDirectlyConnected(m5ss, m3ss):
                        # potential mse

                        for up5ss in sg.upstreamConnectedElements(m3ss):
                            for dn3ss in sg.downstreamConnectedElements(m5ss):

                                if sg.areConnected(dn3ss, up5ss) and not sg.areDirectlyConnected(dn3ss, up5ss):
                                    # currenly only stores the ss between which the multiple skipping occurs
                                    # should also go after the paths, perhaps in a SpliceGraph function: get_paths(el1, el2)
                                    events["%s,%s" % (m5ss,m3ss)] = 1
        return events

    def classify_AFE_experimental(self, sg, maxNrElementsBeforeConverge=3):
        """Identify alternative first exons, defined as two or more TSS converging on a 3'ss
        via a nubmer of elements ( nrElements > 1 <= maxNrElementsBeforeConverge).

        --> dict"""
        
        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            transcript_start_sites = sg.allTypes()["TSS"].keys()
            if len(transcript_start_sites) > 1:

                # multiple TSS in splice graph
                ## print transcript_start_sites

                for tss1 in transcript_start_sites:
                    for tss2 in transcript_start_sites:
                        if tss1 == tss2: continue

                        tss1_downstream_elements = []
                        tss2_downstream_elements = []
                        tss1_recent_downstream_elements = [tss1]
                        tss2_recent_downstream_elements = [tss2]
                        pathsConverged = False

                        for nrElementsInPath in xrange(maxNrElementsBeforeConverge):
                            for testEl in tss1_recent_downstream_elements:
                                for downEl in sg.downstreamConnectedElements(testEl):
                                    #if nrElementsInPath > 1: # not saving the 5'ss of the first exons, since independent AFE is sought for
                                    tss1_downstream_elements.append(downEl)
                                    tss1_recent_downstream_elements.append(downEl)
                                del tss1_recent_downstream_elements[tss1_recent_downstream_elements.index(testEl)]
                            for testEl in tss2_recent_downstream_elements:
                                for downEl in sg.downstreamConnectedElements(testEl):
                                    #if nrElementsInPath > 1:
                                    tss2_downstream_elements.append(downEl)
                                    tss2_recent_downstream_elements.append(downEl)
                                del tss2_recent_downstream_elements[tss2_recent_downstream_elements.index(testEl)]

                            for tss1_downEl in tss1_downstream_elements:
                                for tss2_downEl in tss2_downstream_elements:
                                    # converge and  not directly downstream of TSS 1 and not directly downstream of TSS 2
                                    if tss1_downEl == tss2_downEl and \
                                       not tss1_downEl in sg.downstreamConnectedElements(tss1) and \
                                       not tss2_downEl in sg.downstreamConnectedElements(tss2): 
                                        # paths have converged
                                        pathsConverged = True
                                        if nrElementsInPath > 1:
                                            if not events.has_key("%s,%s" % (sg.downstreamConnectedElements(tss2),
                                                                             sg.downstreamConnectedElements(tss1))):
                                                
                                                events["%s,%s" % (sg.downstreamConnectedElements(tss1),
                                                                  sg.downstreamConnectedElements(tss2))] = \
                                                                  [tss1_downstream_elements,tss2_downstream_elements]
                                        break
                                    
                                if pathsConverged is True:
                                    break
                            if pathsConverged is True:
                                    break
                                
        return events

    def classify_AFE(self, sg):
        """Identify all alternative first exons, defined as the longest exon (in bp TSS to downstream 5'ss)
        --> dict"""
        
        events = {}

        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):

            transcript_start_sites_elements = sg.allTypes()["TSS"].keys()
            if len(transcript_start_sites_elements) > 1:

                # multiple TSS in splice graph
                ## print transcript_start_sites
                alternative_first_exons = {}
                alternative_first_exons_ends = {}
                longest_alternative_first_exons = []

                for tss in transcript_start_sites_elements:
                    downstreamEl = sg.downstreamConnectedElements(tss)[0]
                    if not alternative_first_exons_ends.has_key(downstreamEl):
                        alternative_first_exons_ends[downstreamEl] = {"TSS":[],
                                                                      "longest":""}
                    alternative_first_exons_ends[downstreamEl]["TSS"].append(tss)
                    alternative_first_exons[tss] = downstreamEl

                if len(alternative_first_exons_ends.keys()) == 1:
                    # there is only one region with many TSS mapping to the same TER or 5'ss
                    # not considered to be a AFE event
                    return events

                # only save the longest AFE if multiple TSS map to the same 5'ss
                for end in alternative_first_exons_ends.keys():
                    (chr, coord, strand) = end.split(":")
                    longest = 0
                    for upTSS in alternative_first_exons_ends[end]["TSS"]:
                        if strand == "-":
                            l = int(upTSS.split(":")[1]) - int(coord)
                            if l > longest:
                                alternative_first_exons_ends[end]["longest"] = upTSS
                                longest = l
                        else:
                            l =  int(coord) - int(upTSS.split(":")[1])
                            if l > longest:
                                alternative_first_exons_ends[end]["longest"] = upTSS
                                longest = l
                                
                    # save the first exons
                    longest_alternative_first_exons.append((alternative_first_exons_ends[end]["longest"], end,len(alternative_first_exons_ends[end]["TSS"])))
                    
                events[",".join( "%s::%s::%s" % (e[0],e[1],e[2]) for e in longest_alternative_first_exons)] = 1
            
        return events
            
    

    def prettyString(self, events=None, type=None, gnId=None):
        "return list of strings representing the alternative events"

        list = []

        if events is None and type is not None:
            #return a header string
            if type == 'CE' or type == 'CI':
                list.append("#gnId\tstartElement\tendElement\tcoverage\n")
            # !!!
            elif type == 'CEP' or type == 'CIP':
            	list.append("#gnId\tstartElement\tendElement\tcoverage\tupstreamConnected\tdownstreamConnected\n")
            elif type == 'SE':
                list.append("#gnId\tstartElement\tendElement\tskippedCoverage\tincludedCoverage\n")
            elif type == 'SEP':
                list.append("#gnId\tstartElement\tendElement\tskippedCoverage\tincludedCoverage\tupconnected\tdownconnected\n")
            elif type == 'RI':
                list.append("#gnId\tstartElement\tendElement\tretainedCoverage\tsplicedCoverage\n")
            elif type == 'RIP':
                list.append("#gnId\tstartElement\tendElement\tretainedCoverage\tsplicedCoverage\tupconnected\tdownconnected\n")
            elif type == 'A5SSobsolete':
                list.append("#gnId\tfirst5SS\tsecond5SS\tfirstCoverage\tsecondCoverage\tupstream3SS\tdownstream3SS\n")
            elif type == 'A5SS':
                list.append("#gnId\tnbAlt5SS\talt5SS\tcoverage\tupstream3SS\tdownstream3SS\n")
            elif type == 'A3SSobsolete':
                list.append("#gnId\tfirst3SS\tsecond3SS\tfirstCoverage\tsecondCoverage\tupstream5SS\tdownstream5SS\n")
            elif type == 'A3SS':
                list.append("#gnId\tnbAlt3SS\talt3SS\tcoverage\tupstream5SS\tdownstream5SS\n")
            elif type == 'A3SSP':
                list.append("#gnId\tnbAlt3SS\talt3SS\tcoverage\tupstream5SS\tdownstream5SS\tupExons\n")
            elif type == 'MXEobsolete':
                list.append("#gnId\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd\tfirstCoverage\tsecondCoverage\tupstream5SS\tdownstream3SS\tdirectSplicingCoverage\n")
            elif type == 'MXE':
                list.append("#gnId\tnbMXE\tMXE\tcoverage\tupstreamExon5SS\tdownstreamExon3SS\n")
            elif type == 'MSE':
                list.append("#gnId\tMSE5ss\tMSE3ss\n")
            elif type == 'AFE':
                list.append("#gnId\tTSS::downstream5ss::nrTSS\n")

        elif isinstance(events, dict) and type is not None:
        	
            for event in events.iterkeys():
                if type == 'CE' or type == 'CI':
                    for partner in events[event].keys():
                        list.append("%s\t%s\t%s\t%i\n" % (gnId, event, partner, events[event][partner]))
                        
                elif type == 'CEP' or type == 'CIP':
                    for partner in events[event].keys():
                        list.append("%s\t%s\t%s\t%i\t%s\t%s\n" % (gnId, event, partner, events[event][partner][0], events[event][partner][1], events[event][partner][2]))

                elif type == 'SE' or type == 'RI':
                    for partner in events[event].keys():
                        list.append("%s\t%s\t%s\t%i\t%i\n" % (gnId, event, partner, events[event][partner][0], events[event][partner][1]))
                
                elif type == 'SEP' or type == 'RIP':
                    for partner in events[event].keys():
                        list.append("%s\t%s\t%s\t%i\t%i\t%s\t%s\n" % (gnId, event, partner, events[event][partner][0], events[event][partner][1], events[event][partner][2], events[event][partner][3] ))

                elif type == 'A5SSobsolete':
                    A5SS = event.split(',')
                    list.append("%s\t%s\t%s\t%i\t%i\t%s\t%s\n" % (gnId, A5SS[0], A5SS[1], \
                                                              events[event]["Coverage"][0], events[event]["Coverage"][1], \
                                                              events[event]["Upstream3ss"], events[event]["Downstream3ss"]))
                elif type == 'A5SS':
                    list.append("%s\t%i\t%s\t%s\t%s\t%s\n" % (gnId, len(event.split(',')), event, \
                                                              ','.join(["%s" % c for c in events[event]["Coverage"]]), \
                                                              events[event]["Upstream3ss"], events[event]["Downstream3ss"]))
                elif type == 'A3SSobsolete':
                    A3SS = event.split(',')
                    list.append("%s\t%s\t%s\t%i\t%i\t%s\t%s\n" % (gnId, A3SS[0], A3SS[1], \
                                                              events[event]["Coverage"][0], events[event]["Coverage"][1], \
                                                              events[event]["Upstream5ss"], events[event]["Downstream5ss"]))
                elif type == 'A3SS':
                    list.append("%s\t%i\t%s\t%s\t%s\t%s\n" % (gnId, len(event.split(',')), event, \
                                                              ','.join(["%s" % c for c in events[event]["Coverage"]]), \
                                                              events[event]["Upstream5ss"], events[event]["Downstream5ss"]))
                elif type == 'A3SSP':
                    list.append("%s\t%i\t%s\t%s\t%s\t%s\t%s\n" % (gnId, len(event.split(',')), event, \
                                                              ','.join(["%s" % c for c in events[event]["Coverage"]]), \
                                                              events[event]["Upstream5ss"], events[event]["Downstream5ss"], events[event]["upExons"]))
                elif type == 'MXEobsolete':
                    list.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gnId, events[event]["exon1"][0],events[event]["exon1"][1],
                                                                             events[event]["exon2"][0],events[event]["exon2"][1],
                                                                             events[event]["exon1"][2],events[event]["exon2"][2],
                                                                             event,events[event]["downstream3ss"],events[event]["bothSkipped"]))
                elif type == 'MXE':
                    list.append("%s\t%i\t%s\t%s\t%s\t%s\n" % (gnId, len(event.split(',')), event, \
                                                              ','.join(["%s" % c for c in events[event]["Coverage"]]), \
                                                              events[event]["Upstream5ss"],events[event]["Downstream3ss"]))
                elif type == "MSE":
                    list.append("%s\t%s\t%s\n" % (gnId, event.split(',')[0], event.split(',')[1]))

                elif type == "AFE":
                    list.append("%s\t%s\t\n" % (gnId, event))
        return list


    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            print >> sys.stderr, txt,


def testConstitutiveExons():
    print "testing classify_CE"
    import Classifier

    c = Classifier.Classifier(verbose=0)
    for gnId in ('chr21:46879507:+','chr21:26191805:+'):
        testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % gnId
        events = c.classify_CE(testfilename)
        print "\tgene:", gnId
        for exStart in events.keys():
            for exEnd in events[exStart].keys():
                print "\t\t%s -> %s: constitutive : coverage: %i" % (exStart, exEnd, events[exStart][exEnd])


def testSkippedExons():
    print "testing classify_SE"
    import Classifier
    c = Classifier.Classifier(verbose=0)
    for gnId in ('chr21:22041176:+', 'chr21:10120851:-'):
        testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % gnId
        events = c.classify_SE(testfilename)
        print "\tgene:", gnId
        for exStart in events.keys():
            for exEnd in events[exStart].keys():
                print "\t\t%s -> %s: skipped: %i, included: %i" % (exStart, exEnd, events[exStart][exEnd][0], events[exStart][exEnd][1])

def test5ASS():
    print "testing classify_A5SS"
    import Classifier
    c = Classifier.Classifier(verbose=0)
    for gnId in ('chr21:22041176:+', 'chr21:10120851:-','chr21:46879507:+','chr21:26191805:+'):
        testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % gnId
        print "\tgene:", gnId
        events = c.classify_A5SS(testfilename)
        #print "Testing Classifier of alternative 5' splice sites"
        for key in events.keys():
            A5SS = key.split(',')
            print '\t\t%s , %s : alternative 5ss : coverage: %i, %i' % (A5SS[0],A5SS[1],events[key]["Coverage"][0],
                                                                          events[key]["Coverage"][1])
            print '\t\tupstream 3ss: %s, downstream 3ss:%s' % (events[key]["Upstream3ss"],events[key]["Downstream3ss"]) 

            print
    print "testing classify_A3SS"
    for gnId in ('chr21:22041176:+', 'chr21:10120851:-','chr21:46879507:+','chr21:26191805:+'):
        testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % gnId
        print "\tgene:", gnId
        events = c.classify_A3SS(testfilename)
        #print "Testing Classifier of alternative 3' splice sites"
        for key in events.keys():
            A3SS = key.split(',')
            print '\t\t%s , %s : alternative 3ss : coverage: %i, %i' % (A3SS[0],A3SS[1],events[key]["Coverage"][0],
                                                                          events[key]["Coverage"][1])
            print '\t\tupstream 5ss: %s, downstream 5ss:%s' % (events[key]["Upstream5ss"],events[key]["Downstream5ss"]) 

            print


def testAll():
    import Classifier
    c = Classifier.Classifier(verbose=0)
    for gnId in ('chr21:22041176:+', 'chr21:10120851:-'):
        testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % gnId

        print "\tgene:", gnId


        events = c.classify_MXE(testfilename)
        for exStart in events.keys():
            for entries in events[exStart].keys():
                print exStart, entries, events[exStart][entries]
        #print events
                
        # Alt 5'ss
##         events = c.classify_A5SS(testfilename)
##         print events
##         print "Testing Classifier of alternative 5' splice sites"
##         for exStart in events.keys():
##             for exEnd in events[exStart].keys():
##                 print '\t\t%s -> %s : alternative 5ss : Transcripts through exon: %i, intron: %i' % (exStart,
##                                                                                                      exEnd,
##                                                                                                      events[exStart][exEnd][0],
##                                                                                                      events[exStart][exEnd][1])
##             print

        # SE
        events = c.classify_SE(testfilename)
        for exStart in events.keys():
            for exEnd in events[exStart].keys():
                print "\t\t%s -> %s: skipped: %i, included: %i" % (exStart, exEnd, events[exStart][exEnd][0], events[exStart][exEnd][1])
        
           
            
if __name__ == "__main__":
    import Classifier
    c = Classifier.Classifier(verbose=1)
    print >> sys.stdout, Classifier.__doc__

    #testConstitutiveExons()
    #testSkippedExons()
    #test5ASS()
    #testAll()

    #for gnId in ['chr21:46879507:+']:
    #for gnId in ('chr21:22041176:+', 'chr21:10120851:-','chr21:46879507:+','chr21:26191805:+'):
    iterSG = SpliceGraph.SpliceGraphIterator('/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr21')
    sg = iterSG.next_sg()
    while sg:
        #testfilename = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr21/%s.sg' % gnId
        #print "\tgene:", gnId
        #c.build_connection_matrix(testfilename)
        #ev = c.classify_MSE(testfilename)
        #ev = c.classify_AFE(sg)
        ev = c.classify_all_longest_AFE(sg)
        #ev = c.classify_MSE(sg)
        print sg.name
        for key in ev.iterkeys():
            print "\t", key, ev[key]
        sg = iterSG.next_sg()
