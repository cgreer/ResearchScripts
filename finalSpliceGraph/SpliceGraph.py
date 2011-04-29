#---------------------------------------------------------------------------------------
#   $Id: SpliceGraph.py,v 1.78 2006/03/13 15:58:42 stadler Exp $    
#---------------------------------------------------------------------------------------
"""SpliceGraph: Class for representing a single splice graph

Synopsis:
    import SpliceGraph

    # three exon gene, exon 2 is skipped and included
    sg = SpliceGraph.SpliceGraph(connections={'chrX:100:+':{'chrX:200:+':2},
                                              'chrX:200:+':{'chrX:1000:+':-1,'chrX:2100:+':-1},
                                              'chrX:1000:+':{'chrX:1100:+':1},
                                              'chrX:1100:+':{'chrX:2100:+':-1},
                                              'chrX:2100:+':{'chrX:2200:+':2}},
                                 types={'TSS':{'chrX:100:+':2},
                                        '3ss':{'chrX:1000:+':1,'chrX:2100:+':2},
                                        '5ss':{'chrX:200:+':2,'chrX:1100:+':1},
                                        'TER':{'chrX:2200:+':2}},
                                 transcripts={'chrX:100:+':{'chrX:200:+':['tx1','tx2']},
                                              'chrX:200:+':{'chrX:1000:+':['tx1'],'chrX:2100:+':['tx2']},
                                              'chrX:1000:+':{'chrX:1100:+':['tx1']},
                                              'chrX:1100:+':{'chrX:2100:+':['tx1']},
                                              'chrX:2100:+':{'chrX:2200:+':['tx1','tx2']}},
                                 species='hopschel',
                                 assembly='fischpa1',
                                 )

    sg.store('file.sg')
    sg.visualize('file.ps')

    sg.read('file.sg')

    sg2 = SpliceGraph.SpliceGraph(file='file.sg')

    list  = sg.allElements()

    if sg.is3ss(list[0]) and sg.is5ss(list[1]) and sg.areDirectlyConnected(list[0], list[1]):
        print list[0], '->', list[1], 'is an exons'

    list1 = sg.upstreamElements(list[3])
    list2 = sg.downstreamElements(list[3])
"""

# standard lib
import os
import sys
import copy
import re

# project related modules
import Errors
#import SpliceSites
import configuration
import GenomeFetch

# 3rd party modules
#import pydot
import MySQLdb, _mysql


class SpliceGraph:
    "SpliceGraph: Object that represents a single SpliceGraph"

    INF = sys.maxint - 1

    # how should elements be drawn?
    __defaultShapes = {'TSS'   : 'invtriangle',
                       'TER'   : 'triangle',
                       '5ss'   : 'circle',
                       '3ss'   : 'box',
                       'multi' : 'doubleoctagon'}
    __defaultColors = {'TSS'   : 'lawngreen',
                       'TER'   : 'red',
                       '5ss'   : 'steelblue1',
                       '3ss'   : 'steelblue1',
                       'multi' : 'orange',
                       'emph'  : 'orange',
                       'intron': 'gray',
                       'exon'  : 'navy'}
    __defaultFontsize = 10


    __defaultTypes = ['5ss','3ss','TER','TSS']
    

    # standard genetic code
    __table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
        'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
        'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
        'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
        'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
        'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
        'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
        'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
        'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
        'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
        'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*', }
    __start_codons = {'TTG':1, 'CTG':1, 'ATG':1, }
    __stop_codons  = {'TAA':1, 'TAG':1, 'TGA':1, }

    def __init__(self, name=None, connections=None, types=None, transcripts=None, doChecks=True,
                 filename=None, verbose=False, species="unknown", assembly="unknown", conf=None, log=sys.stderr):
        "instantiate SpliceGraph object based on element connection/types --> SpliceGraph"

        # class attributes
        self.verbose       = verbose
        self.log           = log           # file handle for self.__inform()
        self.name          = name          # equal to name of most upstream element
        self.__connections = connections   # dict[el1][el2] = nb_times
        self.__types       = types         # dict[TYPE][el] = nb_times
        self.__transcripts = transcripts   # dict[el1][el2] = [tx1, tx2, ...]
        self.__elements    = []            # sorted list of all elements (in biological order!)
        self.__strand      = None          # taken from element names
        self.__distances   = {}            # caches the distances calculated by areConnected() / connectedComponents()
        self.__tid2lid     = {}            # caches library id by transcript id
        self.__libkeys     = {}            # caches keywords for libraries by LID
        self.__keysyns     = {}            # caches keyword synonyms by keyword
        self.__path        = {}            # caches element indices corresponding to transcript id (elInxByTranscript)
        self.__tmpdir      = "/tmp"        # used to do temporary pairwise alignments
        self.species       = species
        self.assembly      = assembly
        self.__db          = None          # connection to MySQL, used to retreive library information
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = None #get it later, if needed

        if doChecks:
            #check arguments
            if filename is not None: #init from file
                self.read(filename, name)

            elif ( connections is None or not isinstance(connections, dict) or \
                   types is None or not isinstance(types, dict)             or \
                   transcripts is None or not isinstance(transcripts, dict) ): #init from dictionaries
                self.__inform("ERROR: Missing data/types to generate object from variables in SpliceGraph.__init__()\n")
                raise Errors.ObjectInitError('SpliceGraph.__init__', 'Missing data/types to generate object from variables')

            #make sure that the __defaultTypes keys exist in self.__types
            for type in self.__defaultTypes:
                if not self.__types.has_key(type):
                    self.__types[type] = {}

            #build ordered list of all elements (also sets self.__strand)
            self.__buildElementsList()

            #check if types are available
            for el in self.__elements:
                elementFound = False
                for type in self.__types.iterkeys():
                    if self.__types[type].has_key(el):
                        elementFound = True
                        break
                if not elementFound:
                    self.__inform("ERROR: Missing type for element %s in SpliceGraph.__init__()\n" % el)
                    raise Errors.ObjectInitError('SpliceGraph.__init__', "Missing type for element '%s'" % el)
                    #sys.exit(0)

            #maybe should do some more consistency checking here

        else:
            self.__inform("WARNING: so you don't want checking? then it's all up to you!\n", 1)

        if self.name is None:
            self.name = "unknown"


    def __connect(self):
        "Establish the connection with MySQL"

        if self.conf is None:
            self.conf = configuration.Configuration()

        self.__db = MySQLdb.connect(db     = self.conf[self.species]["Assembly"],
                                    user   = os.getlogin(),
                                    host   = self.conf[self.species]["MySQLHost"],
                                    passwd = '')
        
    def genomicLength(self):
        "return the total length covered (from the most upstream to the most downstream element) --> int"
        start, end = self.genomicRange()
        return end - start + 1


    def genomicRange(self):
        "return the (start,end)-tuple corresonding to extreme genomic coordinates covered by splice graph (start<=end) --> tuple"
        if len(self.__elements) > 0:
            start = int(self.__elements[0].split(':')[1])
            end   = int(self.__elements[-1].split(':')[1])
            if self.__strand == '-':
                start, end = end, start
            return (start, end)
        else:
            return (1,0)


    def fastCopy(self):
        "build and return a stripped-down copy of the object (only containing connections/elements)"

        fastcopy = SpliceGraph(connections={}, doChecks=False, verbose=0)
        for k in self.__connections.iterkeys():
            fastcopy.__connections[k] = {}
            for(l, v) in self.__connections[k].iteritems():
                fastcopy.__connections[k][l] = v
        fastcopy.__elements    = self.__elements[:]

        return fastcopy


    def __buildElementsList(self):
        "build ordered list of elements from connections/types"

        elements = {}

        for el1 in self.__connections.iterkeys():
            if not elements.has_key(el1):
                elements[el1] = 1
            for el2 in self.__connections[el1].iterkeys():
                if not elements.has_key(el2):
                    elements[el2] = 1

        # check if all elements are on same strand
        for el in elements.iterkeys():
            if self.__strand is None:
                self.__strand = el.split(':')[2]
            else:
                if self.__strand != el.split(':')[2]:
                    self.__inform("ERROR: not all elements are on same strand\n")
                    raise Errors.ObjectInitError('SpliceGraph.__buildElementsList', "Not all elements are on same strand")
        
        self.__elements = elements.keys()

        # sort list of all elements (sorting by coord, assuming "chr:coord:strand")
        # reverse the order if the gene is on the minus strand
        self.__elements.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])), \
                             reverse=(self.__strand == '-') )


    def getLibraryID(self, tid):
        "return the library id corresponding to a given transcript --> int or None"

        if self.__tid2lid.has_key(tid):
            return self.__tid2lid[tid]

        else:
            if not isinstance(self.__db, _mysql.connection):
                self.__connect()

            cursor = self.__db.cursor()
            cursor.execute("""SELECT library from gbCdnaInfo WHERE acc = "%s";""" % tid)
            res = cursor.fetchall()
            cursor.close()

            if len(res) == 1:
                self.__tid2lid[tid] = int(res[0][0])
                return int(res[0][0])
            else:
                return None


    def getLibraryKeywords(self, lid):
        "return a dict whose keys are keywords associated to a library in a CGAP annotation file --> dict"

        if self.__libkeys.has_key(lid):
            return self.__libkeys[lid]
        else:
            if not isinstance(self.__db, _mysql.connection):
                self.__connect()

            cursor = self.__db.cursor()
            cursor.execute("""SELECT KEYWORDS FROM LibData WHERE UNIGENE_LIB_ID=%i;""" % lid)
            res = cursor.fetchall()
            cursor.close()

            if len(res) == 1:
                self.__libkeys[lid] = res[0][0].split(', ')
                return self.__libkeys[lid]
            else:
                return {}


    def nbExons(self, internalOnly=False):
        "return the number of exons present in the current SpliceGraph --> int"

        nb = 0
        for exStart in self.allElements():
            if self.is3ss(exStart) or (self.isTSS(exStart) and internalOnly is False):
                for exEnd in self.downstreamConnectedElements(exStart):
                    if self.is5ss(exEnd) or (self.isTER(exEnd) and internalOnly is False):
                        nb += 1
        return nb

    def allExons(self, internalOnly=False):
        "return a list of all exons present in the current SpliceGraph, each entry being (exStart,exEnd) --> list of tuples"
        exons = []
        for exStart in self.allElements():
            if self.is3ss(exStart) or (self.isTSS(exStart) and internalOnly is False):
                for exEnd in self.downstreamConnectedElements(exStart):
                    if self.is5ss(exEnd) or (self.isTER(exEnd) and internalOnly is False):
                        exons.append((exStart, exEnd))

        return exons


    def allExonSeqs(self, internalOnly=False, calcPhases=True):
        "return a list of all exon sequnce objects present in the current SpliceGraph --> list of Exon objects"

        # get phase information
        phases = {}
        if calcPhases:
            phases = self.getExonPhases()

        # generate the Exon objects
        exonObjs = []
        for ex in self.allExons(internalOnly=internalOnly):
            startPhase = -1
            endPhase   = -1
            if phases.has_key(ex[0]) and phases[ex[0]].has_key(ex[1]):
                # if there are several sources of phases available for this exon, us the one from the longest CDS
                maxCdsLen = 0
                for phaseinfo in phases[ex[0]][ex[1]]:
                    if phaseinfo[1] > maxCdsLen:
                        txId, maxCdsLen, startPhase, endPhase = phaseinfo
            exonObjs.append(Exon(ex[0], ex[1], startPhase=startPhase, endPhase=endPhase))

        return exonObjs


    def getExonPhases(self):
        "calculate phases of exons in current splice graph using CDS annotations of known genes --> dict"
        
        chr = self.name.split(':')[0]
        phases = {} # format: startPhases[start][end] = [(txId, cdsLen, startPhase, endPhase), ...]

        if not isinstance(self.__db, _mysql.connection):
            self.__connect()

        cursor = self.__db.cursor()
        cursor.execute("""SELECT name, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds FROM %s WHERE chrom='%s' AND FIND_IN_SET(name,'%s')>0;""" \
                       % (self.conf[self.species]["knownGene"], chr, ','.join(self.allTranscripts()) ) )
        res = cursor.fetchall()
        cursor.close()

        for row in res:
            txId     = row[0]
            cdsStart = int(row[1])+1 # convert to 1-based coord-system
            cdsEnd   = int(row[2])
            nbBlocks = int(row[3])
            starts   = [int(n)+1 for n in row[4].tostring().split(',') if len(n)] # convert to 1-based coord-system
            ends     = [int(n)   for n in row[5].tostring().split(',') if len(n)]
            cumCDS   = 0 # cumulative length of CDS
            # make sure blocks are in biological order
            if self.__strand == '-':
                starts.reverse()
                ends.reverse()
            for b in xrange(nbBlocks):
                # get biological start and end element names ("chr:coord:strand")
                el1 = None
                el2 = None
                if self.__strand == '+':
                    el1 = ':'.join([chr, str(starts[b]), '+'])
                    el2 = ':'.join([chr, str(ends[b]),   '+'])
                else:
                    el1 = ':'.join([chr, str(ends[b]),   '-'])
                    el2 = ':'.join([chr, str(starts[b]), '-'])
                # calculate CDS length in this block
                cdsLen = 0
                if starts[b] >= cdsStart and ends[b] <= cdsEnd:                            # fully coding
                    cdsLen = (ends[b] - starts[b] + 1)
                elif starts[b] < cdsStart and ends[b] > cdsEnd:                            # partially coding, single exon CDS
                    cdsLen = (cdsEnd - cdsStart + 1)
                elif starts[b] >= cdsStart and starts[b] <= cdsEnd and ends[b] > cdsEnd:   # partially coding, start CDS
                    cdsLen = (cdsEnd - starts[b] + 1)
                elif starts[b] < cdsStart and ends[b] >= cdsStart and ends[b] <= cdsEnd:   # partially coding, end CDS
                    cdsLen = (ends[b] - cdsStart + 1)
                # calculate phases for this block
                startPhase = -1
                endPhase   = -1
                if self.__strand == '+':
                    if starts[b] >= cdsStart and starts[b] <= cdsEnd:
                        startPhase = cumCDS % 3
                    if ends[b] >= cdsStart and ends[b] <= cdsEnd:
                        endPhase = (cumCDS + cdsLen) % 3 
                else:
                    if ends[b] <= cdsEnd and ends[b] >= cdsStart:
                        startPhase = cumCDS % 3
                    if starts[b] <= cdsEnd and starts[b] >= cdsStart:
                        endPhase = (cumCDS + cdsLen) % 3
                # update cumCDS
                cumCDS += cdsLen
                # check for ambiguous phase information and store phases
                if phases.has_key(el1):
                    if phases[el1].has_key(el2):
                        self.__inform("WARNING: Multiple sources of phase information (%i,%i) for exon %s:E:%s\n" % (startPhase, endPhase, el1, el2) )
                        phases[el1][el2].append((txId, (cdsEnd - cdsStart + 1), startPhase, endPhase))
                    else:
                        phases[el1][el2] = [(txId, (cdsEnd - cdsStart + 1), startPhase, endPhase)]
                else:
                    phases[el1] = {}
                    phases[el1][el2] = [(txId, (cdsEnd - cdsStart + 1), startPhase, endPhase)]

        return phases


    def nbIntrons(self):
        "return the number of introns present in the current SpliceGraph --> int"

        nb = 0
        for inStart in self.allElements():
            if self.is5ss(inStart):
                for inEnd in self.downstreamConnectedElements(inStart):
                    if self.is3ss(inEnd):
                        nb += 1
        return nb

    def allIntrons(self):
        "return a list of all introns present in the current SpliceGraph, each entry being (inStart,inEnd)"

        introns = []
        for inStart in self.allElements():
            if self.is5ss(inStart):
                for inEnd in self.downstreamConnectedElements(inStart):
                    if self.is3ss(inEnd):
                        introns.append((inStart,inEnd))
        return introns


    def allElements(self):
        "return a sorted list of all elements, lower coordinates first"
        return self.__elements


    def allConnections(self):
        "return a dictionary of dictionaries defining element connections"
        return self.__connections


    def allTypes(self):
        "return a dictionary of dictionaries defining element types"
        return self.__types


    def allTranscripts(self):
        "return a list of transcript ids that provied evidence for this splice graph"
        txs = {}
        for el1 in self.__transcripts.keys():
            for el2 in self.__transcripts[el1].keys():
                for tx in self.__transcripts[el1][el2]:
                    if not txs.has_key(tx):
                        txs[tx] = 1
        return txs.keys()


    def strand(self):
        "return the strand of the elements in the splice graph"
        return self.__strand


    def __removeConnection(self, el1, el2):
        "removes the direct connection between two elements"
        if self.__connections.has_key(el1):
            if self.__connections[el1].has_key(el2):
                del self.__connections[el1][el2]

    def removeTranscript(self, id):
        "complete removes all evidence that is based on a given transcript from the SpliceGraph"

        toBeRemoved = {'transcripts':[], 'types':[], 'elements':[]}

        # walk through transcript evidence, compile remove list
        for el1 in self.__transcripts.keys():

            for el2 in self.__transcripts[el1].keys():

                if id in self.__transcripts[el1][el2]:               # one el1->el2 connection has to go

                    toBeRemoved['transcripts'].append((el1, el2))    #    from .__transcripts / .__connections

                    # PROBLEM:
                    # for elements with more than one type,
                    # there is no way knowing which one
                    # to remove
                    # current implementation: remove first
                    # reported by self.getType
                    for el in [el1, el2]:
                        type = self.getType(el)
                        toBeRemoved['types'].append((type,el))       #    form .__types

        # now remove all that is based on transcript 'id'
        for (el1, el2) in toBeRemoved['transcripts']:                #    from .__transcripts
            if self.__transcripts.has_key(el1) and self.__transcripts[el1].has_key(el2):
                self.__transcripts[el1][el2].remove(id)
                if len(self.__transcripts[el1][el2]) == 0:
                    del self.__transcripts[el1][el2]
            if self.__connections.has_key(el1) and self.__connections[el1].has_key(el2):
                self.__connections[el1][el2] -= 1                    #    from .__connections
                if self.__connections[el1][el2] == 0:
                    del self.__connections[el1][el2]

        for (type, el) in toBeRemoved['types']:
            if self.__types.has_key(type) and self.__types[type].has_key(el):
                self.__types[type][el] -= 1                          #    form .__types
                if self.__types[type][el] == 0:
                    del self.__types[type][el]
                    toBeRemoved['elements'].append(el)

        for el in toBeRemoved['elements']:                           #    from .__elements
            self.removeElement(el)



    def removeElement(self, el):
        "completely removes the element and all associated connections from the SpliceGraph"

        if el in self.__elements:

            for upstreamEl in self.upstreamConnectedElements(el):  # upstream transcript evidence
                if self.__transcripts[upstreamEl].has_key(el):
                    del self.__transcripts[upstreamEl][el]

            if self.__transcripts.has_key(el):                     # downstream transcript evidence
                del self.__transcripts[el]

            for upstreamEl in self.upstreamConnectedElements(el):  # upstream connections
                if self.__connections[upstreamEl].has_key(el):
                    del self.__connections[upstreamEl][el]
            
            if self.__connections.has_key(el):                     # downstream connections
                del self.__connections[el]

            types = self.getAllType(el)                            # type(s)
            for type in types:
                if self.__types[type].has_key(el):
                    del self.__types[type][el]

            if el in self.__elements:                              # element itself
                el_index = self.__elements.index(el)
                del self.__elements[el_index]
        else:
            self.__inform("RemoveElement: tried to remove nonexistent element in elements %s" % el,4)

    def moveJunction(self, old5ss, new5ss, old3ss, new3ss):
        "moves a junction and all its connections, transcript evidence and types to new or existing elements"
        	
        #maybe add a statement here that if it is the chr14 one don't move the junction...
        if old5ss in self.__elements and old3ss in self.__elements:

            tx_list_in_swap = self.__transcripts[old5ss][old3ss]

            # A1. Create transcript evidence for new junction
            if not self.__transcripts.has_key(new5ss):
                self.__transcripts[new5ss] = {}
            if self.__transcripts[new5ss].has_key(new3ss):  # move tx evidence
                for tx in self.__transcripts[old5ss][old3ss]:
                    self.__transcripts[new5ss][new3ss].append(tx) # add to existing transcript list
            else:
                self.__transcripts[new5ss][new3ss] = self.__transcripts[old5ss][old3ss] # create new tx list

            # B1. Create new connection
            if not self.__connections.has_key(new5ss):
                self.__connections[new5ss] = {}               
            if self.__connections[new5ss].has_key(new3ss):
                self.__connections[new5ss][new3ss] += self.__connections[old5ss][old3ss] # add to existing connection
            else:
                self.__connections[new5ss][new3ss] = self.__connections[old5ss][old3ss] # create new connection

            # DELETE the old splice
            # Delete transcript evidence for old junction
            del self.__transcripts[old5ss][old3ss]  # delete old transcript list
            #  Delete old connection
            del self.__connections[old5ss][old3ss]  # delete old connection

            # A2. Update the upstream and downstream transcript evidence for ONLY the transcripts
            # that gave evidence for the moved junction
            for upstreamEl in self.upstreamConnectedElements(old5ss):
                txToSwap = 0 # a counter for the nb of tx in swaped junction is giving evidence to upstream edge
                for tx in self.__transcripts[upstreamEl][old5ss]:
                    if tx in tx_list_in_swap:
                        txToSwap += 1
                        
                if txToSwap > 0:
                                                 
                    if not self.__transcripts[upstreamEl].has_key(new5ss):           
                        self.__transcripts[upstreamEl][new5ss] = []
                    # move transcript evidence of upstreamEl -> 5ss
                    for tx in self.__transcripts[upstreamEl][old5ss]:
                        if tx in tx_list_in_swap:
                            self.__transcripts[upstreamEl][new5ss].append(tx)
                            index = self.__transcripts[upstreamEl][old5ss].index(tx)
                            del self.__transcripts[upstreamEl][old5ss][index] # delete a tx from list
                        
                    if self.__connections[upstreamEl].has_key(new5ss):
                        self.__connections[upstreamEl][new5ss] += txToSwap # add to existing coverage
                    else:
                        self.__connections[upstreamEl][new5ss] = txToSwap  # create new coverage

                    # remove the moved coverage, but first checking if it is a intron or exon
                    # this is needed since some nodes are both 3ss and 5ss, otherwise this will
                    # lead to a non-extinct node that will make an internal inconsistency in
                    # __connections vs. __transcripts.
                    if self.__connections[upstreamEl][old5ss] > 0:
                        self.__connections[upstreamEl][old5ss] -= txToSwap  
                    else:
                        self.__connections[upstreamEl][old5ss] += txToSwap 
                        
                    # Clean up empty connections (i.e. if all evidence was moved)
                    if self.__connections[upstreamEl][old5ss] == 0:
                        del self.__connections[upstreamEl][old5ss]
                    if len(self.__transcripts[upstreamEl][old5ss]) == 0:
                        del self.__transcripts[upstreamEl][old5ss]
                    
            if len(self.__connections[old5ss].keys()) == 0 and \
                   len(self.upstreamConnectedElements(old5ss)) == 0:
                # if the old element doesnot form another splice with a 3ss, remove from
                # and no other incoming exon edge
                self.removeElement(old5ss) # remove
            else:
                # if the element still has evidence,
                # update type-element coverage
                self.recalculateElementCoverage(old5ss,"5ss")
                            
            for downstreamEl in self.downstreamConnectedElements(old3ss):
                txToSwap = 0 # a counter for the nb of tx in swaped junction is giving evidence to downstream edge
                for tx in self.__transcripts[old3ss][downstreamEl]:
                    if tx in tx_list_in_swap:
                        txToSwap += 1
                        
                if txToSwap > 0:
                    if not self.__transcripts.has_key(new3ss):
                        self.__transcripts[new3ss] = {}
                    if not self.__transcripts[new3ss].has_key(downstreamEl):           
                        self.__transcripts[new3ss][downstreamEl] = []
                    # add to existing list
                    for tx in self.__transcripts[old3ss][downstreamEl]:
                        if tx in tx_list_in_swap:
                            self.__transcripts[new3ss][downstreamEl].append(tx)
                            index = self.__transcripts[old3ss][downstreamEl].index(tx)
                            del self.__transcripts[old3ss][downstreamEl][index] # delete a tx from list

                    if not self.__connections.has_key(new3ss):
                        self.__connections[new3ss] = {}
                    if self.__connections[new3ss].has_key(downstreamEl):
                        self.__connections[new3ss][downstreamEl] += txToSwap 
                    else:
                        self.__connections[new3ss][downstreamEl] = txToSwap
                    # Assuming the old3ss is indeed a 3ss and thus downstreamEl a 5ss
                    # However this seems to be causing problems in certain splice graphs
                    # were the same element is both 5ss and 3ss
                    if self.__connections[old3ss][downstreamEl] < 0:
                        self.__connections[old3ss][downstreamEl] += txToSwap
                    else:
                        self.__connections[old3ss][downstreamEl] -= txToSwap                   
                        ## HERE IS THE PROBLEM, IF IT IS BOTH A 5ss and 3ss, how to modify??
   
                    # remove the transcript evidence completely if no more tx support link
                    if len(self.__transcripts[old3ss][downstreamEl]) == 0: # remove tx
                        del self.__transcripts[old3ss][downstreamEl]
                    if self.__connections[old3ss][downstreamEl] == 0:  # remove connections
                        del self.__connections[old3ss][downstreamEl]
                    
            if len(self.downstreamConnectedElements(old3ss)) == 0:
                # no downstream exon edge
                if self.__connections.has_key(old3ss):
                    del self.__connections[old3ss]
                if self.__transcripts.has_key(old3ss):
                    del self.__transcripts[old3ss]

            if len(self.upstreamConnectedElements(old3ss)) == 0 and \
                   len(self.downstreamConnectedElements(old3ss)) == 0    :
                    # non connected
                    self.removeElement(old3ss) # remove
            else:
                self.recalculateElementCoverage(old3ss,"3ss")

            # E. rebuild elements list
            self.__buildElementsList()
                    
            # update type-element coverage
            self.recalculateElementCoverage(new5ss,"5ss")
            self.recalculateElementCoverage(new3ss,"3ss")

            # check for inconsistencies
            if new5ss == "chr14:105610123:-":
            	print 'transcripts: %s' % (self.__transcripts[new5ss].keys())
            	print 'connections: %s' % (self.__connections[new5ss].keys())
            	
            if len(self.__transcripts[new5ss].keys()) != len(self.__connections[new5ss].keys()):
                print "ERROR: Different number of connections and transcripts for new5ss=:", new5ss
            if self.__connections.has_key(old5ss) and len(self.__transcripts[old5ss].keys()) != len(self.__connections[old5ss].keys()):
                print "ERROR: Different number of connections and transcripts for old5ss=:", old5ss
                        
        else:
            self.__inform("SpliceGraph.moveJunction: tried to move nonexistent junction: %s, %s" % (old5ss, old3ss) ,4)

    def moveElement(self, oldEl, newEl):
        "moves an element and all its connections, transcript evidence and types to a new existing or nonexisting element"
        if oldEl in self.__elements:
            for upstreamEl in self.upstreamConnectedElements(oldEl):
                         
                if self.__transcripts[upstreamEl].has_key(newEl):             # move upstream tx evidence
                    # add to existing list
                    for tx in self.__transcripts[upstreamEl][oldEl]:
                        self.__transcripts[upstreamEl][newEl].append(tx)
                else:
                    # make new tx list
                    self.__transcripts[upstreamEl][newEl] = self.__transcripts[upstreamEl][oldEl]
                    
                   
                if self.__connections[upstreamEl].has_key(newEl):             # move upstream connections
                    self.__connections[upstreamEl][newEl] += self.__connections[upstreamEl][oldEl]
                else:
                    self.__connections[upstreamEl][newEl] = self.__connections[upstreamEl][oldEl]

            for downstreamEl in self.downstreamConnectedElements(oldEl):
                if not self.__transcripts.has_key(newEl):                    # move downstream tx evidence
                    self.__transcripts[newEl] = {}    

                if not self.__transcripts[newEl].has_key(downstreamEl):
                    # make new tx list
                    self.__transcripts[newEl][downstreamEl] = self.__transcripts[oldEl][downstreamEl]
                else:
                    # add to existing list               
                    for tx in self.__transcripts[oldEl][downstreamEl]:
                        self.__transcripts[newEl][downstreamEl].append(tx)

                if not self.__connections.has_key(newEl):                    # move downstream connections
                    self.__connections[newEl] = {}

                if self.__connections[newEl].has_key(downstreamEl):
                    self.__connections[newEl][downstreamEl] += self.__connections[oldEl][downstreamEl]
                else:
                    self.__connections[newEl][downstreamEl] = self.__connections[oldEl][downstreamEl]

            # change type information
            types = self.getAllType(oldEl)
            for type in types:
                if self.__types[type].has_key(newEl):
                    self.__types[type][newEl] += self.__types[type][oldEl]   # move type and element counts
                else:
                    self.__types[type][newEl] = self.__types[type][oldEl]

            # finally, remove old element
            self.removeElement(oldEl)

            # and importantly rebuild elements list
            self.__buildElementsList()
            
        else:
            self.__inform("SpliceGraph.moveElement: tried to move nonexistent element: %s" % oldEl,4)
            

    def areDirectlyConnected(self, el1, el2):
        "return nb of observations (==True) if the two elements are DIRECTLY connected in the splice graph"

        if self.__connections.has_key(el1) and self.__connections[el1].has_key(el2):
            return self.__connections[el1][el2]
        else:
            return False


    def areConnectedNotVia(self, el1, el2, vel1, vel2):
        "return True if there is ANY PATH in the splice graph connecting the two elements not using the specified edge"
        
        return self.ListsAreConnectedNotVia([el1],[el2],vel1,vel2)


    def ListsAreConnectedNotVia(self, list1, list2, vel1, vel2):
        "return True if there is ANY PATH in the splice graph connecting the ANY of the two lists of elements not using the specified edge"
        
        # construct a copy of the current splice graph
        #tmp = copy.deepcopy(self) #we do not need to copy everything
        tmp = self.fastCopy()
        
        # remove the edge (vel1, vel2)
        tmp.__removeConnection(vel1, vel2)

        for el1 in list1:
            for el2 in list2:
                # check for path between el1 and el2
                if tmp.areConnected(el1, el2):
                    return True

        return False


    def areConnected(self, el1, el2):
        "return True if there is ANY PATH in the splice graph connecting the two elements (Dijkstra's algorithm)"

        if self.areDirectlyConnected(el1, el2):
            return True #self.areDirectlyConnected(el1, el2)

        else:
            # use Dijkstra's algorithm that solves the single-source shortest path problem
            # for a directed graph with nonnegative edge weights
            self.INF = sys.maxint - 1

            if self.__distances.has_key(el1): #use cashed version
                if self.__distances[el1][el2] == self.INF:
                    return False
                else:
                    return True #self.__distances[el1][el2]

            else:                             #calculate and store in cache
                self.__distances[el1] = {}
                distances             = self.__distances[el1]

                for el in self.__elements:
                    distances[el] = self.INF # Initialize the distances

                distances[el1] = 0 # The distance of the source to itself is 0
                dist_to_unknown = distances.copy() # Select the next vertex to explore from this dict
                last = el1
            
                #while last != el2: #only check a single pair
                while len(dist_to_unknown) > 0: #check all pairs
                    # Select the next vertex to explore, which is not yet fully explored and which 
                    # minimizes the already-known distances.
                    next = min( [(v, k) for (k, v) in dist_to_unknown.iteritems()] )[1]
                    for n in self.downstreamConnectedElements(next): # n is the name of an adjacent vertex
                        distances[n] = min(distances[n], distances[next] + 1)
                        if dist_to_unknown.has_key(n):
                            dist_to_unknown[n] = distances[n]
                    last = next
                    if last in dist_to_unknown: # Delete the completely explored vertex
                        del dist_to_unknown[next]

                if distances[el2] == self.INF:
                    return False
                else:
                    return True #distances[el2]


    def connectedComponents(self):
        "Return a list of lists, each containing the elements of a connected component in the current splice graph --> list"

        # first, solve all-pairs shortest path problem for the splice graph using the Floyd-Warshall algorithm
        elements = self.__elements
        size     = len(elements)
        dist     = []  # dist[i][j] is the shortest distance between every two vertices, i and j,
                       # or self.INF if there is no such path
        pred     = []  # pred[i][j] == i if the shortest path goes directly from i to j
                       # Otherwise the shortest path from i to j is made up of the 
                       # shortest path from i to pred[i][j] + the edge from pred[i][j] to j
        # Initialization
        for i in xrange(size):
            dist.append([])
            pred.append([])
            for j in xrange(size):
              pred[i].append(self.INF)
              if i == j:
                dist[i].append(0)
              else:
                dist[i].append(self.INF)

        # Import graph
        for i in xrange(size):
            for j in xrange(size):
              if self.areDirectlyConnected(elements[i], elements[j]):
                  dist[i][j] = 1
                  pred[i][j] = i

        # Algorithm main loop
        for k in xrange(size):
            for i in xrange(size):
                for j in xrange(size):
                    if dist[i][k] == self.INF or dist[k][j] == self.INF:
                        continue
                    if dist[i][j] > dist[i][k] + dist[k][j]:
                        dist[i][j] = dist[i][k] + dist[k][j]
                        pred[i][j] = pred[k][j]

        #cache the calculated distances in self.__distances
        self.__distances = {}
        for i in xrange(size):
            self.__distances[elements[i]] = {}
            for j in xrange(size):
                self.__distances[elements[i]][elements[j]] = dist[i][j]

        # now, use the distances to extract the components
        inComp = {}
        nbComp = 0
 
        for el1 in self.__elements: #begin at the most upstream element

            if not inComp.has_key(el1):
                # get all connected elements
                els = [ el2 for el2 in self.__elements if el1 != el2 and \
                        (self.__distances[el1][el2] != self.INF or self.__distances[el2][el1] != self.INF) ]

                # check for elements that are already in a component
                comp = None
                for el2 in els:
                    if inComp.has_key(el2):
                        comp = inComp[el2]
                        break

                if comp is None:
                    # start a new component with current_element
                    comp = nbComp
                    nbComp += 1

                inComp[el1] = comp
                for el2 in els:
                    inComp[el2] = comp

        comps  = []
        for i in xrange(nbComp):
            comps.append([])
        for el in self.__elements:
            comps[inComp[el]].append(el)
        return comps


    def bridgeEdges(self):
        """return a list of (start,end)-tubles for all bridge edges in the splice graph
        bridge edges are edge that when removed, will disconnect the graph. they connect articulation points (cut vertices)

        WARNING: not throughoutly tested yet

        --> list"""

        v = { #dictionary of variables for passing by reference to DFS_visit_bridge
            'bridges' : [],
            'color'   : [], # depth-first-searching: white=not visited yet, gray=children being visited, black=visited (including children)
            'p'       : [],
            't_edges' : [],
            'd'       : [],
            'L'       : [],
            'time'    : 0,
            }

        # initialize
        for el in xrange(len(self.__elements)):
            v['color'].append('white')
            v['p'].append(None)
            v['t_edges'].append(0)
            v['d'].append(None)
            v['L'].append(None)

        # helper function for recursive (depth-first) traversal
        def DFS_visit_bridge(el1, v=v):
            v['color'][el1] = 'gray'
            v['time'] += 1
            v['d'][el1] = v['time']
            v['L'][el1] = v['time']
            for el2 in xrange(len(self.__elements)):
                if self.areDirectlyConnected(self.__elements[el1],self.__elements[el2]) \
                       and el2 != v['p'][el1]:  #NOT SURE HERE, originally said:  v X p[u] (X is unknown symbol)
                    if v['color'][el2] == 'white':
                        v['p'][el2] = el1
                        DFS_visit_bridge(el2, v)
                        v['L'][el1] = min(v['L'][el1],v['L'][el2])
                        if v['L'][el2] > v['d'][el1]:
                            v['bridges'].append((self.__elements[el1],self.__elements[el2]))
                    elif v['d'][el2] < v['d'][el1]:
                        v['L'][el1] = min(v['L'][el1],v['d'][el2])
            v['color'][el1] = 'black'
            v['time'] += 1
            #f[el1] = time


        # main loop
        for el in xrange(len(self.__elements)):
            if (v['color'][el] == 'white'):
                DFS_visit_bridge(el, v)

        #print self.name, "has bridges:", v['bridges']
        return v['bridges']
        

    def upstreamElements(self, el):
        "return a sorted list of all elements in splice graph UPSTREAM of a given element"
        try:
            index = self.__elements.index(el)
        except:
            return []
        else:
            return self.__elements[:index]
        

    def downstreamElements(self, el):
        "return a sorted list of all elements in splice graph DOWNSTREAM of a given element"
        try:
            index = self.__elements.index(el)
        except:
            return []
        else:
            return self.__elements[(index+1):]


    def upstreamConnectedElements(self, el):
        "return a sorted list of all elements in splice graph connected DIRECTLY UPSTREAM to a given element"

        connected = []

        for el2 in self.upstreamElements(el):
            if self.areDirectlyConnected(el2, el):
                connected.append(el2)

        return connected
        

    def downstreamConnectedElements(self, el):
        "return a sorted list of all elements in splice graph connected DIRECTLY DOWNSTREAM to a given element"

        connected = []

        for el2 in self.downstreamElements(el):
            if self.areDirectlyConnected(el, el2):
                connected.append(el2)

        return connected


    def is5ss(self, el):
        "return True if element is of 5ss type"
        if self.__types['5ss'].has_key(el):
            return True
        else:
            return False


    def is3ss(self, el):
        "return True if element is of 3ss type"
        if self.__types['3ss'].has_key(el):
            return True
        else:
            return False


    def isTSS(self, el):
        "return True if element is of TSS type"
        if self.__types['TSS'].has_key(el):
            return True
        else:
            return False


    def isTER(self, el):
        "return True if element is of TER type"
        if self.__types['TER'].has_key(el):
            return True
        else:
            return False

    def getType(self, el):
        "returns the first type associated with an element (priority: 5ss->3ss->TSS->TER)"
        if self.is5ss(el):
            return '5ss'
        elif self.is3ss(el):
            return '3ss'
        elif self.isTSS(el):
            return 'TSS'
        elif self.isTER(el):
            return 'TER'
        return None


    def getAllType(self, el):
        "returns all types associated with an element (priority: 5ss->3ss->TSS->TER)"
        types = []
        if self.is5ss(el):
            types.append('5ss')
        elif self.is3ss(el):
            types.append('3ss')
        elif self.isTSS(el):
            types.append('TSS')
        elif self.isTER(el):
            types.append('TER')
        return types

    def recalculateElementCoverage(self, el, type=None):
        "recalculates the type associated element coverage and stores in __types dict --> None"
        if type is None:
            elType = self.getType(el)
        else:
            elType = type
        newCoverage = 0
        if type == "5ss":
            for downstreamEl in self.downstreamConnectedElements(el):
                newCoverage += self.__connections[el][downstreamEl]
            self.__types["5ss"][el] = abs(newCoverage)
        elif type == "3ss":
            for upstreamEl in self.upstreamConnectedElements(el):
                newCoverage += self.__connections[upstreamEl][el]
            self.__types["3ss"][el] = abs(newCoverage)
        else:
            self.__inform("SpliceGraph.recalculateElementCoverage: unsupported type: %s" % type,4)
            
    def coverage(self, el1, el2=None):
        "return coverage for a element (pair), which is the nb of transcripts that evidence it --> int"

        coverage = 0

        if el2 is None:
            if self.is5ss(el1):
                coverage =+ self.__types['5ss'][el1]
            if self.is3ss(el1):
                coverage =+ self.__types['3ss'][el1]
            if self.isTSS(el1):
                coverage =+ self.__types['TSS'][el1]
            if self.isTER(el1):
                coverage =+ self.__types['TER'][el1]
        else:
            if self.__connections.has_key(el1) and self.__connections[el1].has_key(el2):
                coverage = self.__connections[el1][el2]

        return coverage

    def transcript_ids(self, el1, el2=None):
        "returns a list of transcript ids, the transcript evidence, for an element or connection --> list"
        ids = []

        if el2 is None:
            # outgoing transcripts
            if self.__transcripts.has_key(el1):
                for downstream_element in self.__transcripts[el1].keys():
                    for id in self.__transcripts[el1][downstream_element]:
                        if not id in ids:
                            ids.append(id)
            # incoming transcripts
            upstream_elements = self.upstreamConnectedElements(el1)
            for upstream_element in upstream_elements:
                if self.__transcripts.has_key(upstream_element):
                    if self.__transcripts[upstream_element].has_key(el1):
                        for id in self.__transcripts[upstream_element][el1]:
                            if not id in ids:
                                ids.append(id)
        else:
            if self.__transcripts.has_key(el1) and self.__transcripts[el1].has_key(el2):
                ids = self.__transcripts[el1][el2]

        return ids

    def elInxByTranscript(self, txId):
        "return the element indices (in self.__elements) that are supported by a given transcript (id) --> list"
        if not self.__path.has_key(txId):
            pathDict = {}
            for i in xrange(len(self.__elements)):
                if self.__transcripts.has_key(self.__elements[i]):
                    for j in xrange(i+1, len(self.__elements), 1):
                        if self.__transcripts[self.__elements[i]].has_key(self.__elements[j]) and \
                           txId in self.__transcripts[self.__elements[i]][self.__elements[j]]:
                            pathDict[i] = 1
                            pathDict[j] = 1
            self.__path[txId] = pathDict.keys()
            self.__path[txId].sort()
        return self.__path[txId]

    def transcriptLen(self, txId):
        "return the length of a (spliced) transcript in the splice graph --> int"
        path = self.elInxByTranscript(txId)
        l    = 0
        for i in xrange(0, len(path)-1, 2):
            if ( self.is3ss(self.__elements[i])   or self.isTSS(self.__elements[i])   )  and \
               ( self.is5ss(self.__elements[i+1]) or self.isTER(self.__elements[i+1]) ):
                l = l + abs( int(self.__elements[i+1].split(':')[1]) - int(self.__elements[i].split(':')[1]) )
        return l

    def read(self, filename, name=None):
        "read splice graph from file --> boolean"

        if not os.path.exists(filename):
            self.__inform("ERROR reading %s: %s" % (filename, 'File does not exist.'))
            raise Errors.ObjectInitError('SpliceGraph.read', "File '%s' does not exist" % filename)

        try:
            input = open(filename, 'r')
        except:
            self.__inform("ERROR reading %s: %s" % (filename, sys.exc_info()[1]))
            raise Errors.ObjectInitError('SpliceGraph.read', "Could not read file '%s': %s" % (filename, sys.exc_info()[1]))
        else:
            self.name          = name
            self.__connections = {}
            self.__types       = {}
            self.__transcripts = {}
            self.__elements    = []
            self.__strand      = None
            self.species       = "unknown"
            self.assembly      = "unknown"

            pattern = re.compile('^#species,\S+,assembly,\S+$')
            
            while 1:
                line = input.readline()
                if not line:
                    break

                elif pattern.search(line):
                    line = line.strip('\n#')
                    fields = line.split(',')
                    self.species  = fields[1]
                    self.assembly = fields[3]

                elif line.startswith('#'):
                    continue
                            
                else:
                    line = line.strip('\n')
                    (el, types, connections) = line.split('\t')

                    # do not set the name from the first element, but rather get it from the filename
                    # (see below). Otherwise, sg.name may not be consistent with the filename, e.g.
                    # in cases where the first element has been filtered from the graph, but the graph
                    # should still retain its original name
                    #if self.name is None:
                    #    self.name = el

                    for type in types.split(','):
                        (thisType, thisCount) = type.split('::')
                        if not self.__types.has_key(thisType):
                            self.__types[thisType] = {}
                        self.__types[thisType][el] = int(thisCount)


                    for connection in connections.split(','):
                        if len(connection) > 0:
                            (thisConnection, thisList) = connection.split('::')
                            thisList = thisList.split(';')
                            if not self.__connections.has_key(el):
                                self.__connections[el] = {}
                                self.__transcripts[el] = {}
                            self.__connections[el][thisConnection] = len(thisList)
                            self.__transcripts[el][thisConnection] = thisList

            # swap sign for introns
            for el1 in self.__connections.iterkeys():
                for el2 in self.__connections[el1].iterkeys():
                    # some splice graphs doesnt have 5'ss and 3'ss, so I added the following line
                    if self.__types.has_key('5ss') and self.__types.has_key('3ss'):
                        if self.__types['5ss'].has_key(el1) and self.__types['3ss'].has_key(el2):
                            self.__connections[el1][el2] = -self.__connections[el1][el2]
                        

            
            input.close()

            # do I have a valid name?
            if self.name is None:
                # try to get name from filename
                if filename.endswith('.sg'):
                    self.name = os.path.basename(filename)[0:-3]
                else:
                    self.name = "unknown"

            # build ordered list of all elements
            self.__buildElementsList()

            return True


    def store(self, filename=None):
        "store splice graph to a file --> boolean"

        if os.path.exists(filename):
            os.remove(filename)

        try:
            output = open(filename, 'w')
        except:
            self.__inform("ERROR writing %s: %s" % (filename, sys.exc_info()[1]))
            raise IOError('SpliceGraph.store', "Could not write file '%s': %s" % (filename, sys.exc_info()[1]))
        else:
            output.write('#species,'+self.species+',assembly,'+self.assembly+'\n')
            output.write('#'+'\t'.join(['element_id','element_types','connected_to'])+'\n')

            for el in self.__elements:
                types = {}
                for type in self.__types.iterkeys():
                    if self.__types[type].has_key(el):
                        types[type] = self.__types[type][el]
                if self.__connections.has_key(el):
                    sortedKeys = self.__connections[el].keys()
                    sortedKeys.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])), reverse=(self.__strand == '-'))
                    output.write('\t'.join([el,                                                                     \
                                            ",".join("%s::%i" % (k, v) for k, v in types.items()),                  \
                                            ",".join("%s::%s" % (k, ';'.join(self.__transcripts[el][k])) for k in sortedKeys), \
                                            ])+'\n')
                else:
                    output.write('\t'.join([el,                                                                     \
                                            ",".join("%s::%i" % (k, v) for k, v in types.items()),                  \
                                            "",                                                                     \
                                            ])+'\n')

            output.close()
            os.chmod(filename, 0660)

            return True


    def visualize(self, filename=None, emph=None):
        "visualize splice graph to a PS file --> boolean"

        if os.path.exists(filename):
           os.remove(filename)
	
	

        if emph is not None and not isinstance(emph, dict):
           self.__inform("WARNING: Ignoring 'emph'. Need a dictionary for emphasizing edges.\n")
           emph = None

        sg = pydot.Dot(graph_name=self.name.replace(':','_').replace('+','plus').replace('-','minus'))

        #add nodes
        elNodes = {}
        for el in self.__elements:
            elNodes[el] = pydot.Node(name=self.__getName(el))#,                   \
                                     #label=self.__getName(el).replace('_',':'), \
                                     #fontsize=self.__defaultFontsize,           \
                                     #color=self.__getColor(el),                 \
                                     #shape=self.__getShape(el),                 \
                                     #fixedsize='true',                          \
                                     #style=self.__getStyle(el))
            sg.add_node(elNodes[el])

        #add edges
        for Src in self.__connections.iterkeys():
            for Dst in self.__connections[Src].iterkeys():
                sg.add_edge(pydot.Edge(src=elNodes[Src].get_name(),                           \
                                       dst=elNodes[Dst].get_name(),color=self.__getColor(Src, Dst, emph)))                           \
                                      # name=self.__getName(Src, Dst),                         \
                                       #label=self.__getName(Src, Dst, emph).replace('_',':'), \
                                       #fontsize=self.__defaultFontsize                       \
                                       #color=self.__getColor(Src, Dst, emph),                 \
                                       #style=self.__getStyle(Src, Dst),                       \
                                       #arrowType='normal',                                    \
                                       #dir='forward',                                         \
                                       #)
                            #)

        #store PS
        #sg.write-raw(filename+'.dot', prog='dot')
        try:
            sg.write_ps2(filename, prog='dot')
	   
        except:
            raise Errors.ObjectInitError('SpliceGraph.visualize', "Could not write file '%s': %s" % (filename, sys.exc_info()[1]))
        else:
            os.chmod(filename, 0660)


    def __getName(self, el1, el2=None, emph=None):
        "return name for drawing a node/edge"

        name = None

        if el2 is None: #name for a node
            #what type(s) are associated with el1?
            types = []
            for type in self.__types.iterkeys():
                if self.__types[type].has_key(el1):
                    types.append(type)

            (chr, coord, strand) = el1.split(':')
            name = "%s\\n%s:%s\\n%s" % (chr, coord, strand, ",".join(types))
            #name = "%s" % el1

        else:           #name for an edge
            type = None
            label = None
            if emph is not None and emph.has_key(el1) and emph[el1].has_key(el2):
                label = emph[el1][el2]
            elif self.__connections[el1][el2] > 0:
                label = "ex"
            else:
                label = "in"

            name = "%s:%i\\n%i" % (label, abs(self.__connections[el1][el2]),
                                   abs(int(el1.split(':')[1]) - int(el2.split(':')[1]))+1 )

        #replace ':' by '_'
        name = name.replace(':', '_')


        return name


    def __getStyle(self, el1, el2=None):
        "return style for drawing a node/edge"

        style = None

        if el2 is None: #style for a node
            style = 'filled'

        else:           #style for an edge
            if self.__connections[el1][el2] > 0:
                style = 'bold'
            else:
                style = 'solid'

        return style


    def __getShape(self, el):
        "return the shape for drawing a node"

        #what type(s) are associated with el?
        types = []
        for type in self.__types.iterkeys():
            if self.__types[type].has_key(el):
                types.append(type)

        #return default shape
        shape = None
        if len(types) == 1:
            shape = self.__defaultShapes[types[0]]
        elif len(types) > 1:
            shape = self.__defaultShapes['multi']

        return shape


    def __getColor(self, el1, el2=None, emph=None):
        "return the color for drawing a node"

        color = None

        if el2 is None: #node color
            #what type(s) are associated with el1?
            types = []
            for type in self.__types.iterkeys():
                if self.__types[type].has_key(el1):
                    types.append(type)

            #return default shape
            if len(types) == 1:
                color = self.__defaultColors[types[0]]
            elif len(types) > 1:
                color = self.__defaultColors['multi']

        else:           #edge color
            if emph is not None and emph.has_key(el1) and emph[el1].has_key(el2):
                color = self.__defaultColors['emph']
            elif self.__connections[el1][el2] > 0:
                color = self.__defaultColors['exon']
            else:
                color = self.__defaultColors['intron']

        return color


    def __repr__(self):
        "generate a string representation of a SpliceGraph object (automatically called, e.g. by print) --> string"

        return repr({'name':self.name,'species':self.species,'assemby':self.assembly,
                     'connections':self.__connections, 'types':self.__types})


    def __inform(self, txt, level=0):
        "output progress information to self.log if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            self.log.write(txt)

    # SPLICE GRAPH FILTERING METHODS
    ########################################
    def remove_by_coverage_and_ss(self, coverage=None):
        "Copies SpliceGraph and removes elements for which coverage is below set threshold and ss is non-canonical --> SpliceGraph"

        if coverage is None:
            coverage = int(self.conf[self.species]["FilterMinCoverage"])
            #self.__inform("WARNING no coverage information given to SpliceGraph.remove_by_coverage_and_ss. default (%i) used." % coverage,2)
        
        newSpliceGraph = copy.deepcopy(self)
        
        elements = copy.deepcopy(newSpliceGraph.allElements())

        nrElements = len(elements)
        nrRemovedElement = 0

        #print newSpliceGraph.__elements

        for element in elements:
    
            if ( newSpliceGraph.is3ss(element) or newSpliceGraph.is5ss(element) ) and newSpliceGraph.coverage(element) < coverage:
                
                ##test splice site, if canonical returns True otherwise False
                #splicesite = SpliceSites.SpliceSite(species=self.species, assembly=self.assembly, node=element, type=newSpliceGraph.getType(element))
                #if splicesite.isSScanonical() is False:
                splicesite = SS(element, self.getType(element), species=self.species, assembly=self.assembly)
                if not splicesite.isCanonical():
                
                    # remove all edges to and from element
                    upstreamConnectedElements = newSpliceGraph.upstreamConnectedElements(element)
                    downstreamConnectedElements = newSpliceGraph.downstreamConnectedElements(element)

                    # Remove Element
                    newSpliceGraph.removeElement(element)
                    nrRemovedElement += 1
                      
                    for upConnection in upstreamConnectedElements:

                        # check if the now unconnected upConnection became isolated
                        if len( newSpliceGraph.upstreamConnectedElements(upConnection) ) == 0 and \
                           len( newSpliceGraph.downstreamConnectedElements(upConnection) ) == 0:
                            newSpliceGraph.removeElement(upConnection)
                            nrRemovedElement += 1
                            self.__inform("Removed element: %s, because it got isolated (upconnection)\nRemaining elements: %s\n" \
                                          % (upConnection, str(newSpliceGraph.__elements)), 5) 
                        
                    for downConnection in downstreamConnectedElements:

                        # check if the now unconnected element is isolated
                        if len( newSpliceGraph.upstreamConnectedElements(downConnection) ) == 0 and \
                           len( newSpliceGraph.downstreamConnectedElements(downConnection) ) == 0 :
                            newSpliceGraph.removeElement(downConnection)
                            nrRemovedElement += 1
                            self.__inform("Removed element: %s, because it got isolated (downconnection)\nRemaining elements: %s\n" \
                                          % (downConnection, str(newSpliceGraph.__elements)), 5) 

        #print '%i (%.1f %%) elements removed.' % (nrRemovedElement, nrRemovedElement/float(nrElements)*100)
        #print newSpliceGraph
        #print newSpliceGraph.__elements

        return newSpliceGraph


    def tweak_non_canonical(self, nbBases=3):
        """Copies SpliceGraph and tries to optimize elements with non-canonical (non GT->AG or GC->AG)
        splice sites by shifting up to nbBases across the intron

        if a non-canonical splice site can be obtained by this, the copied splice graph is adjusted:
           - flanking elements are renamed (but the associated types and transcripts are retained)
           - all connectivity is adjusted

        --> SpliceGraph
        """

        newSpliceGraph = copy.deepcopy(self)
        introns = copy.deepcopy(self.allIntrons())
       
        position_pairs = [] # to build a position tuple with complmentary positions shifts in 5ss and 3ss, 
        l3ss = SS.SSlen["3ss"] # checks the length of a 3ss from SS class

        for i in range(1,nbBases+1,1):
            if i == 1:
                position_pairs.append((-i,1,l3ss-i-1,l3ss-i)) # Left shift
            else:
                position_pairs.append((-i,-i+1,l3ss-i-1,l3ss-i)) # Left shift
            # right shift
            position_pairs.append((i+1,i+2,l3ss+i-1,l3ss+i))
       
        for intron in introns:

            ss5 = SS(intron[0], "5ss", species=self.species, assembly=self.assembly)
            ss3 = SS(intron[1], "3ss", species=self.species, assembly=self.assembly)

            # sanity check
            # if perhaps the splice sites have already been moved for these ss, skip to next intron
            
            if not intron[0] in newSpliceGraph.allElements() and \
                   not intron[1] in newSpliceGraph.allElements():
                print "skipped intron:", intron
                continue
   
            if not ss5.isCanonical() and not ss3.isCanonical():
                # non canonical intron
                
                for pos5S, pos5E, pos3S, pos3E  in position_pairs:
                    # look for a canonical ss within reach
                    if ( ss5[pos5S:pos5E] == "GT" or ss5[pos5S:pos5E] == "GC" ) and ss3[pos3S:pos3E] == "AG":

                        # the sequence to replace from one exon end to the other must be identical
                        if (pos5S > 0 and  ss5[1:pos5S-1] == ss3[l3ss+1:pos3E] ) or \
                               pos5S < 0 and ss5[pos5S:-1] == ss3[pos3E+1:l3ss]:
                            
                            # 1. Calculate new coordiates for tweaked elements
                            chr, coord5ss, strand = intron[0].split(':')
                            chr, coord3ss, strand = intron[1].split(':')
                            
                            if strand == '+':
                                if pos5S > 0:
                                    new5ss = "%s:%i:%s" % (chr,int(coord5ss) + int(pos5S)-1,strand)
                                else:
                                    new5ss = "%s:%i:%s" % (chr,int(coord5ss) + int(pos5S),strand)
                                new3ss = "%s:%i:%s" % (chr,int(coord3ss) - (l3ss - pos3E),strand)
                                
                            else:
                                if pos5S > 0:
                                    new5ss = "%s:%i:%s" % (chr,int(coord5ss) - int(pos5S)+1,strand)
                                else:
                                    new5ss = "%s:%i:%s" % (chr,int(coord5ss) - int(pos5S),strand)
                                new3ss = "%s:%i:%s" % (chr,int(coord3ss)+ (l3ss - pos3E),strand)

                            # the original elements
                            old5ss, old3ss = intron
                            self.__inform("Tweaked at: %s, %s --> %s, %s\n" % (old5ss,old3ss,new5ss,new3ss),4)
                            #print "Tweaked at: %s, %s --> %s, %s\n" % (old5ss,old3ss,new5ss,new3ss)

                            # This is a quick fix for the pesky chr 14... it skips it to be filtered out later.
                            #if new5ss = X, don't move jxn.
                            newSpliceGraph.moveJunction(old5ss, new5ss, old3ss, new3ss)

               
        return newSpliceGraph

    def remove_non_canonical(self):
        "Copies SpliceGraph and removes elements with non-canonical (non GT->AG or GC->AG) splice sites --> SpliceGraph"

        newSpliceGraph = copy.deepcopy(self)
        elements = copy.deepcopy(newSpliceGraph.allElements())

        nrRemovedElement = 0

        for element in elements:
    
            if newSpliceGraph.is3ss(element) or newSpliceGraph.is5ss(element):

                if not element in newSpliceGraph.__elements:
                    self.__inform("WARNING: element not in __elements list: %s\n" % element , 3)
                    continue
                
                ##test splice site, if canonical returns True otherwise False
                #splicesite = SpliceSites.SpliceSite(species=self.species, assembly=self.assembly, node=element, type=newSpliceGraph.getType(element))
                #if splicesite.isSScanonical() is False:
                splicesite = SS(element, self.getType(element), species=self.species, assembly=self.assembly)
                if not splicesite.isCanonical():
                
                    # remove all edges to and from element
                    upstreamConnectedElements = newSpliceGraph.upstreamConnectedElements(element)
                    downstreamConnectedElements = newSpliceGraph.downstreamConnectedElements(element)

                    # Remove Element
                    newSpliceGraph.removeElement(element)
                    nrRemovedElement += 1
                      
                    for upConnection in upstreamConnectedElements:

                        # check if the now unconnected upConnection became isolated
                        if len( newSpliceGraph.upstreamConnectedElements(upConnection) ) == 0 and \
                           len( newSpliceGraph.downstreamConnectedElements(upConnection) ) == 0:
                            newSpliceGraph.removeElement(upConnection)
                            nrRemovedElement += 1
                            self.__inform("Removed element: %s, because it got isolated (upconnection)\nRemaining elements: %s\n" \
                                          % (upConnection, str(newSpliceGraph.__elements)), 5) 
                        
                    for downConnection in downstreamConnectedElements:

                        # check if the now unconnected element is isolated
                        if len( newSpliceGraph.upstreamConnectedElements(downConnection) ) == 0 and \
                           len( newSpliceGraph.downstreamConnectedElements(downConnection) ) == 0 :
                            newSpliceGraph.removeElement(downConnection)
                            nrRemovedElement += 1
                            self.__inform("Removed element: %s, because it got isolated (downconnection)\nRemaining elements: %s\n" \
                                          % (downConnection, str(newSpliceGraph.__elements)), 5)

                    # update element list
                    newSpliceGraph.__buildElementsList()

        #print '%i (%.1f %%) elements removed.' % (nrRemovedElement, nrRemovedElement/float(len(elements))*100)
        #print newSpliceGraph
        #print newSpliceGraph.__elements

        return newSpliceGraph


    def get_library_keyword_synonyms(self, keyword):
        "return a dictionary with keyword-keys that are hierarchically below a given keyword in an ontology --> dict"

        if self.__keysyns.has_key(keyword):
            return self.__keysyns[keyword]

        else:
            keywords = {keyword:1}

            if self.conf is None: #configuration was not set in __init__ -> get defaults
                self.conf = configuration.Configuration()

            try:
                f = open("%s/vocab.cgap" % self.conf[self.species]["MetaDataStoragePath"], 'r')
            except:
                self.__inform("WARNING: could not open %s/vocab.cgap\n" % self.conf[self.species]["MetaDataStoragePath"])
            else:
                level = None
                for line in f:
                    if line.startswith('Taxonomy of'):
                        continue
                    else:
                        fields = line.split('\t')
                        if len(fields) == 4:
                            if level is None and fields[2].strip(' ') == keyword: #found keyword
                                level = int(fields[0])

                            elif not level is None and level < int(fields[0]):    #found child of keyword
                                keywords[fields[2].strip(' ')] = 1

                            elif not level is None:                               #no more children
                                break
                f.close()

            self.__keysyns[keyword] = keywords
            return keywords


    def remove_by_library_keyword(self, keyword):
        "Copies SpliceGraph and removes elements from libraries that match a certain keyword --> SpliceGraph"

        # find all the synonymous keywords to keyword using vocab.cgap (see MetaDataFiles in configuration.txt)
        keywords = self.get_library_keyword_synonyms(keyword)

        self.__inform('\tkeyword(s) used: %s\n' % ', '.join(keywords.keys()),2)

        # copy splice graph and remove transcripts
        newSpliceGraph = copy.deepcopy(self)
        for tid in newSpliceGraph.allTranscripts():
            lid = newSpliceGraph.getLibraryID(tid)
            for libkeyword in newSpliceGraph.getLibraryKeywords(lid):
                if keywords.has_key(libkeyword):
                    newSpliceGraph.removeTranscript(tid)
                    break

        return newSpliceGraph


    def remove_by_tissue_expression(self, tissue):
        "Copies SpliceGraph and removes elements with available expression data that are not expressed in a given tissue"

    # HEAVEST_BUNDLE ASSOCIATED METHODS (Lee C., Bioinformatics, 2003 and Xing et al., 2004, Genome Research)
    #########################################################################################################
    def HB_generate_consensi(self):
        "generate consensus sequences by heaviest_bundle algorithm (Lee, 2003) --> list"

        S = self.allTranscripts()         # sequences to be bundled
        l = []                            # lengthes of transcripts in S
        for s in S: l.append(self.transcriptLen(s))
        B = []                            # bundling vector: B[seq_index] = bundle_index
        for s in S: B.append(0)
        #print "\tHB_edgeWeights()"
        w = self.HB_edgeWeights(S, l, B)  # get edge weights
        b = 1                             # current consensus (iteration)
        C = []                            # consensi

        while (0 in B):
            #sequences left to be bundled
            #print "B=%s" % str(B)
            #print "w=%s" % str(w)
            #print "starting iteration %i:" % b

            #print "\tHB_heaviest_bundle()"
            P = self.HB_heaviest_bundle(w)

            #print "\tHB_create_sequence_on_path()"
            C.append(self.HB_create_sequence_on_path(P, "%s_consensus%i" % (self.name, b)))

            #print "\tHB_add_sequences_to_bundle()"
            I = self.HB_add_sequences_to_bundle(S, B, P, b)
            if len(I) == 0:
                break
            C[-1]['txs'].extend(I)

            #self.HB_rescale_weights(I, w)
            # isoform generate mode: weights as initial weights, but new template
            #print "\tHB_edgeWeights()"
            w = self.HB_edgeWeights(S, l, B)

            b  = b + 1

        return C

    def HB_edgeWeights(self, S, l, B):
        "return edge weight for current graph --> dict"
        # isoform generation mode:
        # 1. select longes unbundled sequence as template
        # 2. set weight of template and overlapping sequences to 1000

        # select template
        lmax = 0
        tmpl = None
        for i in xrange(len(S)):
            if l[i] > lmax and B[i] == 0:
                lmax = l[i]
                tmpl = S[i]
        #print "\tselected %s (%i nt)" % (tmpl, lmax)

        # calcuate weights
        w = {}
        for n1 in self.__connections.iterkeys():
            w[n1] = {}
            for n2 in self.__connections[n1].iterkeys():
                w[n1][n2] = abs(self.__connections[n1][n2])
                if tmpl in self.__transcripts[n1][n2]:
                    w[n1][n2] = w[n1][n2] * 1000
        return w

    def HB_heaviest_bundle(self, w):
        "dynamic programming reconstruction of the max. likelyhood parse according to weights w --> path (dict)"
        t = [] # traceback information: t[i] = j (predecessor of i is j)
        #print "\t\tHB_traverse_hb()"
        r = self.HB_traverse_hb(w, t)

        if self.__connections.has_key(self.__elements[r]):
            # r has outgoing connections --> internal termination
            #print "\t\tHB_branch_completion()"
            r = self.HB_branch_completion(r, w, t)

        #print "\t\tHB_traceback()"
        return self.HB_traceback(r, t)

    def HB_traverse_hb(self, w, t):
        "traverse graph --> index of last node of traversal (int)"
        s    = []    # dynamic programming scores of nodes (index corresponds to self.__elements)
        smax = 0     # maximum score
        r    = None  # node with smax

        for i in xrange(len(self.__elements)):
            s.append(0)
            t.append(None)
            p = self.HB_best_predecessor(w, i, s)
            if p != None:
                t[i] = p
                s[i] = s[p] + w[self.__elements[p]][self.__elements[i]]
                if s[i] > smax:
                    r = i
        return r

    def HB_best_predecessor(self, w, i, s):
        "return index of best predecessor of i --> int"
        wmax = 0
        smax = 0
        p    = None
        for q in xrange(len(self.__elements)):
            if ( self.__connections.has_key(self.__elements[q])                         ) and \
               ( self.__connections[self.__elements[q]].has_key(self.__elements[i])     ) and \
               ( ( w[self.__elements[q]][self.__elements[i]] > wmax )              or \
                 ( w[self.__elements[q]][self.__elements[i]] == wmax and s[q] > smax )  ):
                # 1. 'q' is a predecessor of 'i'
                # 2. w_qi > w_max  OR  w_qi == w_max AND s_q > s_max
                # 3. as initially: smax=0, just 'q' with s[q]>0 (HB_branch_completion)
                wmax = w[self.__elements[q]][self.__elements[i]]
                smax = s[q]
                p    = q
        return p

    def HB_branch_completion(self, r, w, t):
        "branch completion to resolve internal termination --> index of last node of path (int)"
        s    = []    # dynamic programming scores of nodes (index corresponds to self.__elements)
        smax = 0     # maximum score
        r2   = None  # node with smax

        # exclude all nodes upstream of r from being valid starting nodes (set s[i] to neg. value)
        for i in xrange(r+1):
            s.append(-1)
        s[r] = 0

        # and repeat HB_traverse_hb downstream of r
        for i in xrange(r+1, len(self.__elements), 1):
            s.append(0)
            t[i] = None
            p = self.HB_best_predecessor(w, i, s)
            if p != None:
                t[i] = p
                s[i] = s[p] + w[self.__elements[p]][self.__elements[i]]
                if s[i] > smax:
                    r2 = i
        return r2

    def HB_traceback(self, r, t):
        "global traceback of max. likelyhood path --> list"
        P = [r] # path
        i = r
        while t[i] != None:
            P.append(t[i])
            i = t[i]
        P.reverse()
        return P

    def HB_create_sequence_on_path(self, P, id=None):
        "create sequence corresponding to path P --> dict"
        seqobj = {'seq'  : '', # sequence string
                  'id'   : id, # sequence identifier
                  'desc' : '', # sequence description
                  'path' : [], # path (list of (exonStart, exonEnd)-tuples)
                  'lens' : [], # exon lengthes
                  'txs'  : [], # transcripts (ids) associated to sequence (heaviest bundle)
                  }
        fetcher = GenomeFetch.GenomeFetch(species=self.species,assembly=self.assembly)
        for i in xrange(0, len(P)-1, 2):
            if ( self.is3ss(self.__elements[P[i]])   or self.isTSS(self.__elements[P[i]])   )  and \
               ( self.is5ss(self.__elements[P[i+1]]) or self.isTER(self.__elements[P[i+1]]) ):
                (chr, start, strand) = self.__elements[P[i]].split(':')
                end                  = self.__elements[P[i+1]].split(':')[1]
                seqobj['seq']  = seqobj['seq']  + fetcher.get_seq_from_to(chr, int(start), int(end), strand, fLeaveOpen=True)
                seqobj['desc'] = seqobj['desc'] + "%s:E:%s;" % (self.__elements[P[i]], self.__elements[P[i+1]])
                seqobj['path'].append((self.__elements[P[i]], self.__elements[P[i+1]]))
                seqobj['lens'].append(abs(int(end)-int(start))+1)
        return seqobj

    def HB_add_sequences_to_bundle(self, S, B, P, b):
        "add sequences from S that support path P to bundle b --> sequences in bundle b (list)"
        I = []
        for k in xrange(len(B)):
            if B[k]==0 and self.HB_inclusion_rule_SG(S[k], P)==True:
                B[k] = b
                I.append(S[k])
        return I

    def HB_inclusion_rule_SG(self, txId, consPath):
        "does a transcript fit to path (splice graph based)? --> bool"
        # inclusion rules:
        # 0. n: nb of elements in transcript path txPath
        # 1. if n>3, then x=1,y=1, else x=0,y=0
        # 2. no more than x differing terminal element on each side
        # 3. no more than y differing internal element
        (maxStart, maxEnd, maxInternal) = (1, 1, 1)

        consDict = {}
        for item in consPath: consDict[item] = 1

        txPath = self.elInxByTranscript(txId)
        if len(txPath) <= 3:
            (maxStart, maxEnd, maxInternal) = (0, 0, 0)
        txDict = {}
        for item in txPath: txDict[item] = 1
        (diffStart, diffEnd, diffInternal) = (0, 0, 0)

        common = [item for item in txPath if consDict.has_key(item)]
        if len(common) == 0:
            return False
        else:
            # diffStart
            for t in xrange(0, len(txPath), 1):
                if consDict.has_key(txPath[t]):
                    diffStart = t
                    break
            # diffEnd
            for t in xrange(len(txPath)-1, diffStart-1, -1):
                if consDict.has_key(txPath[t]):
                    diffEnd = len(txPath) - 1 - t
                    break
            # diffInternal
            if len(txPath) > diffStart+diffEnd+1:
                for t in xrange(diffStart, len(txPath)-diffEnd, 1):
                    if not consDict.has_key(txPath[t]): # insert in txPath
                        diffInternal = diffInternal + 1
                for p in xrange(consPath.index(txPath[diffStart]), consPath.index(txPath[len(txPath)-1-diffEnd]), 1):
                    if not txDict.has_key(consPath[p]):     # insert in consPath
                        diffInternal = diffInternal + 1
                diffInternal = diffInternal - diffEnd

            #print "\t\t\tconsPath=%s" % str(consPath)
            #print "\t\t\t  txPath=%s" % str(txPath)
            #print "\t\t\t         diffStart=%i, diffEnd=%i, diffInternal=%i" % (diffStart, diffEnd, diffInternal))

            # alternatively, use the number of bases that overlap, and require minimum fraction
            return (diffStart <= maxStart  and  diffEnd <= maxEnd  and  diffInternal <= maxInternal)

    def HB_inclusion_rule_ALIGN(self, txId, P):
        "does a transcript fit to path (alignment based)? --> bool"
        # get sequences
        txPath  = self.elInxByTranscript(txId)
        txSeq   = self.HB_create_sequence_on_path(txPath)
        consSeq = self.HB_create_sequence_on_path(P)

        # store in temporary files and align
        fh = open("%s/tmp1.fa" % self.__tmpdir, 'w')
        fh.write(">txSeq\n%s\n" % txSeq)
        fh.close()
        fh = open("%s/tmp2.fa" % self.__tmpdir, 'w')
        fh.write(">consSeq\n%s\n" % consSeq)
        fh.close()
        os.system("blat -t=dna -q=dna -tileSize=8 -oneOff=1 -minIdentity=50 -noHead %s/tmp1.fa %s/tmp2.fa %s/tmp.psl &>/dev/null" \
                  % (self.__tmpdir, self.__tmpdir, self.__tmpdir))

        # parse alignment and get quality, then clean up
        (percIdTotal, lenInternalGaps, lenTerBranchesQ, lenTerBranchesT, percIdOverlap) = (None, None, None, None, None)
        fh = open("%s/tmp.psl" % self.__tmpdir, 'r')
        line = fh.readline()
        line = line.strip("\n")
        if line != "":
            fields = line.split("\t")
            # fields[0] = match          fields[9] = Q name         fields[18] = blockSizes
            #        1  = mismatch              10 = Q size                19  = qStarts
            #        2  = repmatch              11 = Q start               20  = tStarts
            #        3  = N's                   12 = Q end
            #        4  = Q gap count           13 = T name
            #        5  = Q gap bases           14 = T size
            #        6  = T gap count           15 = T start
            #        7  = T gap bases           16 = T end
            #        8  = strand                17 = block count
            percIdTotal = int(fields[0]) / ( int(fields[10])<int(fields[14]) and int(fields[10]) or int(fields[14]) )
            lenInternalGaps = int(fields[5])
            lenTerBranchesQ = int(fields[10]) - (int(fields[12])-int(fields[11]))
            lenTerBranchesT = int(fields[14]) - (int(fields[16])-int(fields[15]))
            perIdOverlap = int(fields[0]) / ( int(fields[10])-lenTerBranchesQ < int(fields[14])-lenTerBranchesT and \
                                              int(fields[10])-lenTerBranchesQ or \
                                              int(fields[14])-lenTerBranchesT )
            print "%s" % line
            print "\t\t\t%.2f\t%i\t%i" % (percIdTotal, lenInternalGaps, lenTerBranchesT)
        fh.close()
        os.system("rm %s/tmp1.fa %s/tmp2.fa %s/tmp.psl" % (self.__tmpdir, self.__tmpdir, self.__tmpdir))

        # evaluate alignment quality
        return (percIdTotal >= 0.9 and lenInternalGaps <= 5 and lenTerBranchesT <= 20)

    def HB_rescale_weights(self, I, w):
        "rescale weights for next iteration of HB --> None"
        # remove the bundled sequences (I) from the coverage (w)
        for txId in I:
            path = self.elInxByTranscript(txId)
            for i in xrange(0, len(path)-1, 2):
                w[self.__elements[path[i]]][self.__elements[path[i+1]]] = w[self.__elements[path[i]]][self.__elements[path[i+1]]] - 1

    def HB_get_protein_isoforms(self, minLen=50, checkNMD=True):
        """convenience method that constructs consensi, translates and filters protein isoform sequences
           applied filtering rules:
               1. select longest ORF with valid start/stop codons
               2. filter out peptides shorter than 'minLen' amino acids
               3. filter out NMD candidates (if 'checkNMD' is True)

           --> list"""
        ntC = self.HB_generate_consensi()
        #print "found %i consensus sequences" % len(ntC)
        aaC = []
        for i in xrange(len(ntC)):
            #print "\ttranslating consensus %i..." % (i+1),
            self.HB_translate_longest_ORF(ntC[i])
            #print "\tdone (%i aa, NMD:%s)" % (len(ntC[i]['seq']), (self.HB_NMD_candidate(ntC[i]) and "yes" or "no"))
            if len(ntC[i]['seq'])>=minLen and (checkNMD==False or not self.HB_NMD_candidate(ntC[i])):
                self.HB_calc_exon_phases(ntC[i])
                aaC.append(ntC[i])
        return aaC

    def HB_translate_longest_ORF(self, seq):
        "translate sequence in all three frames and return longest ORF translation starting at Met --> None"
        seq['ntSeq'] = seq['seq']
        seq['seq'] = ''
        for f in [1, 2, 3]:
            frameSeqs = self.HB_translate(seq['ntSeq'], f, True)
            for frameSeq in frameSeqs:
                if len(frameSeq[0])>len(seq['seq']):
                    seq['frame']    = f
                    seq['seq']      = frameSeq[0]
                    seq['CDSstart'] = frameSeq[1]
                    seq['CDSend']   = frameSeq[2]

    def HB_translate(self, seq, frame=1, split=False):
        "translate sequence using standard genetic code --> str (split==False), list (split==True)"
        seq   = seq.upper()
        trans = ''
        for i in xrange(frame-1, len(seq) - (len(seq)-(frame-1))%3, 3):
            # i: starting position of each codon
            codon = seq[i:i+3]
            if self.__table.has_key(codon):
                trans = trans + self.__table[codon]
            else:
                trans = trans + 'X'

        if split==True:
            # return list of translated ORFs, truncated to [Start....Stop)
            # (this will skip ORFs not ending with Stop)
            translist = [] # list elements are 3-tuples: (seqstr, start, stop), start and stop are 1-based,inclusive
            inCDS     = False
            (current, start, end) = ('', None, None)
            for a in xrange(len(trans)):
                codon = seq[((frame-1)+a*3):((frame-1)+a*3+3)]
                if inCDS==False and self.__start_codons.has_key(codon):
                    inCDS = True
                    start = (frame-1)+a*3+1
                elif inCDS==True and self.__stop_codons.has_key(codon):
                    inCDS = False
                    end   = (frame-1)+a*3
                    if len(current)>0:
                        translist.append((current, start, end))
                        (current, start, end) = ('', None, None)
                if inCDS==True:
                    current = current + trans[a]
            return translist
        else:
            return trans

    def HB_NMD_candidate(self, seq):
        "use 50-base-rule to identify NMD candidates --> bool"
        posLastJunction = sum(seq['lens'][0:-1])+1 #1-based,inclusive
        return ( seq['CDSend'] < posLastJunction-50 )

    def HB_calc_exon_phases(self, seq):
        "calculate exon phases based on CDSstart,CDSend and exon lengthes --> None"
        #print "CDSstart: %i, CDSend: %i" % (seq['CDSstart'], seq['CDSend'])
        cumTotLen = 0
        cumCDSLen = 0
        seq['exPhases'] = [] # list of 2-element-tuples: (exStartPhase, exEndPhase)
        for i in xrange(len(seq['lens'])):
            exStart = cumTotLen + 1 # coords are 1-based,inclusive
            exEnd   = exStart + seq['lens'][i] - 1
            startPhase = -1
            endPhase   = -1
            if exEnd >= seq['CDSstart'] and exStart <= seq['CDSend']: # exon overlaps CDS
                partCDSstart = max(seq['CDSstart'], exStart)
                partCDSend   = min(seq['CDSend'],   exEnd)
                partCDSlen   = partCDSend - partCDSstart + 1
                if exStart >= seq['CDSstart']:
                    startPhase = cumCDSLen % 3
                if exEnd <= seq['CDSend']:
                    endPhase = (cumCDSLen + partCDSlen) % 3
                cumCDSLen = cumCDSLen + partCDSlen
                #print "exon: %i-%i, partCDS: %i-%i, phases: (%i, %i)" % (exStart, exEnd, partCDSstart, partCDSend, startPhase, endPhase)
            #else:
                #print "exon: %i-%i, partCDS: ?-?, phases: (%i, %i)" % (exStart, exEnd, startPhase, endPhase)
            seq['exPhases'].append((startPhase, endPhase))
            cumTotLen = cumTotLen + seq['lens'][i]

# SEQUENCE CLASS (HANDLES EXONS/INTRONS/SS
##########################################
class NtSequence:
    "NtSequence: represents a nucleotide sequence defined by a pair of SpliceGraph elements and associated information"

    fetcher         = None   #GenomeFetcher instance
    defaultSpecies  = 'hsa'
    defaultAssembly = 'hg18'
    # # #defaultAssembly = 'hg17'
    width           = 60

    def __init__(self, el1, el2, name=None, desc=None, species=None, assembly=None):
        "instantiate Sequence --> Sequence"

        # get location from elements
        (chr, start, strand) = el1.split(':')
        end = el2.split(':')[1]
        start = int(start)
        end = int(end)

        # make sure that start<end
        if start > end:
            start, end = end, start

        # class attributes
        self.name     = name             # sequence name (used as id in Fasta header)
        self.desc     = desc             # sequence description (used in Fasta header)
        self.chr      = chr              # chromosome
        self.start    = int(start)       # sequence start position (coordinate system is biological, i.e. 1-based inclusive, start<end)
        self.end      = int(end)         # sequence end position
        self.strand   = strand           # sequence strand
        self.species  = species          # species/assembly for fetching of sequence
        self.assembly = assembly
        self.seq      = self.__get_seq(self.chr, self.start, self.end, self.strand) # exon sequence

        if self.species is None:
            self.species = self.defaultSpecies
        if self.assembly is None:
            self.assembly = self.defaultAssembly


    def __get_seq(self, chr, start, end, strand):
        "get the sequence of the exon using GenomeFetch --> str"

        if self.fetcher is None:
            self.fetcher = GenomeFetch.GenomeFetch(species=self.species,assembly=self.assembly)

        seq = self.fetcher.get_seq_from_to(chr, start, end, strand, fLeaveOpen=True)
        return seq.upper()


    def toFasta(self):
        "Fasta formatted string representation of NtSequence object --> str"

        # header: >identifier[ description]
        string = '>%s%s\n' % ( (not self.name is None and self.name or "unknown"),   \
                               (not self.desc is None and (" " + self.desc) or "") )

        for i in xrange(0, len(self.seq), self.width):
            string += "%s\n" % self.seq[i:i+self.width]

        return string


    def __repr__(self):
        "standard representation of NtSequence object --> str"
        return self.toFasta()

    def __str__(self):
        "standard string representation of NtSequence object --> str"
        return self.toFasta()

    def __len__(self):
        "dummy method necessary to trick python not to correct negative [i:j] slice indices --> 0"
        return 0

    def len(self):
        "length of sequence --> int"
        return len(self.seq)


    def __getslice__(self, i, j):
        """get a subsequence of the current sequence:

        seq[i:j]

        i,j are biological coordinates (1-based, inclusive), with no base at position zero (...-2,-1,1,2,...)
        i,j < 0        : indicate positions abs(i),abs(j) bases upstream of sequence start
                         e.g. exon[-2:-1] returns the two bases immediately upstream of the start ('AG' for canonical splice sites)
        i,j > len(seq) : indicate positions i-len(seq),j-len(seq) bases downstream of sequence end
        
        --> str"""

        tmpseq   = self.seq

        # make sure that i<=j and i,j!=0
        if i > j:
            i, j = j, i
        if i == 0:
            i = 1
        if j == 0:
            j = -1

        #print "getting slice from %i --> %i" % (i,j)

        # fetch flanking sequence
        if i < 0 or j < 0 or i > len(self.seq) or j > len(self.seq):
            tmpstart = self.start
            tmpend   = self.end
  
            if self.strand == '+':
                if i < 0:
                    tmpstart = self.start + i
                    if j < 0:
                        tmpend = self.start + j
                        j = tmpend - tmpstart + 1
                    else:
                        j = j - i
                    i = 1
                if j > len(self.seq):
                    tmpend = self.start + j - 1

            elif self.strand == '-':
                if i < 0:
                    tmpend = self.end - i
                    if j < 0:
                        tmpstart = self.end - j
                        j = tmpend - tmpstart + 1
                    else:
                        j = j - i
                    i = 1
                if j > len(self.seq):
                    tmpstart = self.end - j + 1

            #print "adjusted (%i,%i) to (%i,%i), returning range %i --> %i" % (self.start, self.end, tmpstart, tmpend, i-1, j)
            tmpseq = self.__get_seq(self.chr, tmpstart, tmpend, self.strand)

        # convert biological coordinates to array slice (0-based, exclusive) and return seq string
        return tmpseq[i-1:j]


class Exon(NtSequence):
    "Exon: represents an exon (a pair of SpliceGraph elements) and associated information"

    def __init__(self, el1, el2, startPhase=None, endPhase=None, name=None, desc=None, species=None, assembly=None):
        "instantiate Exon --> Exon"
        if name is None:
            name = '%s:E:%s' % (el1, el2)
        NtSequence.__init__(self, el1, el2, name, desc, species, assembly)
        if startPhase is None:
            self.startPhase = -1
        else:
            self.startPhase = startPhase
        if endPhase is None:
            if self.startPhase == -1:
                self.endPhase = -1
            else:
                self.endPhase = (self.startPhase + (self.end - self.start + 1)) % 3
        else:
            self.endPhase = endPhase
        if self.desc is None:
            self.desc = "%i,%i" % (self.startPhase, self.endPhase)


class Intron(NtSequence):
    "Intron: represents an intron (a pair of SpliceGraph elements) and associated information"

    def __init__(self, el1, el2, phase=-1, name=None, desc=None, species=None, assembly=None):
        "instantiate Intron --> Intron"
        if name is None:
            name = '%s:I:%s' % (el1, el2)
        # correct el1,el2 for intron coordinates
        el1I, el2I = (None, None)
        el1parts = el1.split(':')
        el2parts = el2.split(':')
        if el1parts[2] == '+':
            el1I = ':'.join([el1parts[0], str(int(el1parts[1])+1), el1parts[2]])
            el2I = ':'.join([el2parts[0], str(int(el2parts[1])-1), el2parts[2]])
        else:
            el1I = ':'.join([el1parts[0], str(int(el1parts[1])-1), el1parts[2]])
            el2I = ':'.join([el2parts[0], str(int(el2parts[1])+1), el2parts[2]])
        NtSequence.__init__(self, el1I, el2I, name, desc, species, assembly)
        self.phase = phase
        if self.desc is None:
            self.desc = "%i" % (self.phase)

    def isCanonical(self):
        "test if splice sites are canonical (GT->AG or GC->AG) --> bool"
        return ( (self.seq.startswith('GT') or self.seq.startswith('GC')) and self.seq.endswith('AG') )


class SS(NtSequence):
    "SS: represents a splice site (a single SpliceGraph element) and associated information"

    SSlen = {'5ss':8, '3ss':20} # number of nucleotides to represent the splice site

    def __init__(self, el, type, name=None, desc=None, species=None, assembly=None):
        "instantiate SS --> SS"

        if not type in ['5ss','3ss']:
            raise Errors.ObjectInitError('SS.__init__', "'%s' is not a a valid splice site type\n" % type)
        
        if name is None:
            name = el
        # correct el for intron coordinates
        el1, el2 = (None, None)
        elparts  = el.split(':')
        if elparts[2] == '+':
            if type == '5ss':
                el1 = ':'.join([elparts[0], str(int(elparts[1])+1),                 elparts[2]])
                el2 = ':'.join([elparts[0], str(int(elparts[1])+self.SSlen['5ss']), elparts[2]])
            elif type == '3ss':
                el1 = ':'.join([elparts[0], str(int(elparts[1])-self.SSlen['3ss']), elparts[2]])
                el2 = ':'.join([elparts[0], str(int(elparts[1])-1),                 elparts[2]])
        else:
            if type == '5ss':
                el1 = ':'.join([elparts[0], str(int(elparts[1])-1),                 elparts[2]])
                el2 = ':'.join([elparts[0], str(int(elparts[1])-self.SSlen['5ss']), elparts[2]])
            elif type == '3ss':
                el1 = ':'.join([elparts[0], str(int(elparts[1])+self.SSlen['3ss']), elparts[2]])
                el2 = ':'.join([elparts[0], str(int(elparts[1])+1),                 elparts[2]])
        NtSequence.__init__(self, el1, el2, name, desc, species, assembly)
        self.type = type
        if self.desc is None:
            self.desc = self.type

    def isCanonical(self):
        "test if splice sites are canonical (GT->AG or GC->AG) --> bool"
        if self.type == '5ss':
            return self.seq.startswith('GT') or self.seq.startswith('GC')
        elif self.type == '3ss':
            return self.seq.endswith('AG')


class SpliceGraphIterator:
    "SpliceGraphIterator: Iterates over SpliceGraph object in a directory"

    fileExtension = '.sg'

    def __init__(self, dir, doChecks=True, verbose=False, species="unknown", assembly="unknown", conf=None, log=sys.stderr):
        "instantiate SpliceGraphIterator --> SpliceGraphIterator"

        # class attributes
        self.dir      = dir       # directory containing splice graph flat files
        self.verbose  = verbose   # reporting level
        self.files    = None      # splice graph flat files in self.dir
        self.number   = 0         # set to len(self.files)
        self.index    = None      # points to next file to be read by self.next_sg()
        self.log      = log       # passed to SpliceGraph.__init__
        self.doChecks = doChecks  # passed to SpliceGraph.__init__
        self.species  = species   # passed to SpliceGraph.__init__
        self.assembly = assembly  # passed to SpliceGraph.__init__
        self.conf     = None      # passed to SpliceGraph.__init__

        # read configuration
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = configuration.Configuration()
        
        # check dir
        if not os.path.isdir(self.dir):
            self.__inform('SpliceGraphIterator: %s is not a directory\n' % self.dir)
            raise Errors.ObjectInitError('SpliceGraphIterator.__init__', "'%s' is not a directory\n" % self.dir)

        # get input files
        self.files = [ f for f in os.listdir(self.dir) if f.endswith(self.fileExtension) ]
        self.number = len(self.files)
        if self.number > 0:
            self.files.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]),int((y.split(':')[1]))))
            self.index = 0


    def current_index(self):
        "Return the index of the last returned splice graph --> integer|None"
        if not self.index is None:
            return (self.index - 1)
        else:
            return None


    def next_sg(self):
        "Return the next splice graph in the directory --> SpliceGraph, False if there are no (more) files"
        if (not self.index is None) and (self.index < len(self.files)):
            #print self.index, ' of ', len(self.files), '  ', self.files[self.index]
            sg = False
            while (sg == False and self.index < len(self.files)):
                try:
                    sg = SpliceGraph(filename=self.dir+'/'+self.files[self.index], log=self.log, conf=self.conf, \
                                     doChecks=self.doChecks, species=self.species, assembly=self.assembly)
                    self.index += 1
                except:
                    self.__inform('SpliceGraphIterator: ERROR creating SpiceGraph object for %s: %s\n' % (self.files[self.index], sys.exc_info()[1]))
                    sg = False
                else:
                    break
            return sg
        else:
            return False


    def reset(self):
        "Reset iteration and start over"
        if len(self.files) > 0:
            self.index = 0


    def __inform(self, txt, level=0):
        "output progress information to self.log if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            self.log.write(txt)

        
class SpliceGraphIteratorChromosomes:
    "SpliceGraphIteratorChromosomes: Iterates over SpliceGraph object in chromosme directories of a species"

    def __init__(self, dir, prefix='SpliceGraphs_', doChecks=True, verbose=False, species="unknown", assembly="unknown", conf=None, log=sys.stderr):
        "instantiate SpliceGraphIteratorChromosomes --> SpliceGraphIteratorChromosomes"

        # class attributes
        self.topdir    = dir             # iterate over all SpliceGraphs* directories in this dir
        self.dirs      = []              # list of all subdirs with splice graphs
        self.dirnumber = 0               # number of subdirs to iterate over
        self.dirindex  = None            # index of current subdir in self.dirs
        self.diriter   = None            # current SpliceGraphIterator instance (iterates over self.dirs[self.dirindex])
        self.prefix    = prefix          # directory prefix over which to iterate
        self.verbose   = verbose
        self.log       = log

        self.doChecks  = doChecks        # to be passed to SpliceGraphIterator
        self.species   = species         # to be passed to SpliceGraphIterator
        self.assembly  = assembly        # to be passed to SpliceGraphIterator

         # read configuration
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = configuration.Configuration()
        
        # check dir
        if not os.path.isdir(self.topdir):
            self.__inform('SpliceGraphIteratorChromosomes: %s is not a directory\n' % self.topdir)
            raise Errors.ObjectInitError('SpliceGraphIteratorChromosomes.__init__', "'%s' is not a directory\n" % self.topdir)

       # get input directries and self.diriter
        self.dirs = [ "%s/%s" % (self.topdir, d) for d in os.listdir(self.topdir) if d.startswith(self.prefix) and os.path.isdir(self.topdir+'/'+d) ]
        if len(self.dirs) > 0:
            self.dirs.sort()
            self.dirnumber = len(self.dirs)
            self.dirindex  = 0

            # open a SpliceGraphIterator on the first directory (will be stored in self.diriter)
            self.next_dir()


    def next_dir(self):
        "Generate a SpliceGraphIterator() instance on the next directory and store it in self.diriter --> string"
        if self.dirindex < self.dirnumber:
            self.diriter = SpliceGraphIterator(self.dirs[self.dirindex], \
                                               doChecks=self.doChecks, verbose=self.verbose, \
                                               species=self.species, assembly=self.assembly, \
                                               conf=self.conf, log=self.log)
            self.dirindex += 1
            # make sure there are *.sg files in the new directory
            while self.dirindex < self.dirnumber and self.diriter.number == 0:
                self.diriter = SpliceGraphIterator(self.dirs[self.dirindex], \
                                                   doChecks=self.doChecks, verbose=self.verbose, \
                                                   species=self.species, assembly=self.assembly, \
                                                   conf=self.conf, log=self.log)
                self.dirindex += 1

            if self.diriter.number == 0:
                self.diriter = None
                self.dirindex = None

        else:
            self.diriter = None
            self.dirindex = None

        return self.current_dir()


    def number_dirs(self):
        "Return the total number of chromosome directories in current self.topdir --> integer"
        return self.dirnumber


    def current_dir(self):
        "Return the currently opened chromosome directory --> string|None"
        if not self.dirindex is None and (self.dirindex-1) < self.dirnumber:
            return self.dirs[self.dirindex-1]
        else:
            return None


    def current_dirindex(self):
        "Return the index of the next to be opened chromosome directory --> integer|None"
        if not self.dirindex is None:
            return (self.dirindex - 1)
        else:
            return None


    def number(self):
        "Access to self.diriter attribute --> integer"
        return self.diriter.number


    def current_index(self):
        "Access to self.diriter methods --> integer|None"
        return self.diriter.current_index


    def next_sg(self):
        "Return the next splice graph in the current sub-directory --> SpliceGraph, False if there are no (more) files"

        # check if there are any files left
        if isinstance(self.diriter, SpliceGraphIterator) and \
           ( self.diriter.current_index() is None or self.diriter.number == (self.diriter.current_index() + 1) ):
            # time to try next subdir
            self.next_dir()

        if isinstance(self.diriter, SpliceGraphIterator):
            return self.diriter.next_sg()
        else: # no SpliceGraphIterator (maybe no more directories)
            return False


    def __inform(self, txt, level=0):
        "output progress information to self.log if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            self.log.write(txt)



if __name__ == "__main__":
    import SpliceGraph
    print SpliceGraph.__doc__

##     sg = SpliceGraph.SpliceGraph(connections={'chrX:100:+':{'chrX:200:+':2},
##                                               'chrX:200:+':{'chrX:1000:+':-1,'chrX:2100:+':-1},
##                                               'chrX:1000:+':{'chrX:1100:+':1},
##                                               'chrX:1100:+':{'chrX:2100:+':-1},
##                                               'chrX:2100:+':{'chrX:2200:+':2}},
##                                  types={'TSS':{'chrX:100:+':2},
##                                         '3ss':{'chrX:1000:+':1,'chrX:2100:+':2},
##                                         '5ss':{'chrX:200:+':2,'chrX:1100:+':1},
##                                         'TER':{'chrX:2200:+':2}},
##                                  transcripts={'chrX:100:+':{'chrX:200:+':['tx1','tx2']},
##                                               'chrX:200:+':{'chrX:1000:+':['tx1'],'chrX:2100:+':['tx2']},
##                                               'chrX:1000:+':{'chrX:1100:+':['tx1']},
##                                               'chrX:1100:+':{'chrX:2100:+':['tx1']},
##                                               'chrX:2100:+':{'chrX:2200:+':['tx1','tx2']}},
##                                  species="hopschel",
##                                  assembly="fischpa1",
##                                  )
        

##     # testing read/store/visualize
##     sg.store('test.sg')
##     sg.visualize('test.ps')

##     sg.read('test.sg')
##     sg.visualize('test_reread.ps')

##     sg2 = SpliceGraph.SpliceGraph(filename='test.sg')
##     sg2.visualize('test2.ps', {'chrX:1000:+':{'chrX:1100:+':'SE'}})

##     sg3 = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:10120851:-.sg')

##     #sg4 = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:42285183:-.sg')
##     #print >> sys.stderr, sg

##     # testing areConnected/areConnectedNotVia
##     print 'testing simple paths...'
##     for (el1, el2, exp) in (('chr21:10120796:-','chr21:10080032:-','yes'), \
##                             ('chr21:10080032:-','chr21:10120796:-','no'), \
##                             ('chr21:10120796:-','chr21:10048949:-','no') ):
##         if sg3.areConnected(el1, el2):
##             print "%s and %s are connected (expected: %s)" % (el1, el2, exp)
##         else:
##             print "%s and %s are NOT connected (expected: %s)" % (el1, el2, exp)
##     print

##     print "testing 'not-via' paths..."
##     for (el1, el2, vel1, vel2, exp) in (('chr21:10119414:-','chr21:10080194:-','chr21:10111624:-','chr21:10111613:-','yes'), \
##                                         ('chr21:10119414:-','chr21:10080194:-','chr21:10081688:-','chr21:10081597:-','no') ):
##         if sg3.areConnectedNotVia(el1, el2, vel1, vel2):
##             print "%s and %s are still connected (expected: %s)" % (el1, el2, exp)
##         else:
##             print "%s and %s are NOT ANY MORE connected (expected: %s)" % (el1, el2, exp)
##     print

##     # testing remove_by_coverage_and_ss/remove_by_library_keyword

##     # testing connectedComponents
##     print 'testing connected Components...'
##     sg_iter = SpliceGraph.SpliceGraphIterator('/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21')
##     raw = {}
##     sg = sg_iter.next_sg()
##     while ( sg ):
##         nbcomp = len(sg.connectedComponents())
##         if raw.has_key(nbcomp):
##             raw[nbcomp] += 1
##         else:
##             raw[nbcomp] = 1
##         sg = sg_iter.next_sg()
##     sortedkeys = raw.keys()
##     sortedkeys.sort()
##     print 'Connected components in SpliceGraphs_chr21:'
##     print '\n'.join(["\t% 4i : % 6i (% 5.1f%%)" % (k,raw[k],100*raw[k]/sg_iter.number()) for k in sortedkeys])+'\n'

##     sg_iter = SpliceGraph.SpliceGraphIterator('/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr21')
##     raw = {}
##     sg = sg_iter.next_sg()
##     while ( sg ):
##         nbcomp = len(sg.connectedComponents())
##         if raw.has_key(nbcomp):
##             raw[nbcomp] += 1
##         else:
##             raw[nbcomp] = 1
##         sg = sg_iter.next_sg()
##     sortedkeys = raw.keys()
##     sortedkeys.sort()
##     print 'Connected components in SpliceGraphsFiltered_chr21:'
##     print '\n'.join(["\t% 4i : % 6i (% 5.1f%%)" % (k,raw[k],100*raw[k]/sg_iter.number()) for k in sortedkeys])+'\n'

##     test_remove_edges()
##     sys.exit(0)

##     # testing SpliceGraph.NtSequence
##     ex = SpliceGraph.Exon('chr1:100030138:+','chr1:100030280:+')
##     ex2 = SpliceGraph.Exon('chr1:100030280:-','chr1:100030138:-')
##     print "all following pairs of sequences should be reverse complement of each other:"
##     print ex[-3:-1], '\n', ex2[ex2.len()+1:ex2.len()+3], '\n'
##     print ex[-10:4], '\n', ex2[ex2.len()-3:ex2.len()+10], '\n'
##     print ex[1:4], '\n', ex2[ex2.len()-3:ex2.len()], '\n'
##     print ex[1:ex.len()], '\n', ex2[1:ex2.len()], '\n'
##     print ex[ex.len()-3:ex.len()+3], '\n', ex2[-3:4], '\n'

##     # testing method allExonSeqs() and Exon(NtSequence)
##     sg = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:10012791:-.sg')
##     for ex in sg.allExonSeqs():
##         print ex
##     print

##     # testing SS(NtSequence)
##     sg = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr1/chr1:116627878:+.sg')
##     for el in sg.allElements():
##         if sg.getType(el) in ['5ss','3ss']:
##             elSS = SpliceGraph.SS(el, sg.getType(el))
##             print "%s : %s, %scanonical (%s)" % (el, elSS.type, (elSS.isCanonical() and "is " or "is NON-"), elSS.seq)

    # testing heaviest bundle algorithm
    sg = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:10012791:-.sg')
    #C  = sg.HB_generate_consensi()
    #for i in xrange(len(C)):
    #    print ">%s %s\n%s\n" % (C[i]['id'], C[i]['desc'], C[i]['seq'])
    #print
    #sys.exit(0)
    I = sg.HB_get_protein_isoforms()
    for i in xrange(len(I)):
        print ">%s\n%s\n" % (I[i]['id'], I[i]['seq'])
        #print I[i]
    print
    sys.exit(0)

    # testing SpliceGraphIteratorChromosomes
#    iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes('/r100/burge/shared/splice_graphs/hg17/FlatFiles', prefix='SpliceGraphsFiltered_')
#    print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)
    #current_dir = iterChrom.current_dir()
    #nb = 0
    #sg = iterChrom.next_sg()
    #while sg:
    #    nb += 1
    #    sg = iterChrom.next_sg()
    #    if sg != False and not iterChrom.current_dir() is None and current_dir != iterChrom.current_dir():
    #        print "%d splice graphs in %s" % (nb, current_dir)
    #        current_dir = iterChrom.current_dir()
    #        nb = 0
## for  i in xrange(iterChrom.number_dirs()):
##         print "%d splice graphs in %s" % (iterChrom.number(), iterChrom.current_dir())
##         iterChrom.next_dir()                       
##     iterChrom = SpliceGraph.SpliceGraphIteratorChromosomes('/r100/burge/shared/splice_graphs/hg17/FlatFiles', prefix='SpliceGraphs_')
##     #print "there are %d chromosome directories in %s\n" % (iterChrom.number_dirs(), iterChrom.topdir)
    #for i in xrange(iterChrom.number_dirs()):
    #    print "%d splice graphs in %s" % (iterChrom.number(), iterChrom.current_dir())
    #    iterChrom.next_dir()                       

