"""
 SpliceSites.py


 Class used to score splice sites according to different scoring schemes
 Currently suported classification:
    - Canonical Splice Sites (GT->AG or GC->AG)
    - ...
    
    - U12 scoring based on Chris U12Scan

 it fetches sequences using the GenomeFetch.py Module

"""
# std lib
import sys
import os

# project related modules
import GenomeFetch
import SpliceGraph


class SpliceSite:

    sequenceLength = 8   # should be dividable by 2 in order for the script to automatically find exon intron boundaries
    sequence = None  # 
    canonical = None # True/False
    defaultstrand = 1
    fLeaveOpen = False

    
    def __init__(self, species=None, assembly=None, node=None, type=None, splice_graph=None, printit=False):
        # if no species, assembly information is supplied: the GenomeFecth default will be used (hsa,hg17)
        self.species = species
        self.assembly = assembly

        # If a node was supplied to the function
        if not node is None and not type is None:
            self.get_splice_site_for_node(node=node, type=type, printit=printit)

        elif not splice_graph is None:
            self.get_splice_sites_for_splice_graph(sg=splice_graph)

    def print_splice_site(self):
        print '%s|%s' % (self.sequence[:self.sequenceLength/2],self.sequence[self.sequenceLength/2:]),
        print self.type, self.chr, self.position, self.strand, self.canonical


    def isSScanonical(self):
        return self.__annotate_canonical()

    def isU2(self):
        print 'isU2: Not implemented'

    def isU12(self):
        print 'isU12: Not implemented'

    def getSpliceSiteScoreU2(self):
        print 'isU2: Not implemented'

    def getSpliceSiteScoreU12(self):
        print 'isU12: Not implemented'

    def get_splice_site_for_node(self, node=None, type=None, printit=False):
        "retrieves splice site sequence for a node (chr:pos:strand) -> Splice Site Sequence (string)"
        parts = node.split(':')
        self.chr = parts[0]
        self.position = int(parts[1])
        if parts[2] == '+':
            self.strand = 1
        else:
            self.strand = -1
        self.type = type
        self.fetch_splice_site_sequence()
        if printit:
            self.print_splice_site()
        return self.sequence

    def get_splice_sites_for_splice_graph(self, sg, printit=False):
        "Returns a dictionary with nodes [5'ss, 3'ss] as keys and splice site sequences as entries"
        GeneomeFetchInstance = GenomeFetch.GenomeFetch(self.species, self.assembly)
        self.fLeaveOpen = True
        types = sg.allTypes()
        self.splicesites = {}
        self.splicesites['5ss'] = {}
        self.splicesites['3ss'] = {}
        for ss in types['5ss']:
            self.splicesites['5ss'][ss] = self.get_splice_site_for_node(ss,'5ss', printit=printit)           
        for ss in types['3ss']:
            self.splicesites['3ss'][ss] = self.get_splice_site_for_node(ss,'3ss', printit=printit)
        return self.splicesites

    def isSGcanonical(self):
        "Checks all splice sites in self.splicesites if they are canonical, returns a dictionary node: canonical (bool)"
        self.splicesites_canonical = {}
        for type in self.splicesites.keys():
            for ss in self.splicesites[type].keys():
                self.sequence = self.splicesites[type][ss]
                self.splicesites_canonical[ss] = self.__annotate_canonical()
        return self.splicesites_canonical
                


    def __annotate_canonical(self):
        if self.type == '5ss':
            ss = self.sequence[self.sequenceLength/2:self.sequenceLength/2+2]
            if ss == 'GT' or ss == 'GC':
                self.canonical = True
            else:
                self.canonical = False
        elif self.type == '3ss':
            ss = self.sequence[self.sequenceLength/2-2:self.sequenceLength/2]
            if ss == 'AG':
                self.canonical = True
            else:
                self.canonical = False
        return self.canonical

    def fetch_splice_site_sequence(self):
        GeneomeFetchInstance = GenomeFetch.GenomeFetch(self.species, self.assembly)

        if self.type == '5ss':
            if self.strand == 1:
                self.sequence = GeneomeFetchInstance.get_seq_from_to(self.chr,self.position-self.sequenceLength/2+1,
                                                                    self.position+self.sequenceLength/2,self.strand,
                                                                    fLeaveOpen=self.fLeaveOpen)
            else:
                self.sequence = GeneomeFetchInstance.get_seq_from_to(self.chr,self.position-self.sequenceLength/2,
                                                                    self.position+self.sequenceLength/2-1,self.strand,
                                                                    fLeaveOpen=self.fLeaveOpen)                
        elif self.type == '3ss':
            if self.strand == 1:
                self.sequence = GeneomeFetchInstance.get_seq_from_to(self.chr,self.position-self.sequenceLength/2,
                                                                    self.position+self.sequenceLength/2-1,self.strand,
                                                                    fLeaveOpen=self.fLeaveOpen)
            else:
                self.sequence = GeneomeFetchInstance.get_seq_from_to(self.chr,self.position-self.sequenceLength/2+1,
                                                                    self.position+self.sequenceLength/2,self.strand,
                                                                    fLeaveOpen=self.fLeaveOpen)
    
    def profile_splice_sites_in_splice_graph_directory(self, directory):
        'goes through all SpliceGraphs in a directory and gives statistics on splice sites from all elements'
        if not os.path.isdir(directory):
            print sys.stderr.write('WARNING: Not a directory input for SpliceSites.profile_splice_sites_in_splice_graph_directory(): %s' % directory)
        
        # get file list in directory
        files = os.listdir(directory)
        
        NrElements = {}
        
        for file in files:
            if file.split('.')[-1] != 'sg':
                # skip all other files in directory
                continue
            # read Splice Graph
            filename = '%s%s' % (directory,file)
            sg = SpliceGraph.SpliceGraph(filename=filename)
            types = sg.allTypes()

            # check splice sites, for canonical splice sites
            if not types.has_key('5ss') or not types.has_key('3ss'):
                continue
            ss = self.get_splice_sites_for_splice_graph(sg)
            ss = self.isSGcanonical()

            for key in ss:
                if types['5ss'].has_key(key):
                    if not NrElements.has_key(types['5ss'][key]):
                        NrElements[types['5ss'][key]] = {}
                        NrElements[types['5ss'][key]]['True'] = 0
                        NrElements[types['5ss'][key]]['False'] = 0
                    if ss[key] is True:
                        NrElements[types['5ss'][key]]['True'] += 1
                    else:
                        NrElements[types['5ss'][key]]['False'] += 1
                else:

                    if not NrElements.has_key(types['3ss'][key]):
                        NrElements[types['3ss'][key]] = {}
                        NrElements[types['3ss'][key]]['True'] = 0
                        NrElements[types['3ss'][key]]['False'] = 0
                    if ss[key] is True:
                        NrElements[types['3ss'][key]]['True'] += 1
                    else:
                        NrElements[types['3ss'][key]]['False'] += 1
            
            keys = NrElements.keys()
            keys.sort()

        for key in keys:
            if int(key) > 50:
                break
            print key, NrElements[key]['True'],NrElements[key]['False'] 
            
        

def testModule():
    # test splice site module
    # set printit to True to print more information on the splice site, including the sequences
    
    css = SpliceSite(node='chr21:9943805:+', type='3ss')
    print 'chr21:9943805:+ Canonical:', css.isSScanonical()
    css = SpliceSite(node='chr21:46593354:+', type='3ss')
    print 'chr21:46593354:+ Canonical:', css.isSScanonical()

    # test for all splice sites in a SpliceGraph
    print 'Testing a whole splice graph:'
    file = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:46568483:+.sg'
    sg = SpliceGraph.SpliceGraph(filename=file)
    css = SpliceSite(splice_graph=sg)
    res = css.isSGcanonical()
    Can = 0
    Total = len(res.keys())
    print file
    for ss in res.keys():
        #print ss, res[ss]
        if res[ss] is True:
            Can += 1
    print '\t%.1f%% canonical splice sites (%i of %i)' % (Can/float(Total)*100,Can,Total)
    
    
    file = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/chr21:45249070:-.sg'
    sg = SpliceGraph.SpliceGraph(filename=file)
    css = SpliceSite(splice_graph=sg)
    res = css.isSGcanonical()
    Can = 0
    Total = len(res.keys())
    print file
    for ss in res.keys():
        #print ss, res[ss]
        if res[ss] is True:
            Can += 1
    print '\t%.1f%% canonical splice sites (%i of %i)' % (Can/float(Total)*100,Can,Total)
    
    # check all Splice Graphs ina directory
    print 'Testing all Splice Graphs in a directory'
    css = SpliceSite()
    css.profile_splice_sites_in_splice_graph_directory('/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/')
    

def printhelp():
    print """\nSpliceSites Class

    USAGE [OPTIONS]

        -G   splice_graph file
        -N   a node 'chr:pos:strand (also need to submit type)
        -T   type
        -S   species (default hsa)
        -A   assembly (defualt hg17)
        -t   run testModule script
        """
    


from optparse import OptionParser

if __name__ == "__main__":
    opts = OptionParser()
    opts.add_option('-G','--splicegraph', dest='splice_graph_file')
    opts.add_option('-N','--node', dest='node')
    opts.add_option('-T','--type',dest='type')
    opts.add_option('-S','--species',dest='species')
    opts.add_option('-A','--assembly',dest='assembly')
    opts.add_option('-t','--Moduletest',action='store_true', default=False, dest='Test')

    (options, args) = opts.parse_args()
    
    if options.splice_graph_file:
        sg = SpliceGraph.SpliceGraph(filename=options.splice_graph_file)
        css = SpliceSite(splice_graph=sg)
        ss_canonical = css.isSGcanonical()

        Can = 0
        Total = len(ss_canonical.keys())
        for site in ss_canonical.keys():
            print site, ss_canonical[site]
            if ss_canonical[site] is True:
                Can += 1
        print '%f%% canonical splice sites (%i of %i)' % (Can/float(Total)*100,Can,Total)
            
    elif options.node and options.type:
        css = SpliceSite(species=options.species,assembly=options.assembly, node=options.node, type=options.type)
        print options.node, 'Canonical:', css.isSScanonical()
        
    elif options.Test is True:
        testModule()

    else:
        printhelp()
    
