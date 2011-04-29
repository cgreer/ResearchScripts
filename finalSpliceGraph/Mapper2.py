#---------------------------------------------------------------------------------------
#   $Id: Mapper.py,v 1.21 2005/11/02 14:35:48 stadler Exp $    
#---------------------------------------------------------------------------------------
"""Mapper: Class for mapping splice graphs by chromosomal coordinates to known genes

The mapping information is parsed from annotation files that are stored in the MetaData
directory. The file names of these annotation files is as follows:

    species_assebmly_type.gz

e.g.:
    hsa_hg17_ensembl.gz

The class will parse annotation files of the following types:
    ensembl : tab separated list, generated in EnsMART (Structures)
              with the following columns:
              - Chromosome
              - Strand
              - Ensembl Gene ID
              - Ensembl Transcript ID
              - External Gene ID
              - Biotype
              - Ensembl Exon ID (versioned)
              - Exon Start (Chr bp)
              - Exon End (Chr bp)
             

"""

# standard lib
import os
import sys
import re
import gzip

# project related modules
import Errors
import SpliceGraph
import configuration


class Mapper2:
    "Mapper: map a SpliceGraph to a known gene by chromosomal coordinates"

    __defaultSpecies    = 'mmu'
    __defaultAnnotation = 'gnf1m'

    def __init__(self, configFile=None, species=None, parseAnnotation=True, verbose=False, log=sys.stderr):
        "instantiate Mapper object --> Mapper"

        # class attributes
        self.verbose         = verbose
        self.log             = log
        self.species         = ( (species is None) and self.__defaultSpecies or species )
        self.assembly        = None
        self.annotationFiles = {}      #format: self.annotationFiles[type] = file_with_complete_path
        self.exonInfo        = {}      #format: self.exonInfo[type][chr][strand][start][end] = (gene_id, exon_id) #with start<end
        self.ssInfo          = {}      #format: self.ssInfo[type][chr][strand][coord] = (gene_id, exon_id)
        self.geneInfo        = {}      #format: self.geneInfo[type][gene_id][info] = value

        # init configuration
        try:
            self.conf = configuration.Configuration(configFile)
        except:
            raise Errors.ObjectInitError('Mapper.__init__', "could not get configuration from '%s'" % configFile)
        else:
            # check species/assembly
            if self.species not in self.conf.getSpeciesList():
                raise Errors.ObjectInitError('Mapper.__init__', "species %s not found in config file '%s'" % (self.species, self.conf.configFile()) )
            else:
                self.assembly = self.conf[self.species]['Assembly']

            # get available annotation files
            annotFiles = '/r100/burge/shared/splice_graphs/mm6/MetaData/mmu_mm6_gnf1m.gz'
            #[ f for f in os.listdir(self.conf[self.species]['MetaDataStoragePath']) if re.search('^[^_]+_[^_]+_[^_]+\.gz$', f) ]
            for annotFile in annotFiles:
                (spec, ass, type) = annotFile.split('_')
                type = type[0:len(type)-3]
                if spec == self.species and ass == self.assembly:
                    self.annotationFiles[type] = "%s/%s" % (self.conf[self.species]['MetaDataStoragePath'],annotFile)


    def parseAnnotation(self, type):
        "parse and store annotation files from MetaData"

        if type not in self.annotationFiles.keys():
            raise Errors.ArgumentError('Mapper.parseAnnotation', 'no annotation file for type "%s"' % type)

        elif type == 'ensembl' or type == 'vega' or type == 'genrate' or type == "HuExonStv2":
            versionPattern = re.compile('\.\d+$')
            self.geneInfo[type]   = {}
            self.ssInfo[type]     = {}
            self.exonInfo[type]   = {}
            nb = 0

            fh = gzip.GzipFile(self.annotationFiles[type], 'r')

            while 1:
                line = fh.readline()
                if not line:
                    break
                elif line.startswith('#') or line.startswith('Chromosome'):
                    continue
                else:
                    line = line.strip('\n')
                    (chr, strand, gnId, txId, extGnId, biotype, exId, start, end) = line.split('\t')
                    if not chr.startswith('chr'):
                        chr = "chr%s" % chr
                    strand = ((strand == '1') and '+' or '-')
                    exId = versionPattern.sub('', exId)

                    # store in self.geneInfo, self.ssInfo and self.exonInfo
                    if not self.geneInfo[type].has_key(gnId):
                        self.geneInfo[type][gnId] = {'extGnId':extGnId,'biotype':biotype}

                    if not self.ssInfo[type].has_key(chr):
                        self.ssInfo[type][chr] =  {'+':{}, '-':{}}
                    self.ssInfo[type][chr][strand][start] = (gnId, exId)
                    self.ssInfo[type][chr][strand][end]   = (gnId, exId)

                    if not self.exonInfo[type].has_key(chr):
                        self.exonInfo[type][chr] = {'+':{}, '-':{}}
                    if not self.exonInfo[type][chr][strand].has_key(start):
                        self.exonInfo[type][chr][strand][start] = {}
                    self.exonInfo[type][chr][strand][start][end] = (gnId, exId)

                    nb += 1
                    self.__inform("just stored %s -> %s : %s\n" % (start, end, gnId), 5)
                    self.__inform("\t%s,%s,%s,%s\n" % (type, chr, strand,start), 7)

            fh.close()
            self.__inform("parseAnnotation: successfully parsed gene association of %i %s exons\n" % (nb,type))

        elif type == 'U133target' or type == 'gnf1m':
            self.geneInfo[type]   = {}
            self.ssInfo[type]     = {}
            self.exonInfo[type]   = {}
            nb = 0

            fh = gzip.GzipFile(self.annotationFiles[type], 'r')

            while 1:
                line = fh.readline()
                if not line:
                    break
                elif line.startswith('#'):
                    continue
                else:
                    line = line.strip('\n')
                    (probesetId, chr, start, end, strand) = line.split('\t')

                    # store in self.exonInfo
                    if not self.exonInfo[type].has_key(chr):
                        self.exonInfo[type][chr] = {'+':{}, '-':{}}
                    if not self.exonInfo[type][chr][strand].has_key(start):
                        self.exonInfo[type][chr][strand][start] = {}
                    self.exonInfo[type][chr][strand][start][end] = (probesetId, probesetId)

                    nb += 1
                    self.__inform("just stored %s -> %s : %s\n" % (start, end, probesetId), 5)

            fh.close()
            self.__inform("parseAnnotation: successfully parsed gene association of %i %s exons\n" % (nb,type))

        elif type == 'hgu95av2':
            #versionPattern = re.compile('\.\d+$')
            self.geneInfo[type]   = {}
            self.ssInfo[type]     = {}
            self.exonInfo[type]   = {}
            nb = 0

            fh = gzip.GzipFile(self.annotationFiles[type], 'r')

            while 1:
                line = fh.readline()
                if not line:
                    break
                else:
                    line = line.strip('\n')
                    
                blocks = line.split('\t')

                (probeId, map ) = blocks[0].split('::')
                (chr, start, end, strand) = map.split(':')

                if not self.ssInfo[type].has_key(chr):
                        self.ssInfo[type][chr] =  {'+':{}, '-':{}}
                self.ssInfo[type][chr][strand][start] = (probeId, probeId)
                self.ssInfo[type][chr][strand][end]   = (probeId, probeId)

                if not self.exonInfo[type].has_key(chr):
                    self.exonInfo[type][chr] = {'+':{}, '-':{}}
                if not self.exonInfo[type][chr][strand].has_key(start):
                    self.exonInfo[type][chr][strand][start] = {}
                self.exonInfo[type][chr][strand][start][end] = (probeId, probeId)

                
                nb += 1
                self.__inform("just stored %s -> %s : %s\n" % (start, end, probeId), 5)
            fh.close()
            self.__inform("parseAnnotation: successfully parsed gene association of %i %s exons\n" % (nb,type))


    def getGeneInfo(self, gnId, type, info):
        "return additional gene info from annotation, if available"

        if self.geneInfo.has_key(type) and \
           self.geneInfo[type].has_key(gnId) and \
           self.geneInfo[type][gnId].has_key(info):
            return self.geneInfo[type][gnId][info]
        else:
            return "unknown"


    def overlappingGenesExons(self, type, exStart, exEnd, detailed=0):
        """return a (gene-list, exon-dict)-tuple of all genes/exons in annotation 'type'
        that overlap the exon exStart-exEnd on any strand

        gene-list is a list of gene ids
        exon-dict is a dictionary with exon ids as keys, and corresponding gene ids as values
        if detailed flag is set to 1: the exon values will contain a 3-tuple with
        ( corresponding gene ids, exonStart exonEnd )

        --> (list,dict)"""

        genes = {}
        exons = {}

        (exChromosome, exStartCoord, exStrand) = exStart.split(':')
        exEndCoord = exEnd.split(':')[1]

        # exStart, exEnd are biological (not start<end), but we need them start<end --> swap them if exStrand=='-'
        if exStrand == '-':
            exStartCoord, exEndCoord = exEndCoord, exStartCoord

        if self.exonInfo.has_key(type) and \
               self.exonInfo[type].has_key(exChromosome) and \
               self.exonInfo[type][exChromosome].has_key(exStrand):
            annStarts = self.exonInfo[type][exChromosome][exStrand].keys()
            for annStartCoord in annStarts:
                annEnds = self.exonInfo[type][exChromosome][exStrand][annStartCoord].keys()
                for annEndCoord in annEnds:
                    if self.__overlap(exStartCoord, exEndCoord, annStartCoord, annEndCoord):
                        if genes.has_key(self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord]):
                            genes[self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]] += 1
                        else:
                            genes[self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]]  = 1
                        if detailed == 0:
                            exons[self.exonInfo[type][exChromosome][exStrand]\
                                  [annStartCoord][annEndCoord][1]] = self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]
                        elif detailed == 1:
                            exons[self.exonInfo[type][exChromosome][exStrand]\
                                  [annStartCoord][annEndCoord][1]] = (self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0],annStartCoord,annEndCoord)
        return genes.keys(), exons
    
    def withinExons(self, type, exStart, exEnd, detailed=0):
        """return a (gene-list, exon-dict)-tuple of all genes/exons in annotation 'type'
        that overlap the exon exStart-exEnd on any strand

        gene-list is a list of gene ids
        exon-dict is a dictionary with exon ids as keys, and corresponding gene ids as values
        if detailed flag is set to 1: the exon values will contain a 3-tuple with
        ( corresponding gene id, exonStart exonEnd )

        --> (list,dict)"""

        genes = {}
        exons = {}

        (exChromosome, exStartCoord, exStrand) = exStart.split(':')
        exEndCoord = exEnd.split(':')[1]

        # exStart, exEnd are biological (not start<end), but we need them start<end --> swap them if exStrand=='-'
        if exStrand == '-':
            exStartCoord, exEndCoord = exEndCoord, exStartCoord
            
        self.__inform("withinExons: analyzing %s, %s, %s, %s\n" % (exChromosome, exStrand, exStartCoord, exEndCoord),7)
 
        if self.exonInfo.has_key(type) and \
               self.exonInfo[type].has_key(exChromosome) and \
               self.exonInfo[type][exChromosome].has_key(exStrand):
            annStarts = self.exonInfo[type][exChromosome][exStrand].keys()
            for annStartCoord in annStarts:
                annEnds = self.exonInfo[type][exChromosome][exStrand][annStartCoord].keys()
                for annEndCoord in annEnds:
                    self.__inform("withinExons: chr and strand found. now checking annotations",7)
                    if self.__within(exStartCoord, exEndCoord, annStartCoord, annEndCoord):
                        if genes.has_key(self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord]):
                            genes[self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]] += 1
                        else:
                            genes[self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]]  = 1
                        if detailed == 0:
                            exons[self.exonInfo[type][exChromosome][exStrand]\
                                  [annStartCoord][annEndCoord][1]] = self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0]
                        elif detailed == 1:
                            exons[self.exonInfo[type][exChromosome][exStrand]\
                                  [annStartCoord][annEndCoord][1]] = (self.exonInfo[type][exChromosome][exStrand][annStartCoord][annEndCoord][0],annStartCoord,annEndCoord)
        return genes.keys(), exons

    def __overlap(self, start1, end1, start2, end2):
        "return true if two features overlap --> bool"
        return int(end1) >= int(start2) and int(start1) <= int(end2)
    
    def __within(self, start1, end1, start2, end2):
        "return true if feature 2 is fully included within feature 1 (the exon)"
        return int(start1) <= int(start2) and int(end1) >= int(end2)

    def map2gene(self, sg, type, mode="singleSS"):
        """map a SpliceGraph to a gene according to the specified exon matching mode

        returns a (gene-dict,exon-dict)-tuple, where
          gene-dict is a dictionary with mapped gene ids as keys and the number of matching exons as values
          exon-dict is a dictionary with (sg.name, sg.exon_ids)-tuples as keys and the corresponding (gnId, exId)-tuples as values

        --> (dict, dict)
        """

        genes = {} #format: genes[gnId] = nb_of_exons
        exons = {} #format: exons[(sg.name,sg.exon)] = (mapped.gene, mapped.exon) #where sg.exon=element1:E:element2
        
        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # get annotation
            if not self.annotationFiles.has_key(type): # give up
                self.__inform("ERROR in Mapper.map2gene: unknown annotation type '%s'" % type)
                raise Errors.ArgumentError("Mapper.map2gene", "unknown annotation type '%s'" % type)

            elif not self.exonInfo.has_key(type):    # parse it
                self.parseAnnotation(type)

            # for each exon, get associated gene
            for exStart in sg.allElements():

                if sg.is3ss(exStart) or sg.isTSS(exStart):

                    # look for a corresponding exon end (5ss element)
                    for exEnd in sg.downstreamConnectedElements(exStart):

                        # make sure exEnd is a 5ss element
                        if sg.is5ss(exEnd) or sg.isTER(exEnd):
                            
                            # remember: exStart/exEnd are biological, not always exStart<exEnd --> flip
                            (exChromosome, exStartCoord, exStrand) = exStart.split(':')
                            exEndCoord = exEnd.split(':')[1]
                            if exStrand == '-':
                                exStartCoord, exEndCoord = exEndCoord, exStartCoord

                            # if correspondance is found with annotation, store gene and exon IDs
                            if mode == "exon" and \
                                self.exonInfo[type].has_key(exChromosome) and \
                                self.exonInfo[type][exChromosome].has_key(exStrand) and \
                                self.exonInfo[type][exChromosome][exStrand].has_key(exStartCoord) and \
                                self.exonInfo[type][exChromosome][exStrand][exStartCoord].has_key(exEndCoord):
                                if genes.has_key(self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]):
                                    genes[self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]] += 1
                                else:
                                    genes[self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]]  = 1
                                exons[(sg.name, exStart + ":E:" + exEnd)] = self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord]

                            elif mode == "bothSS" and \
                                 self.ssInfo[type].has_key(exChromosome) and \
                                 self.ssInfo[type][exChromosome].has_key(exStrand) and \
                                 self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) and \
                                 self.ssInfo[type][exChromosome][exStrand].has_key(exEndCoord) and \
                                 self.ssInfo[type][exChromosome][exStrand][exStartCoord] == self.ssInfo[type][exChromosome][exStrand][exEndCoord]:
                                if genes.has_key(self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]):
                                    genes[self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]] += 1
                                else:
                                    genes[self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]]  = 1
                                exons[(sg.name, exStart + ":E:" + exEnd)] = self.ssInfo[type][exChromosome][exStrand][exStartCoord]

                            elif mode == "singleSS" and \
                                 self.ssInfo[type].has_key(exChromosome) and \
                                 self.ssInfo[type][exChromosome].has_key(exStrand) and \
                                 ( self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) or \
                                   self.ssInfo[type][exChromosome][exStrand].has_key(exEndCoord) ):
                                ss = self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) and exStartCoord or exEndCoord
                                if genes.has_key(self.ssInfo[type][exChromosome][exStrand][ss][0]):
                                    genes[self.ssInfo[type][exChromosome][exStrand][ss][0]] += 1
                                else:
                                    genes[self.ssInfo[type][exChromosome][exStrand][ss][0]]  = 1
                                exons[(sg.name, exStart + ":E:" + exEnd)] = self.ssInfo[type][exChromosome][exStrand][ss]

                            elif mode == "overlap":
                                (overlapGns, overlapExons) = self.overlappingGenesExons(type, exStart, exEnd)
                                for overlapGn in overlapGns:
                                    if genes.has_key(overlapGn):
                                        genes[overlapGn] += 1
                                    else:
                                        genes[overlapGn]  = 1
                                for overlapExon in overlapExons.iterkeys():
                                    exons[(sg.name, exStart + ":E:" + exEnd)] = (overlapExons[overlapExon], overlapExon)
                                    
                            elif mode == "within":
                                (overlapGns, overlapExons) = self.withinExons(type, exStart, exEnd)
                                for overlapGn in overlapGns:
                                    if genes.has_key(overlapGn):
                                        genes[overlapGn] += 1
                                    else:
                                        genes[overlapGn]  = 1
                                for overlapExon in overlapExons.iterkeys():
                                    exons[(sg.name, exStart + ":E:" + exEnd)] = (overlapExons[overlapExon], overlapExon)
                                
                            elif mode not in ["bothSS", "singleSS", "overlap", "exon","within"]:
                                self.__inform("map2gene: Ignoring unknown mode %s\n" % mode)

        else:
            self.__inform("ERROR in map2gene: 1st argument must either be filename or instance of SpliceGraph")
            raise Errors.ArgumentError("Mapper.map2gene", "1st argument must either be filename or instance of SpliceGraph")

        return genes, exons
    
    def mapEvents(self, event, event_type, type, outFile=None,mode="within"):
        """maps a list of events (e.g. CE or SE) in the raw text format stored in *_events files 
        to annotation type using the specifid mode, by default within (to use for microarray probeset regions"""
        
        if isinstance(event,basestring):
            if os.path.exists(event):
                f = open(event,"r")
                event = f.readlines()
                f.close()
        if not isinstance(event, list):
            raise Errors.ArgumentError("Mapper.mapEvents", "event argument must either be a filename or list of events (strings)")
            
        # get annotation
        if not self.annotationFiles.has_key(type): # give up
            self.__inform("ERROR in Mapper.map2gene: unknown annotation type '%s'" % type)
            raise Errors.ArgumentError("Mapper.map2gene", "unknown annotation type '%s'" % type)

        elif not self.exonInfo.has_key(type):    # parse it
            self.parseAnnotation(type)
            
        genes = {} #format: genes[gnId] = nb_of_exons
        exons = {} #format: exons[(sg.name,sg.exon)] = (mapped.gene, mapped.exon) #where sg.exon=element1:E:element2
        
        # open outFile fh
        if outFile is None:
            outFile = "%s.mapped_%s" % (event_type, type)
        fh = open(outFile, "w")
        
        
        if event_type == "SE" or event_type == "CE" or event_type == "CI" or event_type == "RI":
            fh.write("#gnId\texStart\texEnd\tExon_id\tExon_coordinates\n")
            for ev in event:
                if ev.startswith("#"): continue
                line = ev.strip().replace('\n','')
                if event_type == "SE":
                    gnId, exStart, exEnd, skCoverage, incCoverage = line.split('\t')
                elif event_type == "CE" or event_type == "CI":
                    gnId, exStart, exEnd, Coverage  = line.split('\t')
                elif event_type == "RI":
                    gnId, exStart, exEnd, retainedCoverage, splicedCoverage = line.split('\t')
                elif event_type == "A5SS" or event_type == "A3SS":
                    gnId, nbASS, altSS, Coverage, anotherSS, yetanotherSS = line.split('\t')
                    if nbASS == 2:
                        exStart, exEnd = altSS.split(',')
                    else:
                        continue # try next event from the list
                else:
                    sys.stderr.write("Error. Mapper.MapEvents. event not implemented: %s" % event_type)
                    
                # remember: exStart/exEnd are biological, not always exStart<exEnd --> flip
                (exChromosome, exStartCoord, exStrand) = exStart.split(':')
                exEndCoord = exEnd.split(':')[1]
                if exStrand == '-':
                    exStartCoord, exEndCoord = exEndCoord, exStartCoord
                
                # if correspondance is found with annotation, store gene and exon IDs
                if mode == "exon" and \
                    self.exonInfo[type].has_key(exChromosome) and \
                    self.exonInfo[type][exChromosome].has_key(exStrand) and \
                    self.exonInfo[type][exChromosome][exStrand].has_key(exStartCoord) and \
                    self.exonInfo[type][exChromosome][exStrand][exStartCoord].has_key(exEndCoord):
                    if genes.has_key(self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]):
                        genes[self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]] += 1
                    else:
                        genes[self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord][0]]  = 1
                    exons[(gnId, exStart + ":E:" + exEnd)] = self.exonInfo[type][exChromosome][exStrand][exStartCoord][exEndCoord]

                elif mode == "bothSS" and \
                     self.ssInfo[type].has_key(exChromosome) and \
                     self.ssInfo[type][exChromosome].has_key(exStrand) and \
                     self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) and \
                     self.ssInfo[type][exChromosome][exStrand].has_key(exEndCoord) and \
                     self.ssInfo[type][exChromosome][exStrand][exStartCoord] == self.ssInfo[type][exChromosome][exStrand][exEndCoord]:
                    if genes.has_key(self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]):
                        genes[self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]] += 1
                    else:
                        genes[self.ssInfo[type][exChromosome][exStrand][exStartCoord][0]]  = 1
                    exons[(gnId, exStart + ":E:" + exEnd)] = self.ssInfo[type][exChromosome][exStrand][exStartCoord]

                elif mode == "singleSS" and \
                     self.ssInfo[type].has_key(exChromosome) and \
                     self.ssInfo[type][exChromosome].has_key(exStrand) and \
                     ( self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) or \
                       self.ssInfo[type][exChromosome][exStrand].has_key(exEndCoord) ):
                    ss = self.ssInfo[type][exChromosome][exStrand].has_key(exStartCoord) and exStartCoord or exEndCoord
                    if genes.has_key(self.ssInfo[type][exChromosome][exStrand][ss][0]):
                        genes[self.ssInfo[type][exChromosome][exStrand][ss][0]] += 1
                    else:
                        genes[self.ssInfo[type][exChromosome][exStrand][ss][0]]  = 1
                    exons[(gnId, exStart + ":E:" + exEnd)] = self.ssInfo[type][exChromosome][exStrand][ss]

                elif mode == "overlap":
                    (overlapGns, overlapExons) = self.overlappingGenesExons(type, exStart, exEnd)
                    for overlapGn in overlapGns:
                        if genes.has_key(overlapGn):
                            genes[overlapGn] += 1
                        else:
                            genes[overlapGn]  = 1
                    for overlapExon in overlapExons.iterkeys():
                        exons[(gnId, exStart + ":E:" + exEnd)] = (overlapExons[overlapExon], overlapExon)
                        fh.write("%s\t%s\t%s\t%s\t%s\n" % (gnId, exStart, exEnd, overlapExons[overlapExon], overlapExon))
                    
                elif mode == "within":
                    (withinGns, withinExons) = self.withinExons(type, exStart, exEnd, detailed = 1)
                    for withinGn in withinGns:
                        if genes.has_key(withinGn):
                            genes[withinGn] += 1
                        else:
                            genes[withinGn]  = 1
                    for withinExon in withinExons.iterkeys():
                        exons[(gnId, exStart + ":E:" + exEnd)] = (withinExons[withinExon], withinExon)
                    if len(withinExons.keys()) > 0:
                        fh.write("%s\t%s\t%s\t%s\t%s\n" % (gnId, exStart, exEnd,\
                        #",".join(str(withinExons[withinExon]) for withinExon in withinExons),\
                        ",".join("%s:%s:%s" % withinExons[withinExon] for withinExon in withinExons),\
                        ",".join(withinExon for withinExon in withinExons)))
                        #",".join("%s:%s:%s" % withinExons[withinExon] for withinExon in withinExons),\
                   
                    
                          
                elif mode not in ["bothSS", "singleSS", "overlap", "exon", "within"]:
                    self.__inform("map2gene: Ignoring unknown mode %s\n" % mode)
                
            
        elif event_type == "A5SS" or event_type == "A3SS":
            pass
        elif event_type == "MXE":
            pass
        
        fh.close()


    def map2region(self, sg, type):
        """get annotations for the regions spaned by a SpliceGraph
        returns the annotation dictionary for all element within the specified region
        with complete information, such as chr, strand, start, stop etc...
        --> (dict)
        """
        
        if not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)

        if isinstance(sg, SpliceGraph.SpliceGraph):
            
            # get annotation
            if not self.annotationFiles.has_key(type): # give up
                self.__inform("ERROR in map2region: unknown annotation type [%s]" % type)
                raise Errors.ArgumentError("map2region", "unknown annotation type [%s]" % type)

            elif not self.exonInfo.has_key(type):    # parse it
                self.parseAnnotation(type)

            # get region
            #regStart, regEnd = sg.genomicRange()
            ex = sg.allExons()
            regStart = ex[0][0]
            regEnd = ex[-1][1]

            # get overlap
            genes = {}
            (overlapGns, overlapExons) = self.overlappingGenesExons(type,regStart,regEnd,detailed=1)

            for overlapGn in overlapGns:
                if genes.has_key(overlapGn):
                    genes[overlapGn] += 1
                else:
                    genes[overlapGn]  = 1
                    
        return genes, overlapExons


    def mapAllInDir(self, indir, outfile, type, method):
        "Map all *.sg files in 'indir' to 'type' by 'method', write to 'outfile' --> bool"

        if not os.path.isdir(indir):
            self.__inform("mapAllInDir: skipping %s: directory is not found.\n" % indir)
            return False

        else:
            # get file list in directory
            infiles = [ f for f in os.listdir(indir) if f.endswith('.sg')]

            #output files with distances of neighboring sites
            try:
                if os.path.exists("%s.exons"%outfile):
                    os.remove("%s.exons"%outfile)
                if os.path.exists(outfile):
                    os.remove(outfile)
                outFH = open(outfile,'w')
                outFHexons = open("%s.exons"%outfile,"w")
                
            except:
                self.__inform("mapAllInDir: skipping %s: outfile '%s' problem: %s\n" % (indir, outfile, sys.exc_info()[1]))
                return False
            else:
                outFH.write("#mapping for SpliceGraph gene models, generated by Mapper.py, indir: %s\n" % indir)
                outFH.write("#exon match mode: %s\n" % method)
                outFH.write("#mappting to: %s\n" % type)
                outFH.write("#gnId\tnbExons\tnbMappings\tMappingIds\tMappingTypes\tnbMatchingExons\tmaxNbMatchingExons\tnbConnectedComponents\n")

                outFHexons.write("#mapping for SpliceGraph gene models, generated by Mapper.py, indir: %s\n" % indir)
                outFHexons.write("#exon match mode: %s\n" % method)
                outFHexons.write("#mappting to: %s\n" % type)
                outFHexons.write("#gnId\texonId\t%s_gnId\t%s_exonId\n" % (type, type))
                
                for i in xrange(len(infiles)):

                    infile = infiles[i]

                    sg = SpliceGraph.SpliceGraph(filename=indir+'/'+infile, name=infile[0:-3])

                    gnIds, exons   = self.map2gene(sg, type, method)
                    gnInfo  = {}
                    maxNbEx = 0
                    for k in gnIds.iterkeys():
                        gnInfo[k] = self.getGeneInfo(k, type, 'biotype')
                        if gnIds[k] > maxNbEx:
                            maxNbEx = gnIds[k]

                    outFH.write('\t'.join([sg.name,                                              \
                                           str(sg.nbExons(internalOnly=False)),                  \
                                           str(len(gnIds.keys())),                               \
                                           ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                           ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                           ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                           str(maxNbEx),                                         \
                                           str(len(sg.connectedComponents())),                   \
                                           ])+'\n')
                    for exon in exons.iterkeys():
                        outFHexons.write('\t'.join([exon[0], exon[1], exons[exon][0], exons[exon][1]]) + '\n')
                outFHexons.close()
                outFH.close()
                os.chmod("%s.exons"%outfile, 0660)
                os.chmod(outfile, 0660)
                return True

    def __repr__(self):
        "generate a string representation of a Mapper object (automatically called, e.g. by print) --> string"

        return repr({'species':self.species,'assemby':self.assembly,'types':self.exonInfo.keys()})


    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            self.log.write(txt)
            self.log.flush()


if __name__ == "__main__":
    import Mapper, SpliceGraph
    print Mapper.__doc__

    sgNames = ['chr21:39906621:-','chr21:45184188:-','chr21:9839812:+','chr21:30234222:-',
               'chr21:25725884:-','chr21:35239207:-','chr21:42401630:-','chr21:46703021:-',
               'chr21:32555767:-','chr21:38550534:+','chr21:44803071:-','chr21:25856376:+',
               'chr21:35333593:-','chr21:42492868:+','chr21:46703031:-','chr21:32573247:-',
               'chr21:36181239:+','chr21:38564543:+','chr21:42527933:-','chr21:44818034:+',
               'chr21:25863299:-','chr21:33364320:+','chr21:45540940:-',
               ]
    
    mapper  = Mapper.Mapper(verbose=1)

    for sgName in sgNames:
        sg = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr21/%s.sg' % sgName)

        print '%s correspondances (%i exons):' % (sgName, sg.nbExons(internalOnly=False))
        
        gnIds, exons = mapper.map2gene(sg, 'ensembl', 'exon')
        gnInfo = {}
        for k in gnIds.iterkeys():
            gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
        print '\tmatching bothSS  : ensembl gene(s) %s' % \
              (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
        print "\texon information"
        print '\n'.join("\t\t%s <- %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 
        
        gnIds, exons = mapper.map2gene(sg, 'ensembl', 'bothSS')
        gnInfo = {}
        for k in gnIds.iterkeys():
            gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
        print '\tmatching bothSS  : ensembl gene(s) %s' % \
              (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
        print "\texon information"
        print '\n'.join("\t\t%s <- %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 
        
        gnIds, exons = mapper.map2gene(sg, 'ensembl', 'singleSS')
        gnInfo = {}
        for k in gnIds.iterkeys():
            gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
        print '\tmatching singleSS: ensembl gene(s) %s' % \
              (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
        print "\texon information"
        print '\n'.join("\t\t%s <- %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 

        gnIds, exons = mapper.map2gene(sg, 'ensembl', 'overlap')
        gnInfo = {}
        for k in gnIds.iterkeys():
            gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
        print '\tmatching overlaps: ensembl gene(s) %s' % \
              (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
        print "\texon information"
        print '\n'.join("\t\t%s -> %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 
        
        
##         gnIds, exons = mapper.map2gene(sg, 'vega', 'bothSS')
##         gnInfo = {}
##         for k in gnIds.iterkeys():
##             gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
##         print '\tmatching bothSS  : vega    gene(s) %s' % \
##               (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
##         print "\texon information"
##         print '\n'.join("\t\t%s -> %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 

##         gnIds, exons = mapper.map2gene(sg, 'vega', 'singleSS')
##         gnInfo = {}
##         for k in gnIds.iterkeys():
##             gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
##         print '\tmatching singleSS: vega    gene(s) %s' % \
##               (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
##         print "\texon information"
##         print '\n'.join("\t\t%s -> %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 

##         gnIds, exons = mapper.map2gene(sg, 'vega', 'overlap')
##         gnInfo = {}
##         for k in gnIds.iterkeys():
##             gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
##         print '\tmatching overlaps: vega    gene(s) %s' % \
##               (','.join("%s:%s[%i]" % (k, gnInfo[k], gnIds[k]) for k in gnIds.iterkeys()))
##         print "\texon information"
##         print '\n'.join("\t\t%s -> %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 

        gnIds, exons = mapper.map2gene(sg, 'hgu95av2', 'overlap')
        print '\tmatching overlap : hgu95av2(s)     %s' % \
              (','.join("%s:[%i]" % (k, gnIds[k]) for k in gnIds.iterkeys()))
        print "\texon information"
        print '\n'.join("\t\t%s -> %s" % (exon[1], exons[exon]) for exon in exons.iterkeys() ) 

##     #work, but takes sometime to run 
##     res = raw_input("test the mousemapper? (y/n): ")
##     if res.startswith('yes'):
##         mousemapper = Mapper.Mapper(species="mmu",verbose=1)
##         indir = "/r100/burge/shared/splice_graphs/mm6/FlatFiles/SpliceGraphsFiltered_chr1"
##         outfile = "mm6.chr1.genrate.map"
##         mousemapper.mapAllInDir( indir, outfile, "genrate", method="overlap")
        
