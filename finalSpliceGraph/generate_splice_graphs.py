"""Framework for generating splice graphs

Synopsis:

"""

# standard lib
import os
import sys
import time
import shelve
from array import array
from stat import S_ISGID, S_IRWXU, S_IRWXG

# 3rd party libraries
import MySQLdb

# project related libraries
import configuration
import SpliceGraph
import Classifier


class GenerateSpliceGraphs:
    "GenerateSpliceGraphs: Object that generates the splice-graphs"

    defaultSpecies = "hsa"
    currentChromosome = None

    #standardEvents = ['CE','SE','A5SS','A3SS','CI','RI','MXE','MSE']
    #classification of standard events is now called in do_complete_pipeline.py

    def __init__(self, MySQL_information=None, conf=None, species=None, verbose=False, useShelve=False, log=sys.stderr):
        #set verbose mode
        self.verbose = verbose
        self.log     = log

        #store data structures to files?
        self.__useShelve = useShelve

        #init configuration
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = configuration.Configuration()

        #define species
        if species is None:
            self.species = self.defaultSpecies
        elif species in self.conf.getSpeciesList():
            self.species = species
        else:
            self.__inform("ERROR: species %s is not in configuration file %s\n" % (species, self.conf.configFile()))

        #store MySQL
        if isinstance(MySQL_information, dict):
            self.MySQL_information = MySQL_information
            self.connect_to_MySQL()
        elif MySQL_information is None:
            #get MySQL login from configuration file for default species
            if self.species in self.conf.getSpeciesList():
                self.MySQL_information = {"user": os.getlogin(),
                                          "DB": self.conf[self.species]["Assembly"],
                                          "PASSWORD": '',
                                          "host": self.conf[self.species]["MySQLHost"],
                                          "knownGene": self.conf[self.species]["knownGene"],
                                          "mRNA": self.conf[self.species]["mRNA"],
                                          "MetaTable": self.conf[self.species]["MetaTable"],
                                          "Chromosomes": self.conf[self.species]["Chromosomes"].split(","),
                                          }
                self.connect_to_MySQL()
            else:
                self.__inform("ERROR: Unable to get MySQL information in GenerateSpliceGraphs.__init__()\n")
                sys.exit(0)


    def printMySQL_information(self):
        "Prints the information in GenerateSpliceGraphs.MySQL_information --> None"

        print >> sys.stdout, 'MySQL information:'
        print >> sys.stdout, "\n".join(["%15s = %s" % (k, self.MySQL_information[k]) for k in self.MySQL_information.keys()])
        print >> sys.stdout


    def connect_to_MySQL(self):
        "Establish the connection with MySQL --> Boolean"

        self.db = MySQLdb.connect(db = self.MySQL_information['DB'],
                                  user = self.MySQL_information['user'],
                                  host=self.MySQL_information['host'],
                                  passwd=self.MySQL_information['PASSWORD'])
        
        
    def generate_splice_graphs(self):
        "Method that generates all splice graphs"
        # initiate MySQL connection
        self.cursor = self.db.cursor()
        
        self.mSS_files = []
        self.geneModels_files = []
        self.mSS_transcript_ids_files = []

        for chromosome in self.MySQL_information['Chromosomes']:
            self.currentChromosome = chromosome
            self.generate_splice_graphs_for_chromosome(chromosome)

        # end MySQL connection
        self.cursor.close()


    def generate_splice_graphs_for_chromosome(self, chromosome):
        "Method that generates all splice graphs for a chromosome"


        # Three main data objects
        self.mSS = {}
        self.elements = {'TSS':{},
                         '5ss':{},
                         '3ss':{},
                         'TER':{}}
        self.genes = {}

        # save information on transcript already used, redundancy check
        self.transcript_ids_encountered = {}
        nrRedundantTranscripts = 0
        # save all transcript ids that use a 'element -> element link'
        self.mSS_transcript_ids = {}
        
        self.__inform('Generating splice graphs for chromosome: %s\n' % chromosome)


        # Data set 1. KnownGenes
        self.__inform('Retreiving KnownGenes...\n')
        self.cursor.execute("""SELECT strand, exonCount, exonStarts, exonEnds, name from %s WHERE chrom = "%s";""" % \
                            (self.conf[self.species]["knownGene"], chromosome))
        res = self.cursor.fetchall()
        self.__inform('\tNumber of knownGenes read: %d\n' % len(res))

        # Remove contaminated EST libraries, if YES, first retrieve all accession id from the conataminated libraries
        if self.conf[self.species]["RemoveContaminatedESTlibraries"] == 'YES':
            self.RemoveContaminatedESTlibraries = True
            self.cursor.execute("""SELECT name from %s, Transcript2library NATURAL JOIN MetaTable \
            WHERE %s.name = Transcript2library.accession AND %s.chrom = "%s" AND MetaTable.contaminated = 1;""" % \
                                (self.conf[self.species]["knownGene"],self.conf[self.species]["knownGene"],
                                 self.conf[self.species]["knownGene"], chromosome))
            self.contaminated = {}
            contaminated = self.cursor.fetchall()
            for name in contaminated:
                self.contaminated[name[0]] = 1
            self.__inform('\t%i knownGene(s) from contaminated libraries excluded\n' % len(self.contaminated),3)
        else:
            self.contaminated = {}
            self.RemoveContaminatedESTlibraries = False
            

        knownGenes = []

        for tx in res:

            chr = chromosome
            strand = tx[0]
            exonCount = int(tx[1])
            name = tx[4]

            # redundancy check
            if self.transcript_ids_encountered.has_key(name):
                self.transcript_ids_encountered[name] += 1
                nrRedundantTranscripts += 1
                continue
            else:
                self.transcript_ids_encountered[name] = 1

            if self.RemoveContaminatedESTlibraries:
                if self.contaminated.has_key(name):
                    # then, dont use this contaminated sequence
                    continue
            
            if isinstance(tx[2], array):
                # for newer MySQLdb versions that return an array.array object
                exonStarts = tx[2].tostring().split(',')
            else:
                # older version that return a string
                exonStarts = tx[2].split(',')
                
            if isinstance(tx[3], array):
                exonEnds = tx[3].tostring().split(',')
            else:
                exonEnds = tx[3].split(',')

            knownGenes.append((exonStarts, exonEnds, exonCount, strand, name))
            
        self.__inform('\tNumber of redundant knownGenes: %i\n\n' % nrRedundantTranscripts,2)
        self.__inform("\tPreprocessing knownGenes\n")
        knownGenes = self.preprocess_ESTs(knownGenes)
        self.__inform("\n\tAdding knownGenes to mSS and elements...\n")
        self.add_transcripts_to_splice_matrix(knownGenes, chromosome)
        self.__checkTranscriptsVersusSplicingSiteMatrix(knownGenes,self.mSS,self.elements,chromosome,checkEnds=True)

        '''
        # Data set 2. mRNAs
        self.__inform('Retreiving mRNAs...\n')
        self.cursor.execute("""SELECT strand, blockCount, tStarts, blockSizes, Qname from %s WHERE tName = "%s";""" % \
                            (self.conf[self.species]["mRNA"], chromosome))
        res = self.cursor.fetchall()
        self.__inform('\tNumber of mRNAs read: %d\n' % len(res))
        
        # Remove contaminated EST libraries, if YES, first retrieve all accession id from the contaminated libraries
        if self.RemoveContaminatedESTlibraries is True:
            self.cursor.execute("""SELECT qName from %s, Transcript2library NATURAL JOIN MetaTable \
            WHERE %s.qName = Transcript2library.accession AND %s.tName = "%s" AND MetaTable.contaminated = 1;""" % \
                                (self.conf[self.species]["mRNA"],self.conf[self.species]["mRNA"],
                                 self.conf[self.species]["mRNA"], chromosome))
            self.contaminated = {}
            contaminated = self.cursor.fetchall()
            for name in contaminated:
                self.contaminated[name[0]] = 1
            self.__inform('\t%i  mRNA(s) from contaminated libraries excluded\n' % len(self.contaminated),3)
        else:
            self.contaminated = {}

        mRNAs = []
        nrRedundantTranscripts = 0

        for tx in res:

            chr = chromosome
            strand = tx[0]
            exonCount = int(tx[1])
            name = tx[4]

            # redundancy check
            if self.transcript_ids_encountered.has_key(name):
                self.transcript_ids_encountered[name] += 1
                nrRedundantTranscripts += 1
                continue
            else:
                self.transcript_ids_encountered[name] = 1

            if self.RemoveContaminatedESTlibraries:
                if self.contaminated.has_key(name):
                    # then, dont use this contaminated sequence
                    continue
            
            if isinstance(tx[2], array):
                # for newer MySQLdb versions that return an array.array object
                exonStarts = tx[2].tostring().split(',')
            else:
                # older version that return a string
                exonStarts = tx[2].split(',')
                
            if isinstance(tx[3], array):
                exonSizes = tx[3].tostring().split(',')
                exonEnds = []
                for exon in xrange(exonCount):
                    exonEnds.append(int(exonStarts[exon])+int(exonSizes[exon]))
            else:
                exonSizes = tx[3].split(',')
                exonEnds = []
                for exon in xrange(exonCount):
                    exonEnds.append(int(exonStarts[exon])+int(exonSizes[exon]))

            mRNAs.append((exonStarts, exonEnds, exonCount, strand, name))
            
        self.__inform('\tNumber of redundant mRNAs: %i\n\n' % nrRedundantTranscripts,2)
        self.__inform("\tPreprocessing mRNAs\n")
        mRNAs = self.preprocess_ESTs(mRNAs)
        self.__inform("\n\tAdding mRNAs to mSS and elements...\n")
        self.add_transcripts_to_splice_matrix(mRNAs, chromosome)
        self.__checkTranscriptsVersusSplicingSiteMatrix(mRNAs,self.mSS,self.elements,chromosome,checkEnds=True)
        
        # Data set 3. intronESTs
        self.__inform('Retreiving intronESTs...\n')
        qNumInsert = 10
        tNumInsert = 10
        self.cursor.execute("""SELECT strand, blockCount, tStarts, blockSizes, Qname from %s%s;""" % \
                            (chromosome, self.conf[self.species]["ESTsuffix"]))
        res = self.cursor.fetchall()
        self.__inform("\tNumber of intronEST's read: %d\n" % len(res))

        # Remove contaminated EST libraries, if YES, first retrieve all accession id from the conataminated libraries
        if self.RemoveContaminatedESTlibraries is True:
            self.cursor.execute("""SELECT qName from %s%s, Transcript2library NATURAL JOIN MetaTable \
            WHERE %s%s.qName = Transcript2library.accession AND MetaTable.contaminated = 1;""" % \
                                (chromosome,self.conf[self.species]["ESTsuffix"],chromosome,self.conf[self.species]["ESTsuffix"]))
            self.contaminated = {}
            contaminated = self.cursor.fetchall()
            for name in contaminated:
                self.contaminated[name[0]] = 1
            self.__inform("\t%i intronEST(s) from contaminated libraries excluded\n" % len(self.contaminated),3)
        else:
            self.contaminated = {}
 
        intronEsts = []
        nrRedundantTranscripts = 0
        
		####EDIT### This is where we remove transcripts not wanted
        rmtxfile = open('TranscriptExcludeList.txt')
        rmtx = rmtxfile.read().strip()
        removeList = rmtx.split('\t')
		
        for tx in res:

            chr = chromosome
            strand = tx[0]
            exonCount = int(tx[1])
            name = tx[4]

            # redundancy check
            if self.transcript_ids_encountered.has_key(name):
                self.transcript_ids_encountered[name] += 1
                nrRedundantTranscripts += 1
                continue
            else:
                self.transcript_ids_encountered[name] = 1

            if self.RemoveContaminatedESTlibraries:
                if self.contaminated.has_key(name):
                    # then, dont use this contaminated sequence
                    continue
            
            ###EDIT### remove transcipts
            if name in removeList:
            	print 'this transcript was not included in splice graph: %s' % name
            	continue
            
            if isinstance(tx[2], array):
                # for newer MySQLdb versions that return an array.array object
                exonStarts = tx[2].tostring().split(',')
            else:
                # older version that return a string
                exonStarts = tx[2].split(',')
                
            if isinstance(tx[3], array):
                exonSizes = tx[3].tostring().split(',')
                exonEnds = []
                for exon in xrange(exonCount):
                    exonEnds.append(int(exonStarts[exon])+int(exonSizes[exon]))
            else:
                exonSizes = tx[3].split(',')
                exonEnds = []
                for exon in xrange(exonCount):
                    exonEnds.append(int(exonStarts[exon])+int(exonSizes[exon]))

            intronEsts.append((exonStarts, exonEnds, exonCount, strand, name))
            
        self.__inform('\tNumber of redundant intronEsts: %i\n\n' % nrRedundantTranscripts,2)
        self.__inform("\tPreprocessing intronESTs\n")
        intronEsts = self.preprocess_ESTs(intronEsts)
        self.__inform("\n\tAdding intronESTs to mSS and elements...\n")
        self.add_transcripts_to_splice_matrix(intronEsts, chromosome, skip_ends=True)
        self.__checkTranscriptsVersusSplicingSiteMatrix(intronEsts,self.mSS,self.elements,chromosome,checkEnds=False)
        '''
        self.__checkElements()
         
        self.__inform("Generating gene models\n")
        self.generateGeneModels()
        self.__inform("\ngenerated %i gene models\n" % len(self.genes.keys()))
        self.__checkGeneModels()

        self.__inform("Storing orignial splice graphs...")
        self.storeSpliceGraphs() #will store flat files of raw data based splice graphs
        self.__inform("finished.\n")

        # classification will be called from do_complete_pipeline.py
        #self.__inform("Classifying standard events (%s) in original splice graphs..." % ', '.join(self.standardEvents) )
        #self.classifyStandardEvents(indir= "%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "SpliceGraphs", self.currentChromosome),
        #                            outdir="%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "StandardEvents", self.currentChromosome))
        #self.__inform("finished.\n")
        #self.__inform("Classifying standard events (%s) in filtered splice graphs..." % ', '.join(self.standardEvents) )
        #self.classifyStandardEvents(indir= "%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "SpliceGraphsFiltered", self.currentChromosome),
        #                            outdir="%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "StandardEventsFiltered", self.currentChromosome))
        #self.__inform("finished.\n")

        shelveHandle = None
        if self.__useShelve:
            
            self.mSS_files.append('%smSS_%s' % (self.conf[self.species]["FlatFilesStoragePath"],chromosome))
            if os.path.exists(self.mSS_files[-1]):
                os.remove(self.mSS_files[-1])
            shelveHandle = shelve.open(self.mSS_files[-1], writeback=True)
            os.chmod('%smSS_%s' % (self.conf[self.species]["FlatFilesStoragePath"], \
                                   chromosome),S_ISGID | S_IRWXU | S_IRWXG)

            shelveHandle['mSS'] = {}
            for key in self.mSS.keys():
                shelveHandle['mSS'][key] = self.mSS[key]            
            shelveHandle['elements'] = {}
            for el in self.elements.keys():
                shelveHandle['elements'][el] = self.elements[el]
            shelveHandle.close()

            # save the geneModels in a shelve file
            self.geneModels_files.append('%sgeneModels_%s' % (self.conf[self.species]["FlatFilesStoragePath"],chromosome))
            if os.path.exists(self.geneModels_files[-1]):
                os.remove(self.geneModels_files[-1])
            shelveHandle = shelve.open(self.geneModels_files[-1], writeback=True)
            os.chmod(self.geneModels_files[-1],S_ISGID | S_IRWXU | S_IRWXG)
            for id in self.genes.keys():
                shelveHandle[id] = self.genes[id]
            shelveHandle.close()

            # save the transcript ids in mSS to shelve file
            self.mSS_transcript_ids_files.append('%smSS_transcript_ids_files_%s' % (self.conf[self.species]["FlatFilesStoragePath"],chromosome))
            if os.path.exists(self.mSS_transcript_ids_files[-1]):
                os.remove(self.mSS_transcript_ids_files[-1])
            shelveHandle = shelve.open(self.mSS_transcript_ids_files[-1], writeback=True)
            os.chmod(self.mSS_transcript_ids_files[-1],0770)
            
            for element1 in self.mSS_transcript_ids.keys():
                shelveHandle[element1] = self.mSS_transcript_ids[element1]
            shelveHandle.close()
 


    def __checkTranscriptsVersusSplicingSiteMatrix(self, transcripts, mSS, elements, chromosome, checkEnds=True):
        correctTranscripts = 0
        falseTranscripts = 0
        self.__inform("\tCheck Transcripts Versus mSS\n",2)
        for transcript in transcripts:
            res = self.__checkTranscriptVersusSplicingSiteMatrix(transcript, mSS, elements,chromosome, checkEnds)
            if res is True:
                correctTranscripts += 1
            else:
                falseTranscripts += 1
        self.__inform("\tTranscripts with full path in mSS: %i\n\tTranscripts where path is broken: %i\n\n" % (correctTranscripts, falseTranscripts),2)
        
    def __checkTranscriptVersusSplicingSiteMatrix(self, transcript, mSS, elementsDict,chromosome, checkEnds):
        "Internal consistency: checking that transcript splice sites have a path through mSS -> Boolean"
        strand = transcript[3]
        exonCount =  int(transcript[2])
        exonStarts = transcript[0]
        exonEnds = transcript[1]
        

        elements = []
        # generating a list of elements for transcript
        for exon_nr in xrange(exonCount):
            element1 = '%s:%i:%s' % (chromosome, int(exonStarts[exon_nr])+1,strand)
            element2 = '%s:%i:%s' % (chromosome, int(exonEnds[exon_nr]),strand)
            if exon_nr == 0:
               if checkEnds is True: 
                   elements.append(element1)
               else:
                   pass
               elements.append(element2)
            elif exon_nr == exonCount - 1:
                elements.append(element1)
                if checkEnds is True: 
                    elements.append(element2)
                else:
                    pass
            else:
                elements.append(element1)
                elements.append(element2)
        nrElements = len(elements)

        if strand == '+':
            sequence = range(nrElements)
        else:
            sequence = range(nrElements)
            sequence.reverse()
     
        self.__inform("NrExons: %s, %s\n" % (exonCount, nrElements),5)
        # checking that the path through transcript exist in mSS
        for element_id in sequence:
            self.__inform("Element: %s\n" % elements[element_id],5)
            if strand == '+' and element_id < nrElements -1:
                if mSS.has_key(elements[element_id]) is False:
                    self.__inform("Element in transcript not found in mSS: %s (%i of %i)\n" % \
                                  (elements[element_id],element_id,nrElements), 3)

                    incrementor = 0
                    self.__inform("Name: %s\n" % transcript[4],3)
                    for el in elements:
                        try:
                            self.__inform("%i %s %s\n" % (incrementor, el, mSS[el]),3)
                        except:
                            self.__inform("%i %s\n" % (incrementor, el),3)
                        incrementor += 1
                    return False
                if mSS[elements[element_id]].has_key(elements[element_id+1]) is False:
                    self.__inform("Map between element1 and element2 not found in mSS: %s, %s = %i,%i\n" \
                                  % (elements[element_id], elements[element_id+1],\
                                                       element_id,element_id+1), 3)
                    incrementor = 0
                    self.__inform("Name: %s\n" % transcript[4],3)
                    for el in elements:
                        try:
                            self.__inform("%i %s %s\n" % (incrementor, el, mSS[el]),3)
                        except:
                            self.__inform("%i %s\n" % (incrementor, el),3)
                        incrementor += 1
                    return False
            elif strand == '-' and element_id > 0:
                if mSS.has_key(elements[element_id]) is False:
                    self.__inform("Element in transcript not found in mSS: %s (%i of %i)\n" % \
                                  (elements[element_id],element_id,nrElements), 3)
                    incrementor = 0
                    self.__inform("Name: %s\n" % transcript[4],3)
                    
                    for el in elements:
                        try:
                            self.__inform("%i %s %s\n" % (incrementor, el, mSS[el]),3)
                        except:
                            self.__inform("%i %s\n" % (incrementor, el),3)

                        incrementor += 1
                    return False
                if mSS[elements[element_id]].has_key(elements[element_id-1]) is False:
                    self.__inform("Map between element1 and element2 not found in mSS: %s, %s = %i,%i\n" \
                                  % (elements[element_id], elements[element_id-1],\
                                                       element_id,element_id-1), 3)
                    incrementor = 0
                    self.__inform("Name: %s\n" % transcript[4],3)
                    for el in elements:
                        try:
                            self.__inform("%i %s %s\n" % (incrementor, el, mSS[el]),3)
                        except:
                            self.__inform("%i %s\n" % (incrementor, el),3)
                        incrementor += 1
                    return False
        return True
                
  
    def preprocess_ESTs(self, ESTs):
        "extend blocks over gaps smaller than a set length (default= 32 bp, see MergeBlockThreshold in configuration.txt)"

        merge_block_threshold = int(self.conf[self.species]["MergeBlockThreshold"])
        min_est_end_length = int(self.conf[self.species]["MinESTEndLength"])

        new_ESTs = []
        for tx in ESTs:
            # input
            strand = tx[3]
            exonCount =  int(tx[2])
            exonStarts = tx[0]
            exonEnds = tx[1]
            name = tx[4]
            if name == 'ENSMUST00000103283': print tx

            # new output
            tStarts = [] # exonStarts
            tEnds = []
            blockCount = 0 # exonCount

            exon_nr = 0
            while exon_nr < exonCount:
                j = 1
                cumulativeExonEnd = int(exonEnds[exon_nr])
                while exon_nr+j < exonCount and \
                          int(exonStarts[exon_nr+j]) - cumulativeExonEnd  < merge_block_threshold:
                    cumulativeExonEnd = int(exonEnds[exon_nr+j])
                    j += 1
                blockCount += 1
                tStarts.append(exonStarts[exon_nr])
                tEnds.append(cumulativeExonEnd)
                exon_nr += j

                if abs(int(exonStarts[exon_nr-j]) - cumulativeExonEnd ) < 2:
                    self.__inform('\tPreprocessing encountered short exon fragment: %s-%s in %s\n' % \
                                  (exonStarts[exon_nr-j], cumulativeExonEnd, name ),3 )
                    self.__inform('ExonStarts and ExonEnds\n',4)
                    self.__inform(exonStarts,4)
                    self.__inform("\n",4)
                    self.__inform(exonEnds,4)
                    self.__inform("\n",4)
                    self.__inform(tStarts,4)
                    self.__inform("\n",4)
                    self.__inform(tEnds,4)
                    self.__inform('\n',4)
                    
                if exon_nr-j-1 >= 0:
                    if abs(int(exonStarts[exon_nr-j]) - int(exonEnds[exon_nr-j-1]) ) < 32:
                        self.__inform('\tPreprocessing encountered short intron: %s-%s in %s\n' % \
                                      (int(exonEnds[exon_nr-j-1]), exonStarts[exon_nr-j], name ),3 )
                        self.__inform('ExonStarts and ExonEnds\n',4)
                        self.__inform(exonStarts,4)
                        self.__inform("\n",4)
                        self.__inform(exonEnds,4)
                        self.__inform("\n",4)
                        self.__inform(tStarts,4)
                        self.__inform("\n",4)
                        self.__inform(tEnds,4)
                        self.__inform('\n',4)
                 
            
            # checks that ends are longer than MinESTEndLength as defined in configuration.txt
            # otherwise exclude ends
            try:
                x = int(tEnds[0]) - int(tStarts[0])
            except:
            	continue
            
            if int(tEnds[0]) - int(tStarts[0]) < min_est_end_length:
                self.__inform("EST end shorter than minimal end length (threshold=%i): %i\n" % (min_est_end_length,
                                                                                              int(tEnds[0]) - int(tStarts[0])),3)

                del tEnds[0]
                del tStarts[0]
                blockCount -= 1
            
            try:
                x = int(tEnds[-1]) - int(tStarts[-1])
            except:
            	continue
            
            if int(tEnds[-1]) - int(tStarts[-1]) < min_est_end_length:
                self.__inform("EST end shorter than minimal end length (threshold=%i): %i\n" % (min_est_end_length,
                                                                                              int(tEnds[-1]) - int(tStarts[-1])),3)
                del tEnds[-1]
                del tStarts[-1]
                blockCount -= 1

            # collects the new intronEst exon blocks
            new_ESTs.append((tStarts, tEnds, blockCount, strand, name))
        
        return new_ESTs


    def add_transcripts_to_splice_matrix(self, transcripts, chromosome, skip_ends=False):
        # Input transcripts is a tuple of [[exonStarts],[exonEnds], exonCount, strand]
        elements = self.elements
        mSS = self.mSS
        mSS_transcript_ids = self.mSS_transcript_ids

        for tx in transcripts:

            strand = tx[3]
            exonCount =  int(tx[2])
            exonStarts = tx[0]
            exonEnds = tx[1]
            name = tx[4]
                         
            for exon_nr in xrange(exonCount): 
                element1 = '%s:%i:%s' % (chromosome, int(exonStarts[exon_nr])+1,strand)
                element2 = '%s:%i:%s' % (chromosome, int(exonEnds[exon_nr]),strand)

                #1 SAVING ELEMENTS information
                #element1

                if exon_nr == 0:
                    if skip_ends:
                        # for ESTs where we don't care about start and end position
                        pass
                    else:
                        if strand == '+':
                            if elements['TSS'].has_key(element1):
                                elements['TSS'][element1] += 1
                            else:
                                elements['TSS'][element1] = 1
                        else:
                            if elements['TER'].has_key(element1):
                                elements['TER'][element1] += 1
                            else:
                                elements['TER'][element1] = 1
                else:
                    if strand == '+':
                        if elements['3ss'].has_key(element1):
                            elements['3ss'][element1] += 1
                        else:
                            elements['3ss'][element1] = 1
                    else:
                        if elements['5ss'].has_key(element1):
                            elements['5ss'][element1] += 1
                        else:
                            elements['5ss'][element1] = 1
                #element2
                if exon_nr == exonCount - 1:
                    if skip_ends:
                        pass
                    else:
                        if strand == '+':
                            if elements['TER'].has_key(element2):
                                elements['TER'][element2] += 1
                            else:
                                elements['TER'][element2] = 1
                        else:
                            if elements['TSS'].has_key(element2):
                                elements['TSS'][element2] += 1
                            else:
                                elements['TSS'][element2] = 1
                else:
                    if strand == '+':
                        if elements['5ss'].has_key(element2):
                            elements['5ss'][element2] += 1
                        else:
                            elements['5ss'][element2] = 1
                    else:
                        if elements['3ss'].has_key(element2):
                            elements['3ss'][element2] += 1
                        else:
                            elements['3ss'][element2] = 1
                            
                #2 Saving mSS information, exons
                if ( exon_nr == 0 and skip_ends ) or ( exon_nr == exonCount-1 and skip_ends ) :
                    pass
                else:
                    if strand == '+':
                        if mSS.has_key(element1):
                            if mSS[element1].has_key(element2):
                                mSS[element1][element2] += 1
                                mSS_transcript_ids[element1][element2].append(name) ###
                            else:
                                mSS[element1][element2] = 1
                                mSS_transcript_ids[element1][element2] = [name] ###
                        else:
                            mSS[element1] = {element2:1}
                            mSS_transcript_ids[element1] = {element2 : [name]} ###
                    else:
                        if mSS.has_key(element2):
                            if mSS[element2].has_key(element1):
                                mSS[element2][element1] += 1
                                mSS_transcript_ids[element2][element1].append(name) ###
                            else:
                                mSS[element2][element1] = 1
                                mSS_transcript_ids[element2][element1] = [name] ###
                        else:
                            mSS[element2] = {element1:1}
                            mSS_transcript_ids[element2] = {element1 : [name]} ###
                                
                #3 Saving mSS information, introns
                if exon_nr <= exonCount - 2:
                    next_exon_element1 = '%s:%i:%s' % (chromosome, int(exonStarts[exon_nr+1])+1,strand)

                    #enforce minimum intron length
                    if abs(int(element2.split(':')[1]) - int(next_exon_element1.split(':')[1])) >= int(self.conf[self.species]["MinIntronLength"]):
                        if strand == '+':
                            if mSS.has_key(element2):
                                if mSS[element2].has_key(next_exon_element1):
                                    mSS[element2][next_exon_element1] -= 1
                                    mSS_transcript_ids[element2][next_exon_element1].append(name) ###
                                else:
                                    mSS[element2][next_exon_element1] = -1
                                    mSS_transcript_ids[element2][next_exon_element1] = [name] ###
                            else:
                                mSS[element2] = {next_exon_element1:-1}
                                mSS_transcript_ids[element2] = {next_exon_element1 : [name]} ###
                        else:
                            if mSS.has_key(next_exon_element1):
                                if mSS[next_exon_element1].has_key(element2):
                                    mSS[next_exon_element1][element2] -= 1
                                    mSS_transcript_ids[next_exon_element1][element2].append(name) ###
                                else:
                                    mSS[next_exon_element1][element2] = -1
                                    mSS_transcript_ids[next_exon_element1][element2] = [name] ###
                            else:
                                mSS[next_exon_element1] = {element2:-1}
                                mSS_transcript_ids[next_exon_element1] = {element2 : [name]} ###
                        
        self.__inform('\tNr of keys in mSS: %d\n' % len(mSS.keys()),3)

    def get_transcripts_using_link(self,element1,element2):
        if self.mSS_transcript_ids.has_key(element1):
            if self.mSS_transcript_ids[element1].has_key(element2):
                return self.mSS_transcript_ids[element1][element2]
            else:
                self.__inform("Element2 (%s) is not mapped from element1 in mSS_transcript_ids.\n" % element2)
                return 0
        else:
            self.__inform("Element1 (%s) not in mSS_transcript_ids.\n" % element1)
            return 0




    def __checkElements(self, elements=None):
        "Internal consistency checking for self.elements"

        if elements is None: elements = self.elements

        self.__inform("\nchecking elements\n", 2)
        ErrorsDetected = False

        #check uniqueness of types
        if True:
            self.__inform("\tchecking uniqueness of element types: ", 2)
            thisFailed = False
            nonUniqueElements = {}
            Types = elements.keys()
            nonAllowedTypes = {'TSS': ['5ss','TER'],
                               'TER': ['3ss','TSS'],
                               '5ss': ['3ss','TSS'],
                               '3ss': ['5ss','TER']}
            for thisType in Types:
                #otherTypes = filter(lambda x: x != thisType, Types)
                otherTypes = nonAllowedTypes[thisType]
                thisElements = elements[thisType].keys()

                for thisElement in thisElements:
                    for otherType in otherTypes:
                        if elements[otherType].has_key(thisElement):
                            thisFailed = True
                            if nonUniqueElements.has_key(thisElement):
                                if thisType not in nonUniqueElements[thisElement]:
                                    nonUniqueElements[thisElement].append(thisType)
                                if otherType not in nonUniqueElements[thisElement]:
                                    nonUniqueElements[thisElement].append(otherType)
                            else:
                                nonUniqueElements[thisElement] = [thisType, otherType]
                                
            if thisFailed:
                self.__inform("NOT ok", 2)

                ids = nonUniqueElements.keys()
                ids.sort(cmp=lambda x,y: cmp(x.split(":")[1], y.split(":")[1]))
                self.__inform("\n", 3)
                self.__inform("\n".join("\t\t%20s is %s" % (id, ", ".join(nonUniqueElements[id])) for id in ids), 3)

                self.__inform("\n", 2)
                ErrorsDetected = True
            else:
                self.__inform("ok\n", 2)

        # summary
        if ErrorsDetected:
            self.__inform("elements NOT all ok\n\n", 2)
        else:
            self.__inform("elements all ok\n\n", 2)


    def generateGeneModels(self):
        "Walk through splice site matrix mSS, starting at TSS elements, and collect all genes"

        mSS      = self.mSS
        elements = self.elements
        genes    = self.genes

        existing = {} # stores gene ids for elements that are already in a gene model

        # sort elements by chromosomal coordinate (all on same chromosome)
        #all_sorted = mSS.keys()
        #all_sorted.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])))
        TSS_sorted = elements["TSS"].keys()
        TSS_sorted.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])))
        
        # begin a new gene at each TSS, extend as far as possible
        for tss in TSS_sorted:
            newGene     = {}    # keys: elements associated to that gene
            elStack     = [tss] # stack of elements that need to be extended
            fuseGeneIds = {}    # gene ids of genes in self.genes to be fused with newGene
            (chr, coord, strand) = tss.split(':')
            coord       = int(coord)

            self.__inform("\textending TSS %s: \n" % tss, 5)

            while len(elStack) > 0:

                # take first item from elStack and associate it with newGene
                thisElement = elStack.pop(0)
                (thisChr, thisCoord, thisStrand) = thisElement.split(':')
                thisCoord   = int(thisCoord)

                if newGene.has_key(thisElement):
                    newGene[thisElement] += 1
                else:
                    newGene[thisElement] = 1
                
                    # if element thisElement is already part of another gene, store id of that gene
                    # for later fusing with newGene/removing from self.genes
                    if existing.has_key(thisElement):
                        if fuseGeneIds.has_key(existing[thisElement]):
                            fuseGeneIds[existing[thisElement]] += 1
                        else:
                            fuseGeneIds[existing[thisElement]] = 1
            
                    # where do we go from here?
                    #...get continuening elements and append to elStack
                    if mSS.has_key(thisElement):
                        nextElements = mSS[thisElement].keys()
                        for nextElement in nextElements:
                            (nextChr, nextCoord, nextStrand) = nextElement.split(':')
                            nextCoord = int(nextCoord)
                            if thisStrand == '+' and nextStrand == '+' and nextCoord - thisCoord > 0: 
                                elStack.append(nextElement)
                            elif thisStrand == '-' and nextStrand == '-' and nextCoord - thisCoord < 0:
                                elStack.append(nextElement)
                            else:
                                self.__inform("WARNING: Ignoring problematic extension: %s --> %s\n" % \
                                              (thisElement, nextElement), 5)

            self.__inform("\t%i elements in newGene, %i genes in fuseGeneIds\n" % (len(newGene.keys()), len(fuseGeneIds.keys())), 5)

            # when no extension is possible any more, add that gene to list of collected genes
            # first, fuse all fuseGeneIds into newGene
            for gId in fuseGeneIds.keys():
                for fuseElement in genes[gId]:
                    if newGene.has_key(fuseElement):
                        newGene[fuseElement] += genes[gId][fuseElement]
                    else:
                        newGene[fuseElement] = genes[gId][fuseElement]

            # remove the fuseGeneIds from self.genes
            for gnId in fuseGeneIds:
                discardedGene = genes.pop(gnId)

            # now, add newGene to self.genes and update existing
            sortedElements = newGene.keys()
            sortedElements.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])), \
                                reverse=(sortedElements[0].split(':')[2] == '-'))
            genes[sortedElements[0]] = newGene
            for newElement in newGene:
                existing[newElement] = sortedElements[0]

            self.__inform("found %i gene models so far (last one was fused with %i existing models)\r" % \
                      (len(genes), len(fuseGeneIds)), 2)
                

    def __checkGeneModels(self, genes=None, mSS=None, elements=None):
        "Internal consistency checking for self.genes"

        if genes    is None: genes    = self.genes
        if mSS      is None: mSS      = self.mSS
        if elements is None: elements = self.elements

        self.__inform("\nchecking gene models\n", 2)
        ErrorsDetected = False

        #for each 5'SS in a gene, get all connected 3'SS elements in the same gene and check the intron lengthes
        if True:
            self.__inform("\tchecking minimum intron length (%i bp): " % int(self.conf[self.species]["MinIntronLength"]), 2)
            thisFailed = False
            shortIntrons = {}
            for gene in genes.keys():
                fiveSSs  = filter(elements['5ss'].has_key, genes[gene].keys())
                threeSSs = filter(elements['3ss'].has_key, genes[gene].keys())

                for fiveSS in fiveSSs:
                    for threeSS in threeSSs:
                        if mSS.has_key(fiveSS) and mSS[fiveSS].has_key(threeSS):
                            if abs(int(fiveSS.split(':')[1]) - int(threeSS.split(':')[1])) < \
                                   int(self.conf[self.species]["MinIntronLength"]):
                                thisFailed = True
                                shortIntrons["%s to %s" % (fiveSS, threeSS)] = abs(int(fiveSS.split(':')[1]) - int(threeSS.split(':')[1]))
                                #break
                    #if thisFailed:
                        #break

            if thisFailed:
                self.__inform("NOT ok", 2)

                ids = shortIntrons.keys()
                ids.sort(cmp=lambda x,y: cmp(x.split(":")[1], y.split(":")[1]))
                self.__inform("\n", 3)
                self.__inform("\n".join("\t\t%s is %i bp" % (id, shortIntrons[id]) for id in ids), 3)

                self.__inform("\n", 2)
                ErrorsDetected = True
            else:
                self.__inform("ok\n", 2)

        #for each element, check whether it belongs to a single gene
        if True:
            self.__inform("\tchecking whether elements are associated to a single gene: ", 2)
            thisFailed = False
            foundIN    = {}
            ambiguous  = {}
            for gene in genes.keys():
                for el in genes[gene].keys():
                    if foundIN.has_key(el):
                        thisFailed = True
                        if not ambiguous.has_key(el):
                            ambiguous[el] = []
                        ambiguous[el].extend(gene, foundIN[el])
                        #break
                    else:
                        foundIN[el] = gene
                #if thisFailed:
                    #break
            foundIN.clear()
        
            if thisFailed:
                self.__inform("NOT ok", 2)
                self.__inform(": %s" % ambiguous, 3)
                self.__inform("\n", 2)
                ErrorsDetected = True
            else:
                self.__inform("ok\n", 2)

        #for each gene, check if all the elements are on the same strand
        if True:
            self.__inform("\tchecking whether elements in a gene are all on same strand: ", 2)
            thisFailed   = False
            problemGenes = {}
            for gene in genes.keys():
                gnElements = genes[gene].keys()
                sameStrand = filter(lambda x: x.split(':')[2] == gnElements[0].split(':')[2], gnElements)
                if len(gnElements) != len(sameStrand):
                    thisFailed = True
                    problemGenes[gene] = gnElements
                    #break
        
            if thisFailed:
                self.__inform("NOT ok", 2)
                self.__inform(": %s" % problemGenes, 3)
                self.__inform("\n")
                ErrorsDetected = True
            else:
                self.__inform("ok\n", 2)

        # summary
        if ErrorsDetected:
            self.__inform("gene models NOT all ok\n", 2)
        else:
            self.__inform("gene models all ok\n", 2)


    def storeSpliceGraphs(self):
        "store splice graphs as flat files suitable for visualization"

        #if genes              is None: genes              = self.genes
        #if mSS                is None: mSS                = self.mSS
        #if elements           is None: elements           = self.elements
        #if mSS_transcript_ids is None: mSS_transcript_ids = self.mSS_transcript_ids

        self.__inform("\nstoring splice graphs\n", 2)

        storagePath = "%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "SpliceGraphs", self.currentChromosome)
        #storagePathFiltered = "%s/%s_%s" % (self.conf[self.species]['FlatFilesStoragePath'], "SpliceGraphsFiltered", self.currentChromosome)
        if not os.path.exists(storagePath):
            os.mkdir(storagePath)
            os.chmod(storagePath, 0770)
        #if not os.path.exists(storagePathFiltered):
        #    os.mkdir(storagePathFiltered)
        #    os.chmod(storagePathFiltered, 0770)

        gnIds = self.genes.keys()
        gnIds.sort(cmp=lambda x,y: cmp(x.split(":")[1], y.split(":")[1]))
        for g in xrange(len(gnIds)):
            try:
                sg = self.__makeSpliceGraph(gnIds[g])
            except:
                self.__inform("WARNING: Could generate SpliceGraph object for gene %s: %s\n" % (gnIds[g], sys.exc_info()[1]) )
            else:
                try:
                    sg.store("%s/%s.sg" % (storagePath, gnIds[g]))
                except:
                    self.__inform("WARNING: Could not store gene %s: %s\n" % (gnIds[g], sys.exc_info()[1]) )

                
                # filtering / visualization will be called from do_complete_pipeline.py
                #try:
                #    sg.visualize("%s/%s.ps" % (storagePath, gnIds[g]))
                #except:
                #    self.__inform("WARNING: Could not visualize gene %s: %s\n" % (gnIds[g], sys.exc_info()[1]) )
                #
                ## do the same thing for the filtered splice graph
                #newsg = sg.remove_by_coverage_and_ss(coverage=self.conf[self.species]["FilterMinCoverage"])
                #
                #try:
                #    newsg.store("%s/%s.sg" % (storagePathFiltered, gnIds[g]))
                #except:
                #    self.__inform("WARNING: Could not store filtered gene %s: %s\n" % (gnIds[g], sys.exc_info()[1]) )
                #
                #try:
                #    newsg.visualize("%s/%s.ps" % (storagePathFiltered, gnIds[g]))
                #except:
                #    self.__inform("WARNING: Could not visualize filtered gene %s: %s\n" % (gnIds[g], sys.exc_info()[1]) )

            self.__inform("%i done...\r" % (g+1), 2)

        self.__inform("all done             \n", 2)


    def __makeSpliceGraph(self, gnId):
        "make SpliceGraph object from gene id and data in self"

        gnElements = self.genes[gnId].keys()

        sgConnections = {}
        sgTypes       = {}
        sgTranscripts = {}

        for gnElement1 in gnElements:
            
            if self.mSS.has_key(gnElement1):

                if not sgConnections.has_key(gnElement1):
                    sgConnections[gnElement1] = {}
                    sgTranscripts[gnElement1] = {}

                for gnElement2 in self.mSS[gnElement1].keys():
                    sgConnections[gnElement1][gnElement2] = self.mSS[gnElement1][gnElement2]
                    sgTranscripts[gnElement1][gnElement2] = self.mSS_transcript_ids[gnElement1][gnElement2]

            for type in self.elements.keys():

                if not sgTypes.has_key(type):
                    sgTypes[type] = {}

                if self.elements[type].has_key(gnElement1):
                    sgTypes[type][gnElement1] = self.elements[type][gnElement1]

        sg = SpliceGraph.SpliceGraph(name=gnId, connections=sgConnections, types=sgTypes, transcripts=sgTranscripts,
                                     species=self.species, assembly=self.conf[self.species]['Assembly'])

        return sg

# classification of standard events is now called by do_complete_pipeline.py
##     def classifyStandardEvents(self, indir, outdir):
##         "classify standard events for current genes and store as flat files"

##         self.__inform("\nclassifying standard events\n", 2)

##         # check existing outdir
##         if not os.path.exists(outdir):
##             os.mkdir(outdir)
##             os.chmod(outdir, 0770)

##         # check existing indir
##         if os.path.isdir(indir):
        
##             # get classifier methods
##             c = Classifier.Classifier()
##             classifyMethod = {}
##             for event in self.standardEvents:
##                 classifyMethod[event] = getattr(c, "classify_%s" % event)

##             # open filehandles
##             filehandle = {}
##             for event in self.standardEvents:
##                 filehandle[event] = open("%s/%s_events" % (outdir, event), "w")
##                 for string in c.prettyString(None, event):
##                     filehandle[event].write(string)

##             # get list of input files
##             infiles = [ f for f in os.listdir(indir) if f.endswith('.sg')]

##             # classify events in genes
##             for g in xrange(len(infiles)):
##                 infile = infiles[g]
##                 try:
##                     sg = SpliceGraph.SpliceGraph(filename=indir+'/'+infile)
##                 except:
##                     self.__inform('ERROR creating SpiceGraph object for %s: %s\n' % (infile, sys.exc_info()[1]))
##                 else:
##                     self.__inform("doing %i of %i (%s): " % (g+1, len(infiles), infile[0:-3]), 2)

##                     for event in self.standardEvents:
##                         self.__inform("%s " % (event), 2)
##                         events = classifyMethod[event](sg)

##                         for string in c.prettyString(events, event, infile[0:-3]):
##                             filehandle[event].write(string)

##                     self.__inform("\r", 2)

##                 for event in self.standardEvents:
##                     filehandle[event].flush()

##             self.__inform("all done             \n", 2)

##             # close filehandles, make group writable
##             for event in self.standardEvents:
##                 filehandle[event].close()
##                 os.chmod("%s/%s_events" % (outdir, event), 0660)

##         else:
##             self.__inform("classifyStandardEvents: skipping nonexisting indir '%s'\n" % indir)


    def Test_Read_Shelve(self):
        if self.__useShelve is False:
            return 
        shelveHandle = shelve.open('%s'%self.mSS_files[0], 'r')
        mSS = shelveHandle['mSS']
        keys = mSS.keys()
      
        # Sorting, based on chr and then coordinates
        #keys.sort(cmp=lambda x,y: cmp(x.split(':')[0], y.split(':')[0]) or cmp(int(x.split(':')[1]), int(y.split(':')[1])))
        # Sorting on coordinates directly
        keys.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])))

        #for i in xrange(30):
        #    print keys[i], mSS[keys[i]]
            
        elements = shelveHandle['elements']

        #self.printDictionary(elements)
        self.__inform('\nNr of elements per type in mSS:\n\n')
        TSS_elements = len(elements['TSS'].keys())
        self.__inform('Nr of TSS elements: %d\n' % TSS_elements)
        v5SS_elements = len(elements['5ss'].keys())
        self.__inform('Nr of 5ss elements: %d\n' % v5SS_elements)
        v3SS_elements = len(elements['3ss'].keys())
        self.__inform('Nr of 3ss elements: %d\n' % v3SS_elements)
        TER_elements = len(elements['TER'].keys())
        self.__inform('Nr of TER elements: %d\n' % TER_elements)

        TSS_elements = elements['TSS'].keys()
        v5ss_elements = elements['5ss'].keys()
        v3ss_elements = elements['3ss'].keys()
        TER_elements = elements['TER'].keys()

        self.__inform('\nFinding the first TSS and printing the following 29 elements and their maps:\n\n')

        TSS_elements.sort(cmp=lambda x,y: cmp(int(x.split(':')[1]), int(y.split(':')[1])))
        start_index = keys.index(TSS_elements[0])
        for i in xrange(start_index,start_index+30,):
            print keys[i], mSS[keys[i]], keys[i] in TSS_elements,  keys[i] in v5ss_elements,
            print keys[i] in v3ss_elements, keys[i] in TER_elements


    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            #print >> sys.stderr, txt,
            self.log.write(txt)
            self.log.flush()


    def printDictionary(self, dict):
        print >> sys.stdout, "\n".join(["%3s = %s" % (k, dict[k]) for k in dict.keys()])


if __name__ == "__main__":
    print 'Testing GenerateSpliceGraphs Class\n'
    test_dict = {"user": os.getlogin(),
                 "DB": 'mm6',
                 "PASSWORD":'',
                 "host": 'retro.mit.edu',
                 #"EST":['chr21_intronEst'], #these are generated automatically based on Chromosomes/ESTsuffix
                 "knownGene":'knownGene',
                 "mRNA":'all_mrna',
                 "MetaTable":'metatable',
                 "Chromosomes":['chr14'],
                 }
    


    generate_splice_graphs = GenerateSpliceGraphs(MySQL_information=test_dict, species="mmu", verbose=3, useShelve=False)
    
    print generate_splice_graphs.__doc__
    generate_splice_graphs.printMySQL_information()
    generate_splice_graphs.generate_splice_graphs()
    #generate_splice_graphs.classifyStandardEvents()
    #generate_splice_graphs.Test_Read_Shelve()
    
