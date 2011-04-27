##This is for taking Ensembl Transcipt data and making gene structure with them
import bioLibCG
import compareData as compare
import GenomeFetch
import cgGenes3
import cgSeqMod

class transcript:
	'''Remember, exons does NOT mean Protein.  The UTRs can and WILL overlap with exons'''
	def __init__(self):
		self.tcc = ''
		self.id = ''
                self.strand = ''
                self.chromosome = ''
                self.utr5 = []
                self.exonList = [] #[ [], [], ]
                self.intronList = []
                self.utr3 = []
		self.parent = '' #ID of the parent...
                self.cds5Stat = ''
                self.cds3Stat = ''
                self.tType = None
                self.gType = None

        def getOverlappingElements(self, tcc):
                '''Given region, Which element (INTRON, EXON, 5UTR, 3UTR)'''
                overlappingElements = []
                try:
                        for utrSegment in self.utr5: 
                                utr5Tcc = bioLibCG.makeTcc(self.chromosome, self.strand, utrSegment[0], utrSegment[1])
                                if bioLibCG.tccOverlap(utr5Tcc, tcc):
                                        overlappingElements.append([utrSegment, '5UTR'])
                except IndexError:
                        pass

                for exon in self.exonList:
                        exonTcc = bioLibCG.makeTcc(self.chromosome, self.strand, exon[0], exon[1])
                        #print '@ ', exonTcc, tcc, 'EXON'
                        if bioLibCG.tccOverlap(exonTcc, tcc):
                                overlappingElements.append([exon, 'EXON'])

                for intron in self.intronList:
                        intronTcc = bioLibCG.makeTcc(self.chromosome, self.strand, intron[0], intron[1])
                        #print '@ ', intronTcc, tcc, 'INTRON'
                        if bioLibCG.tccOverlap(intronTcc, tcc):
                                overlappingElements.append([intron, 'INTRON'])

                try:
                        for utrSegment in self.utr3: 
                                utr3Tcc = bioLibCG.makeTcc(self.chromosome, self.strand, utrSegment[0], utrSegment[1])
                                if bioLibCG.tccOverlap(utr3Tcc, tcc):
                                        overlappingElements.append([utrSegment, '3UTR'])
                except IndexError:
                        pass

                #!!!Eventually add a way to find if overlapping EXON_UTR as well
                if 'EXON' in overlappingElements and 'INTRON' in overlappingElements:
                        overlappingElements.append('EXON_INTRON')
                        overlappingElements.remove('EXON')
                        overlappingElements.remove('INTRON')


                return overlappingElements

        def getRelativePositionMRNA(self, position, coding = True):
                coding = bool(coding)

                if self.strand == '1':
                        
                        try:
                                cdsStart = self.utr5[-1][1] + 1 #end of utr + 1
                        except IndexError:
                                cdsStart = int(self.tcc.split(':')[2])
                        
                        try:
                                cdsEnd = self.utr3[-1][0] - 1 #Check if this is how you added Elements to the UTR
                        except IndexError:
                                cdsEnd = int(self.tcc.split(':')[3])

                        if not coding:
                                cdsStart = int(self.tcc.split(':')[2])
                                cdsEnd = int(self.tcc.split(':')[3])
                        
                        theTcc = bioLibCG.makeTcc(self.chromosome, self.strand, position, position + 1)
                        overlapTypes = self.getOverlappingElements(theTcc)
                        
                        if not len(overlapTypes) == 1:
                                return -1

                        if not overlapTypes[0][1] == 'EXON':
                                return -1

                        coordinates = []
                        coordinates.append(position)
                        coordinates.append(cdsStart)

                        #Collect all the coords in between cdsStart and position
                        for exon in self.exonList:
                                for exonCoord in exon:
                                        if cdsStart < exonCoord < position:
                                                coordinates.append(exonCoord)

                        #now start at cdsStart and make pairs for each region till cdsEnd
                        coordinates.sort()
                        #print coordinates
                        starts = [(coordinates[i*2]) for i in range(len(coordinates)//2)] 
                        ends = [coordinates[1 + i*2] for i in range(len(coordinates)//2)]
                        pairs = zip(starts, ends)
                        #print pairs

                        totalLength = 0
                        for pair in pairs:
                                #print pair[0], pair[1], pair[1] - pair[0] + 1
                                totalLength += pair[1] - pair[0] + 1

                        #print totalLength

                        #return (-1) because we are looking for its position (0 based), not its length
                        return totalLength - 1

                if self.strand == '-1':
                        
                        try:
                                cdsStart = self.utr5[-1][0] - 1
                        except IndexError:
                                cdsStart = int(self.tcc.split(':')[3])
                        
                        try:
                                cdsEnd = self.utr3[-1][1] + 1
                        except IndexError:
                                cdsEnd = int(self.tcc.split(':')[2])

                        if not coding:
                                cdsStart = int(self.tcc.split(':')[3])
                                cdsEnd = int(self.tcc.split(':')[2])

                        theTcc = bioLibCG.makeTcc(self.chromosome, self.strand, position, position + 1)
                        overlapTypes = self.getOverlappingElements(theTcc)
                        
                        if not len(overlapTypes) == 1:
                                return -1

                        if not overlapTypes[0][1] == 'EXON':
                                return -1
                        
                        coordinates = []
                        coordinates.append(position)
                        coordinates.append(cdsStart)
                        for exon in self.exonList:
                                for exonCoord in exon:
                                        if cdsStart > exonCoord > position:
                                                coordinates.append(exonCoord)

                        #now start at cdsStart and make pairs for each region till cdsEnd
                        coordinates.sort()
                        starts = [(coordinates[i*2])for i in range(len(coordinates)//2)] 
                        ends = [coordinates[1 + i*2] for i in range(len(coordinates)//2)]
                        pairs = zip(starts, ends)

                        totalLength = 0
                        for pair in pairs:
                                totalLength += pair[1] - pair[0] + 1 #0-based.  len = y - x + 1

                        return totalLength - 1



        def getMRNA(self, coding = False):

                gf = GenomeFetch.GenomeFetch('hg19')                
                if self.strand == '1':
                        mRNAList = []
                       
                        #Get CDS START/END
                        try:
                                cdsStart = self.utr5[-1][1] + 1 #start of coding is right after utr
                        except IndexError:
                                cdsStart = int(bioLibCG.tccSplit(self.tcc)[2])
                        try:
                                cdsEnd = self.utr3[-1][0] - 1 #Check if this is how you added Elements to the UTR
                        except IndexError:
                                cdsEnd = int(bioLibCG.tccSplit(self.tcc)[3])
                        if not coding:
                                chrom, strand, start, end = bioLibCG.tccSplit(self.tcc)
                                cdsStart = start
                                cdsEnd = end
                        
                        coordinates = []
                        coordinates.extend([cdsStart, cdsEnd])

                        for exon in self.exonList:
                                for exonCoord in exon:
                                        if cdsStart < exonCoord < cdsEnd:
                                                coordinates.append(exonCoord)

                        #now start at cdsStart and make exon pair for each region till cdsEnd
                        coordinates.sort()
                        tccStarts = [(coordinates[i*2]) for i in range(len(coordinates)//2)]
                        tccEnds = [coordinates[1 + i*2] for i in range(len(coordinates)//2)]
                        tccPairs = zip(tccStarts, tccEnds)

                        for tccPair in tccPairs:
                                mRNAList.append(gf.get_seq_from_to(self.chromosome, tccPair[0] + 1, tccPair[1] + 1)) #gf is 1-based...

                        mRNA = ''.join(mRNAList)
                        return mRNA.replace('T', 'U')

                if self.strand == '-1':
                        '''This this differs in that: 1) the UTR coordinates have to be retrieved differently, 2) The 
                        exon coordinates that will be kept have to be checked with signs reversed, 3) have to revComp it'''

                        mRNAList = []
                        
                        try:
                                cdsStart = self.utr5[-1][0] - 1
                        except IndexError:
                                cdsStart = int(bioLibCG.tccSplit(self.tcc)[3])
                        try:                                
                                cdsEnd = self.utr3[-1][1] + 1 #Check if this is how you added Elements to the UTR
                        except IndexError:
                                cdsEnd = int(bioLibCG.tccSplit(self.tcc)[2])
                        if not coding:
                                chrom, strand, start, end = bioLibCG.tccSplit(self.tcc)
                                cdsStart = end
                                cdsEnd = start
                        coordinates = []
                        coordinates.extend([cdsStart, cdsEnd])

                        for exon in self.exonList:
                                for exonCoord in exon:
                                        if cdsStart > exonCoord > cdsEnd: #switched signs for negative
                                                coordinates.append(exonCoord)

                        #now start at cdsStart and make tccs for each region till cdsEnd
                        coordinates.sort()

                        #print coordinates
                        #print len(coordinates)
                        tccStarts = [(coordinates[i*2]) for i in range(len(coordinates)//2)] 
                        tccEnds = [coordinates[1 + i*2] for i in range(len(coordinates)//2)]
                        tccPairs = zip(tccStarts, tccEnds)

                        #print tccPairs
                        for tccPair in tccPairs:
                                mRNAList.append(gf.get_seq_from_to(self.chromosome, tccPair[0] + 1, tccPair[1] + 1)) #gf is 1-based

                        mRNA = ''.join(mRNAList)
                        mRNA = cgSeqMod.reverseComplementSequence(mRNA, False)
                        return mRNA.replace('T', 'U')
                        

class gene:
	
	def __init__(self, id, transcriptList):
		self.id = id
		self.transcripts = transcriptList
		

class geneSet:
	
	def __init__(self, geneList):
		self.set = {}
		
		for gene in geneList:
			self.set[gene.id] = gene #id: [ts, etc]
		
		self.genes = self.set.values()
		self.gIDs = self.set.keys()
	        
                tList = []
                for gene in self.genes:
                        for transcript in gene.transcripts:
                                tList.append(transcript)

                self.transcripts = tList
	
	def getOverlappingTranscriptList(self, tccs):
		'''given tccs, return the transcripts that overlap with them in a single list'''
		if not isinstance(tccs, type([])):
			#print 'transcript overlaps: NEED TCC LIST, not a single tcc!'
			return 1
		
		#print 'num of tccs being compared', len(tccs)
		#gather all transcript tccs(make tcc --> id).
		tccDict = {}
		for gene in self.genes:
			for transcript in gene.transcripts:
		                tccDict[transcript.tcc] = 1

		#print 'num of tcc transcript tccs', len(tccDict.keys())
		
		
		overlapped = compare.compareTwoTcc(tccDict.keys(), tccs, 1)
		
		tList = []
		for gene in self.genes:
			for transcript in gene.transcripts:
					if transcript.tcc in overlapped:
						tList.append(transcript)
		
		return tList
	
	def getGene(self, gID):
		return self.set[gID]
		
	def getTranscriptTccsFromGIDs(self, idList):
		'''give a tcc list from a transcript list THIS GIVES ALL TRANSCRIPT IDS...'''
		
		tccList = []
		for gID in self.gIDs:
			if gID in idList:
				for transcript in self.getGene(gID).transcripts:
					tccList.append(transcript.tcc)
		
		return tccList
	
	def getTccsFromTIDs(self, idList):
		'''give a tcc list from a transcript list'''
		
		tccList = []
		
		for gene in self.genes:
			for transcript in gene.transcripts:
				if transcript.id in idList:
					tccList.append(transcript.tcc)
		
		return tccList
	
	def getParentGenesFromTranscripts(self, transcriptList):
		'''Given a list of transcripts, give a list of all the genes that those transcripts belong to'''
		
		myGIDs = []
		for transcript in transcriptList:
			if transcript.parent not in myGIDs:
				myGIDs.append(transcript.parent)
		
		geneList = []
		for gID in myGIDs:
			geneList.append(self.getGene(gID))
		
		return geneList
	
	def getOverlappingGeneList(self, tccs):
                '''Given tccs, find the genes that overlap with these tccs'''
		overlappingTranscripts = self.getOverlappingTranscriptList(tccs)
		#print 'num of overlapping transcripts', len(overlappingTranscripts)
		overlappingGenes = self.getParentGenesFromTranscripts(overlappingTranscripts)
		#print 'num of overlapping genes', len(overlappingGenes)
		
		return overlappingGenes
				
def createGeneSetEditing(fN):
        '''Read transcript info into a gene set from our editing file...'''	
	file = open(fN, 'r')
	file.readline() #header
	
	#collect all transcripts
	allTranscripts = []
	for line in file:
		
                #parse info
                ls = line.strip().split('\t')
		tID = ls[0]
                chrom = ls[1]
                strand = ls[2]
                if strand == '+':
                        strand = '1'
                else:
                        strand = '-1'
                tStart = int(ls[3])
                tEnd = int(ls[4]) - 1 #0 based
                cStart = int(ls[5])
                cEnd = int(ls[6]) - 1 #0 based
                numExons = int(ls[7])
                exonStarts = [int(x) for x in ls[8][:-1].split(',')]
                exonEnds = [int(x) - 1 for x in ls[9][:-1].split(',')] # 0 based
                geneName = ls[10]
                cStartStat = ls[11]
                cEndStat = ls[12]
                tType = ls[13]
                gType = ls[15]

                #make transcript
                t = transcript()
                t.id = tID
                t.parent = geneName
                t.strand = strand
                t.chromosome = chrom
                t.tcc = '%s:%s:%s:%s' % (chrom, strand, tStart, tEnd)  
                t.tType = tType
                t.gType = gType
                
                #Exons
                for i, eStart in enumerate(exonStarts):
                        t.exonList.append([eStart, exonEnds[i]])
                #Introns
                for i, eStart in enumerate(exonStarts):
                        if i == 0: continue
                        t.intronList.append([exonEnds[i - 1] + 1, eStart - 1])
                
                #Set UTR stats.
                if strand == '1':
                        if cStart == tStart or cStart == tEnd: 
                                t.cds5Stat = 'INC'
                        else:
                                t.cds5Stat = 'COMP'

                        if cEnd == tEnd or cEnd == tStart:
                                t.cds3Stat = 'INC'
                        else:
                                t.cds3Stat = 'COMP'
                else:
                        if cStart == tStart or cStart == tEnd: 
                                t.cds3Stat = 'INC'
                        else:
                                t.cds3Stat = 'COMP'

                        if cEnd == tEnd or cEnd == tStart:
                                t.cds5Stat = 'INC'
                        else:
                                t.cds5Stat = 'COMP'

                #only do UTR if complete!!! Will be empty list if not complete
                if (cStart != tStart): #update 3' or 5' depending on which strand...
                        if strand == '1':
                                utr = t.utr5
                        else:
                                utr = t.utr3

                        for exon in t.exonList:
                                if bioLibCG.simpleOverlap(exon[0], exon[1], cStart, cStart + 1):
                                        utr.append([exon[0], cStart - 1]) #-1 because utr ends right before coding begins...
                                        break
                                else:
                                        utr.append(exon)
                                
                
                if (cEnd != tEnd): #again, based on strand
                        if strand == '1':
                                utr = t.utr3
                        else:
                                utr = t.utr5
                        
                        for exon in reversed(t.exonList): #don't permanently reverse the list...
                                if bioLibCG.simpleOverlap(exon[0], exon[1], cEnd, cEnd + 1):
                                        utr.append([cEnd + 1, exon[1]])
                                        break
                                else:
                                        utr.append(exon)

                
                
		allTranscripts.append(t)
		
	#put transcripts into genes
	genes = {}
	for t in allTranscripts:
		if t.parent in genes:
			genes[t.parent].append(t)
		else:
			genes[t.parent] = [t]
		
	allGenes = []
	for gID in genes:
		g = gene(gID, genes[gID])
		allGenes.append(g)
	
	##print 'testing gene creation', allGenes[0].transcripts, allGenes[0].id
	
	#return a gene set
	return geneSet(allGenes)

def createGeneSetFromFile(fN, type = 'Ensembl'):
	'''Note that the files must be in a specific format!!!'''
	
	file = open(fN, 'r')
	file.readline() #header
	
	#collect all transcripts
	allTranscripts = []
	for line in file:
		s = line.strip().split('\t')
		gID, tID, chr, strand, start, end = s[0], s[1], s[2], s[3], s[4], s[5]
		tcc = 'chr%s:%s:%s:%s' % (chr, strand, start, end)
		t = transcript(tcc, tID, gID)
		allTranscripts.append(t)
	##print 'testing transcript creation', allTranscripts[0].tcc, allTranscripts[0].parent, allTranscripts[0].id
		
	#put transcripts into genes
	genes = {}
	for t in allTranscripts:
		if t.parent in genes:
			genes[t.parent].append(t)
		else:
			genes[t.parent] = [t]
		
	allGenes = []
	for gID in genes:
		g = gene(gID, genes[gID])
		allGenes.append(g)
	
	##print 'testing gene creation', allGenes[0].transcripts, allGenes[0].id
	
	#return a gene set
	return geneSet(allGenes)
	
