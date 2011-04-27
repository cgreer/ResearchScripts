##This is for taking Ensembl Transcipt data and making gene structure with them
import bioLibCG as cg
import compareData as compare


class transcript:
	'''upgrade to exons and introns...'''
	def __init__(self, tcc, id, parent, type = None):
		self.tcc = str(tcc)
		self.id = str(id)
		self.parent = str(parent)
		self.type = type

class gene:
	
	def __init__(self, id, transcriptList, gType):
		self.id = id
		self.transcripts = transcriptList
		self.type = gType
		

class geneSet:
	
	def __init__(self, geneList):
		self.set = {}
		
		for gene in geneList:
			self.set[gene.id] = gene #id: [ts, etc]
		
		self.genes = self.set.values()
		self.gIDs = self.set.keys()
		
	
	#given a tcc or tccs, which transcripts overlap with it in this gene set?
	def transcriptOverlaps(self, tccs):
		'''return list of overlapping transcripts'''
		if not isinstance(tccs, type([])):
			print 'transcript overlaps: NEED TCC LIST, not a single tcc!'
			return 1
		
		print 'num of tccs being compared', len(tccs)
		#gather all transcript tccs(make tcc --> id).
		tccDict = {}
		for gene in self.genes:
			for transcript in gene.transcripts:
				if transcript.tcc in tccDict:
					tccDict[transcript.tcc].append(transcript.id)
				else:
					tccDict[transcript.tcc] = [transcript.id]
					
		print 'num of tcc transcript tccs', len(tccDict.keys())
		
		
		overlapped = compare.compareTwoTcc(tccDict.keys(), tccs, 1)
		
		tList = []
		for gene in self.genes:
			for transcript in gene.transcripts:
					if transcript.tcc in overlapped:
						tList.append(transcript)
		
		return tList
	
	def getGene(self, gID):
		return self.set[gID]
		
	def getTccsFromGIDs(self, idList):
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
	
	def getParentGenes(self, transcriptList):
		'''Given a list of transcripts, give a list of all the genes that those transcripts belong to'''
		
		myGIDs = []
		for transcript in transcriptList:
			if transcript.parent not in myGIDs:
				myGIDs.append(transcript.parent)
		
		geneList = []
		for gID in myGIDs:
			geneList.append(self.getGene(gID))
		
		return geneList
	
	def geneOverlaps(self, tccs):
		overlappingTranscripts = self.transcriptOverlaps(tccs)
		print 'num of overlapping transcripts', len(overlappingTranscripts)
		overlappingGenes = self.getParentGenes(overlappingTranscripts)
		print 'num of overlapping genes', len(overlappingGenes)
		
		return overlappingGenes
				
			
def createGeneSetFromFile(fN, type = 'Ensembl'):
	'''Note that the files must be in a specific format!!!'''
	
	file = open(fN, 'r')
	file.readline() #header
	
	#collect all transcripts
	allTranscripts = []
	for line in file:
		s = line.strip().split('\t')
		try:
			gID, tID, chr, strand, start, end, tType = s[0], s[1], s[2], s[3], s[4], s[5], s[6]
		except IndexError:
			gID, tID, chr, strand, start, end, tType = s[0], s[1], s[2], s[3], s[4], s[5], "None"
			
		tcc = 'chr%s:%s:%s:%s' % (chr, strand, start, end)
		t = transcript(tcc, tID, gID, tType)
		allTranscripts.append(t)
	print 'testing transcript creation', allTranscripts[0].tcc, allTranscripts[0].parent, allTranscripts[0].id, allTranscripts[0].type
		
	#put transcripts into genes
	genes = {}
	for t in allTranscripts:
		if t.parent in genes:
			genes[t.parent].append(t)
		else:
			genes[t.parent] = [t]
		
	allGenes = []
	for gID in genes:
		#get gene type
		for t in genes[gID]:
			gType = t.type
			break
		
		
		g = gene(gID, genes[gID], gType)
		allGenes.append(g)
	
	#print 'testing gene creation', allGenes[0].transcripts, allGenes[0].id
	
	#return a gene set
	return geneSet(allGenes)
	
		
		
	
	
	
	
		
	
