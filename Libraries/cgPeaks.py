#module for analyzing stretches of expression signals.
import stepVectorScan as svs
import bioLibCG as cg
import cgConfig as c

def returnPeaks(pointsDict, span = 2, minLevel = 5):
	lowest = pointsDict.keys()
	lowest.sort()
	peaks = []
	
	for i in range(int(lowest[0]) + span + 1, int(lowest[-1]) - span - 1):
		val = pointsDict[i]
		if val < minLevel: #minimum
			continue
			
		#grab all values of density around peak
		checkList = []
		j = 1
		while j <= span:
			checkList.extend([pointsDict[i - j], pointsDict[i + j]])
			j += 1
		
		
		#check if lower than val
		peakFlag = True
		for check in checkList:
			if val >= check: #greater or equal...
				pass
			else:
				peakFlag = False
			
		#output
		if peakFlag:
			#print checkList, val
			peaks.append(i)
		
	return peaks

def returnContBlocks(profile, tcc, minLevel = 5):
	chrom, strand, start, end = cg.tccSplit(tcc)
	
	pCoords = profile.keys()
	pCoords.sort()
	
	inBlock = False
	bStart = None
	
	blocks = []	
	for pCoord in pCoords:
		if int(profile[pCoord]) > minLevel: #expression's high enough
			if inBlock:
				continue
			else: #block Start
				bStart = pCoord
				inBlock = True
		else: #not high enough
			if inBlock: #end the block
				blocks.append('%s:%s:%s:%s' % (chrom, strand, bStart, pCoord - 1))
				inBlock = False
			else:
				continue
	
	return blocks				
				
		
class stretch:
	
	def __init__(self, tcc, config):
		
		
		#get the absolute profile for this stretch
		self.tcc = tcc
		#this uses the merged wig files...
		self.profile = svs.svCoord([tcc], config)
		
	
	def createPeaks(self, span = 2, minVal = 5):
		self.peaks = returnPeaks(self.profile, span, minVal)
	
	def createContBlocks(self, minVal = 5):
		self.blocks = returnContBlocks(self.profile, self.tcc, minVal)
		
	def getHighestLevel(self, coord = False):
		highest = 0
		c = None
		for i in self.profile:
			if self.profile[i] > highest:
				highest = self.profile[i]
				c = i
		
		if coord:
			return c
		else:
			return highest
			
	def getHighestPeak(self):
		highest = None
		hVal = 0
		for pCoord in self.peaks:
			if self.profile[pCoord] > hVal:
				hVal = self.profile[pCoord]
				highest = pCoord
			
		return highest
	
	def getSumOfLevels(self):
		sum = 0
		for pCoord in self.profile:
			sum += self.profile[pCoord]
		
		return sum
	

	
	
		

		
