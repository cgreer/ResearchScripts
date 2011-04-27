#for each intron prediction, get noise distribution
import bioLibCG as cg
from bioLibCG import ss
import cgConfig as c
import getHairpins
import wigValue
import compareData as compare
import stepVectorScan

def intronNoisy(cName = None):
	mConf = c.cgConfig('Main.conf')
	conf = c.getConfig(cName)
	
	#init
	cHairs = getHairpins.getHairpins(conf.conf['resultsIntrons']) #CID: HAIRPIN
	organism = conf.conf['organism']
	exonList = compare.tccFileToList('%sExons.tcc' % organism, 0)
	slide = 1000
	
	#make prediction overlap hitmap
	predMap = {}
	predList = []
	for CID in cHairs:
		hPin = cHairs[CID]
		predList.append(hPin)
	
	#collapse Overlaps
	print ' collapsing predictions'
	predList = compare.collapseOverlaps(predList)
	print ' collapsing exons'
	exonList = compare.collapseOverlaps(exonList)
	
	
	#collect levels for each hairpin region
	cidLevels = {}
	for CID in cHairs:
		print CID
		hPin = cHairs[CID]
		chrom = ss(hPin, ':')[0]
		strand = ss(hPin, ':')[1]
		start = int(ss(hPin, ':')[2])
		end = int(ss(hPin, ':')[3])
		
		scanStart = start - slide
		scanEnd = end + slide
		
		scanRange = []
		scanRange.append('%s:%s:%s:%s' % (chrom, strand, scanStart, start))
		scanRange.append('%s:%s:%s:%s' % (chrom, strand, end, scanEnd))
		
		print scanRange
		scanRange = compare.subtractTwoTccLists(scanRange, predList)
		scanRange = compare.subtractTwoTccLists(scanRange, exonList)
			
		levels = []
		
		print '  Retrieving Expression levels:', cg.getTccListTotalLength(scanRange)
		levels = []
		
		
		hPinLevels = stepVectorScan.scanVectorsHist(scanRange, cName)
		for hPin in hPinLevels:
			levels.extend(hPinLevels[hPin])
		
			
		cidLevels[CID] = levels
		
	#output levels to file
	
	#find longest
	longest = 0
	for CID in cidLevels:
		length = len(cidLevels[CID])
		if length > longest:
			longest = length
	
	sortedKeys = cidLevels.keys()
	sortedKeys.sort()
	
	newLines = []
	for j in range(0, longest): #how many lines are there
		newLine = []
		for CID in sortedKeys:
			if len(cidLevels[CID]) > j:# add it
				newLine.append(str(cidLevels[CID][j]))
			else:
				newLine.append('NA')
	
		newLines.append('\t'.join(newLine) + '\n')
	
	outFileN = conf.conf['intronNoiseData']
	outFile = open(outFileN, 'w')
	outFile.write('\t'.join(sortedKeys) + '\n')
	outFile.writelines(newLines)
	outFile.close()
		
if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		intronNoisy(sys.argv[1])
	else:
		intronNoisy()
	
	
	
	
	
