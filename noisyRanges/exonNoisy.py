#for each exon prediction, get noise levels of all other transcripts in genes it overlaps.
import bioLibCG as cg
from bioLibCG import ss
import cgConfig as c
import getHairpins
import wigValue
import compareData as compare
import cgGenes
import stepVectorScan

def exonNoisy(cName = None):
	#init
	mConf = c.cgConfig('Main.conf')
	conf = c.getConfig(cName)
	cHairs = getHairpins.getHairpins(conf.conf['resultsExons']) #CID: HAIRPIN
	organism = conf.conf['organism']
	geneSetFolder = mConf.conf['geneSets%s' % organism]
	
	#make prediction overlap hitmap
	print 'Making prediction list'
	predList = [] 
	for CID in cHairs:
		hPin = cHairs[CID]
		predList.append(hPin)
	
	if compare.checkIfOverlaps(predList):
		predList = compare.collapseOverlaps(predList)
	
	
	#make genes for Ensemble/make list of tccs for exons.
	print 'Creating gene set'
	
	ensGenes = cgGenes.createGeneSetFromFile(geneSetFolder + '/ensemblAllExons.tsv')
	print '  loaded # genes:', len(ensGenes.set)
	
	
	#collect levels for each haipin region
	print '[Checking all levels]'
	cidLevels = {}
	for CID in cHairs:
		print CID
		hPin = cHairs[CID]
			
		#for each hairpin, --> find overlapping transcripts in same gene
		overlappingGenes = ensGenes.geneOverlaps([hPin])
		if len(overlappingGenes) > 0:
			gIDs = [gene.id for gene in overlappingGenes]
			allTccs = ensGenes.getTccsFromGIDs(gIDs)
			if compare.checkIfOverlaps:
				print '  Overlaps...collapsing'
				allTccs = compare.collapseOverlaps(allTccs)
		else:
			print 'NO GENE OVERLAPS!!!!!', CID, hPin
		
		
		#filter out my predictions.
		print '  Filtering out predictions'
		checkList = compare.subtractTwoTccLists(allTccs, predList)
			
		
		#Get Expression level for gene.
		print '  Retrieving Expression levels:', cg.getTccListTotalLength(checkList)
		levels = []
		
		
		hPinLevels = stepVectorScan.scanVectorsHist(checkList, cName)
		for hPin in hPinLevels:
			levels.extend(hPinLevels[hPin])
		
			
		cidLevels[CID] = levels
		
	
	
	
	#output levels to file
	print 'Outputting to file'
	#find longest
	longest = 0
	for CID in cidLevels:
		length = len(cidLevels[CID])
		if length > longest:
			longest = length
	
	sortedKeys = cidLevels.keys()
	sortedKeys.sort()
	#print sortedKeys
	
	newLines = []
	for j in range(0, longest): #how many lines are there
		newLine = []
		for CID in sortedKeys:
			if len(cidLevels[CID]) > j:# add it
				newLine.append(str(cidLevels[CID][j]))
			else:
				newLine.append('NA')
	
		newLines.append('\t'.join(newLine) + '\n')
	
	outFileN = conf.conf['exonNoiseData']
	outFile = open(outFileN, 'w')
	outFile.write('\t'.join(sortedKeys) + '\n')
	outFile.writelines(newLines)
	outFile.close()
		
if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		exonNoisy(sys.argv[1])
	else:
		exonNoisy()

	
	
	
	
	
