###This will calculate the MAX of hits per microRNA --> update 
##the newResults file with them --> put the value of the highest mature --> CID

import cgConfig
import bioLibCG as cg
import cgPeaks

def updateReadDensity(tType, cName):
	#go through wig each chromosome and check the mature seqs
	mainConf = cgConfig.cgConfig('Main.conf')
	conf = cgConfig.getConfig(cName)
	organism = conf.conf['organism']
	wigFolder = mainConf.conf['wig%s' % organism]	
	newLines = []
	
	
	#Differentiate between exon or intron...
	if tType == 'E':
		pFileName = conf.conf['resultsExons']
	elif tType == 'I':
		pFileName = conf.conf['resultsIntrons']
	else:
		print 'READ UPDATE FAIL'

	print '  Updating Read Density:', tType

	
	#get read density for each line...
	print '  calculating hits for mature seqs'
	#calculate total hits per mature
	mirFile = open(pFileName, 'r')
	for line in mirFile:
		mTcc = line.strip().split('\t')[1]
		mirID = line.strip().split('\t')[0]
		
		tccStretch = cgPeaks.stretch(mTcc, cName)
		highestHit = 0
		for i in range(int(mTcc.split(':')[2]), int(mTcc.split(':')[3])):
			if i in tccStretch.profile:
				if tccStretch.profile[i] > highestHit:
					highestHit = tccStretch.profile[i]		
		
		newLines.append(cg.appendToLine(line, str(highestHit), 11))
	
	mirFile.close()

	print 'Writing New File'
	#write new results file
	outFile = open(pFileName, 'w')
	for line in newLines:
		outFile.write(line)
	outFile.close()

	####NOW UPDATE HIGHEST HIT PER CLUSTER####

	clusterCount = {}

	pFile = open(pFileName, 'r')
	for line in pFile:
		predictionCount = int(line.strip().split('\t')[11])
		CID = line.strip().split('\t')[7]
		if CID in clusterCount:
			if clusterCount[CID] < predictionCount:
				clusterCount[CID] = predictionCount
		else:
			clusterCount[CID] = predictionCount
	pFile.close()

	#update the file --> cluster small count
	newLines = []
	predFile = open(pFileName, 'r')
	for line in predFile:
		CID = line.strip().split('\t')[7]
		numMax = clusterCount[CID]
		newLines.append(cg.appendToLine(line, str(numMax), 12))
	predFile.close()

	#sort newLines by clusterID
	sortDict = {}
	CIDs = []
	for line in newLines:
		CID = int(line.strip().split('\t')[7])
		if CID not in CIDs:
			CIDs.append(CID)
		if CID in sortDict:
			sortDict[CID].append(line)
		else:
			sortDict[CID] = [line]
		
	CIDs.sort()

	newLines = []
	for CID in CIDs:
		for line in sortDict[CID]:
			newLines.append(line)

	#write new File
	newFile = open(pFileName, 'w')
	for line in newLines:
		newFile.write(line)
	newFile.close()

if __name__ == "__main__":
	import sys
	
	if len(sys.argv) > 2:
		updateReadDensity(sys.argv[1], sys.argv[2])
	else:
		updateReadDensity(sys.argv[1])

		
			
				
		
	
	
				
