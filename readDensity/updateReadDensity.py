###This will calculate the total number of hits per microRNA --> update 
##the newResults file with them --> put the value of the highest mature --> CID

import cgConfig
import bioLibCG as cg

#go through wig each chromosome and check the mature seqs
mainConf = cgConfig.cgConfig('Main.conf')
conf = cgConfig.cgConfig()
newLines = []
pFileName = 'newResults.results'


for wigFileN in cg.recurseDir(mainConf.conf['wigMouse'], end = '.wig'):
	
	
	#init
	chrom = wigFileN.strip().split('.')[-2]
	strand = wigFileN.strip().split('.')[-4]
	wigFile = open(wigFileN, 'r')
	mirFile = open(conf.conf['results'], 'r')
	print wigFileN
	
	#get rid of header
	wigFile.readline()
	
	print '  populating hitmap'
	#populate hitmap
	wigMap = {}
	for line in wigFile:
		value = int(line.strip().split('\t')[3].split('.')[0])
		if value > 0:
			start = int(line.strip().split('\t')[1])
			end = int(line.strip().split('\t')[2])
			for i in range(start, end):
				wigMap[i] = value
	wigFile.close()
	
	print '  calculating hits for mature seqs'
	#calculate total hits per mature
	for line in mirFile:
		mTcc = line.strip().split('\t')[1]
		mirID = line.strip().split('\t')[0]
		if (mTcc.split(':')[0] == chrom) and (mTcc.split(':')[1] == strand):
			#if mirID == '26477.30.106643972': print 'Starting Total Count'
			totalHits = 0
			for i in range(int(mTcc.split(':')[2]), int(mTcc.split(':')[3])):
				#if mirID == '26477.30.106643972': print '  ', i 
				if i in wigMap:
					totalHits += wigMap[i]
					#if mirID == '26477.30.106643972': print '    ', i, totalHits, wigMap[i]
		
			newLines.append(cg.appendToLine(line, str(totalHits), 11))
	
	mirFile.close()

print 'Writing New File'
#write new results file
outFile = open(pFileName, 'w')
for line in newLines:
	outFile.write(line)
outFile.close()

####NOW UPDATE MAX HITS PER CLUSTER####

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



		
			
				
		
	
	
				
