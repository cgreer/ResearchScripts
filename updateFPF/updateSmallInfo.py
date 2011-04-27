import bioLibCG as cg
import cgConfig

#init
conf = cgConfig.cgConfig()
mainConf = cgConfig.cgConfig('Main.conf')
pFileName = conf.conf['results']
smallFileName = conf.conf['smallAnalysis']
organism = conf.conf['organism']

#make list of acceptable files for organism
organismFileList = [] #NOTE: THESE ARE THE NAKED BASE FILE NAMES!
metaDict = cg.getMetaFileDict(mainConf.conf['metaFileName'])
for baseFName in metaDict:
	if metaDict[baseFName][1] == organism:
		organismFileList.append(baseFName)

#put small hits for each prediction in dictionary
pCount = {}
smallFile = open(smallFileName, 'r')

currID = 'NONE'
for line in smallFile:
	if '\t' not in line: #This is the line with the id in it -->  store another ID
		currID = line.strip()
	else: #this line contains library and count info  --> add 
		lib = line.strip().split('\t')[0]
		count = int(line.strip().split('\t')[1])
		
		if cg.getBaseFileName(lib, naked = True) in organismFileList:
			if currID in pCount:
				pCount[currID] = pCount[currID] + count
			else:
				pCount[currID] = count

#update the file --> any line with kmer on it give it the count
newLines = []
predFile = open(pFileName, 'r')
for line in predFile:
	kmer = line.strip().split('\t')[0].split('.')[0]
	if kmer in pCount:
		numSmall = pCount[kmer]
	else:
		numSmall = 0
	newLine = line.strip().split('\t')
	newLine[9] = str(numSmall)
	newLine = '\t'.join(newLine) + '\n'
	newLines.append(newLine)
predFile.close()

#write new File
newFile = open(pFileName, 'w')
for line in newLines:
	newFile.write(line)
newFile.close()


#################Now update cluster Count#######################

clusterCount = {}
pFile = open(pFileName, 'r')
for line in pFile:
	predictionCount = int(line.strip().split('\t')[9])
	CID = line.strip().split('\t')[7]
	if CID in clusterCount:
		clusterCount[CID] = clusterCount[CID] + predictionCount
	else:
		clusterCount[CID] = predictionCount
pFile.close()

#update the file --> cluster small count
newLines = []
predFile = open(pFileName, 'r')
for line in predFile:
	CID = line.strip().split('\t')[7]
	numSmall = clusterCount[CID]
	newLine = line.strip().split('\t')
	newLine[10] = str(numSmall)
	newLine = '\t'.join(newLine) + '\n'
	newLines.append(newLine)
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

################REORDER THE FILE BY FRAMESTART###################



