#for cid, check if passed pval and update it.
import cgStats
import cgConfig as c
import bioLibCG as cg


def updateNoise(pType, cName = None):
	
	#init
	mainConf = c.cgConfig('Main.conf')
	conf = c.getConfig(cName)

	if pType == 'E':
		predName = conf.conf['resultsExons']
	else:
		predName = conf.conf['resultsIntrons']


	#populate cid: exon dist
	print 'Populating CID/INtron/exon distribution data'
	if pType == 'E':
		noiseFN = conf.conf['exonNoiseData']
		f = open(noiseFN, 'r')
	else:
		noiseFN = conf.conf['intronNoiseData']
		f = open(noiseFN, 'r')
		
	exonDists = {} #cid: [exon dist]
	header = f.readline()
	order = {} # num:CID
	for i,CID in enumerate(header.strip().split('\t')):
		order[i] = CID
		exonDists[CID] = []
		
	for line in f:
		data = line.strip().split('\t')
		for i, dataPoint in enumerate(data):
			if dataPoint == 'NA' or dataPoint == '':
				continue
			else:
				dataPoint = float(dataPoint)
				CID = order[i]
				exonDists[CID].append(dataPoint)

	#get highest expression level for each cluster
	print 'Populating highest expression levels'
	predExpression = {} # CID; highest level
	exonFile = open(predName, 'r')
	for line in exonFile:
		CID = line.strip().split('\t')[7]
		hDensity = line.strip().split('\t')[12]
		
		predExpression[CID] = hDensity



	#get pVals for each CID
	print 'Getting pvals for each cluster'
	pVals = {} # CID; [lam,pVal]
	for CID in exonDists:
		if not len(exonDists[CID]) > 0: #no data in 2kb range.
			lam = 'NA'
			pVal = 'NA'
		else:
			lam = cgStats.getLam(exonDists[CID])
			pVal = cgStats.getPValExp(predExpression[CID], lam)
			
		pVals[CID] = [lam, pVal] #lam gives a good approximation of noise levels in region...

	print 'Updating the file'
	#update file...
	predFile = open(predName, 'r')
	newLines = []
	for line in predFile:
		CID = line.split('\t')[7]
		newLine = cg.appendToLine(line, pVals[CID][0], 14)
		newLine = cg.appendToLine(newLine, pVals[CID][1], 15)
		newLines.append(newLine)
	predFile.close()

	predFile = open(predName, 'w')
	predFile.writelines(newLines)
	predFile.close()
	

if __name__ == "__main__":
	import sys
	
	if len(sys.argv) > 2:
		updateNoise(sys.argv[1], sys.argv[2])
	else:
		updateNoise(sys.argv[1])	

	
		
	
			
