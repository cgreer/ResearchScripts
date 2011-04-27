#Sort the updated cluster file to display the best results
from operator import itemgetter
import cgConfig

def sortResults(cName = None):
	#INIT
	conf = cgConfig.getConfig(cName)
	
	pFileName = conf.conf['results']
	minDensity = 4
	minSmallHits = 1
	
	pFile = open(pFileName, 'r')
	fileLines = [] #This will hold the lists to be sorted...
	for line in pFile:
		fileLines.append(line.strip().split('\t'))		
	pFile.close()
	
	#highest prediction density
	densityDict = {} #CID: highest density --> used to sort out clusters without proper density
	for line in fileLines:
		CID = line[7]
		pDensity = int(line[5])
		
		if CID in densityDict:
			if pDensity > densityDict[CID]:
				densityDict[CID] = pDensity
		else:
			densityDict[CID] = pDensity
	
	
	#take out clusters that didn't make the cut
	CIDpassed = []
	keptLines = []
	for line in fileLines:
		CID = line[7]
		#smallClusterHits = int(line[10]) Not using this metric anymore...
		if (densityDict[CID] >= minDensity): #Density 
			if line[8] == '0': #doesn't overlap with anything known (the cluster doesn't...that is)
				keptLines.append(line)
				if CID not in CIDpassed:
					CIDpassed.append(CID)
	
	#remake keptLines with integer in field ten
	#at this point just sort by cluster density...
	sID = 5 
	for line in keptLines:
		line[sID] = int(line[sID])
		
	sortedData = sorted(keptLines, key=itemgetter(sID), reverse = True) #sort by small RNA hits
	
	for line in sortedData:
		line[sID] = str(line[sID])
	
	#output
	sortedFile = open(conf.conf['results'] + '.sorted', 'w')
	for line in sortedData:
		sortedFile.write('\t'.join(line) + '\n')
	sortedFile.close()
	
	#Now output stats
	statFile = open('statFile.data', 'w')
	statFile.write('Total Clusters: %s\n' % len(densityDict))
	statFile.write('Passed: %s\n' % len(CIDpassed))

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		sortResults(sys.argv[1])
	else:
		sortResults()
