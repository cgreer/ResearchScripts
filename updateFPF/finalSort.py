#Sort the updated cluster file to display the best results
from operator import itemgetter
import cgConfig

	
def finalSort(pType, cName = None):
	
	#INIT
	conf = cgConfig.getConfig(cName)

	if pType == "E":
		pFileName = conf.conf['resultsExons']
	else:
		pFileName = conf.conf['resultsIntrons']
		
	minDensity = 4
	maxPVal = float(.05)

	pFile = open(pFileName, 'r')
	fileLines = [] #This will hold the lists to be sorted...
	for line in pFile:
		fileLines.append(line.strip().split('\t'))		
	pFile.close()


	keptLines = []
	for line in fileLines:
		CID = line[7]
		cDensity = int(line[12])
		pVal = float(line[15])
		
		#print pVal, cDensity
		if cDensity > 5 and pVal < maxPVal:
			#print '  kept'
			keptLines.append(line)

	print len(keptLines)
	#remake keptLines with float(pVal)
	i = 12
	for line in keptLines:
		line[i] = float(line[i])
		
	sortedData = sorted(keptLines, key=itemgetter(i), reverse = True) #sort by i

	for line in sortedData:
		line[i] = str(line[i])

	#print len(sortedData)
	#output
	sortedFile = open(pFileName + '.sorted', 'w')
	for line in sortedData:
		sortedFile.write('\t'.join(line) + '\n')
	sortedFile.close()

		
if __name__ == "__main__":
	import sys
	
	if len(sys.argv) > 2:
		finalSort(sys.argv[1], sys.argv[2])
	else:
		finalSort(sys.argv[1])	
