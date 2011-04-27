import bioLibCG as cg
import cgConfig

def addPeriods(cName = None):
	
	#init
	conf = cgConfig.getConfig(cName) #gets the current configuration instructions
	
	pFileName = conf.conf['resultsRaw']
	nFileName = conf.conf['results']
	
	pFile = open(pFileName, 'r')
	nFile = open(nFileName, 'w')
	newLines = []
	for line in pFile:
		newLine = '\t'.join(line.strip().split('\t')[0:5]) + '\t.\t.\t.\t.\t.\t.\n'
		nFile.write(newLine)

	
if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		addPeriods(sys.argv[1])
	else:
		addPeriods()
