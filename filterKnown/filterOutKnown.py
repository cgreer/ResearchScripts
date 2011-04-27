#Get a list of mature predictions that overlap already known microRNAs...

import bioLibCG as cg
import cgConfig
import compareData as compare

def filterOut(cName = None):
	
	#Init
	conf = cgConfig.getConfig(cName)
	
	predictionList = compare.tccFileToList(conf.conf['resultsRaw'], 1)
	#predictionList = compare.tccFileToList(conf.conf['resultsRaw'], 0)
	overlapped = compare.filterOutTccs(predictionList, conf.conf['knownDirectory'], True) #True gives me the filtered out ones instead of the list without filtered out
	
	
	matureOverlaps = open(conf.conf['matureOverlaps'], 'w')
	for tcc in overlapped:
		matureOverlaps.write(tcc + '\n')

if __name__ == "__main__":
	import sys
	import sys
	if len(sys.argv) > 1:
		filterOut(sys.argv[1])
	else:
		filterOut()



