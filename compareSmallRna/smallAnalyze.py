#counts the amount of unique hits (not clusters) in small RNA libraries and prints the ID and AMOUNT of hits for each unique hit
'''#For every ID: output amount of hits per library for that ID
#ID
 Lib # hits
 Lib # hits
'''

from optparse import OptionParser
import bioLibCG as cg

def smallAnalyze(inFile = None):
	
	#Caste and defaults
	if not inFile:
		print 'need inFile'
		return 1
	
	# start Timer
	timer = cg.cgTimer()
	timer.start()
	
	# for every id in outfile, count how many matches there are
	countDict = {} 
	countFile = open(inFile, 'r')
	
	for line in countFile:
		(id, library) = (line.strip().split(':')[1], line.strip().split(':')[0])
		if id not in countDict:
			countDict[id] = {}
		else:
			if library not in countDict[id]:
				countDict[id][library] = 1
			else:
				countDict[id][library] = countDict[id][library] + 1
	#print 'Time for counting lib hits: ', timer.split()
	
	sortList = []
	for id in countDict:
		sortList.append(id)
	sortList.sort()
	
	for id in sortList:
		print '%s' % id
		for lib in countDict[id]:
			print '%s\t%s' % (lib, countDict[id][lib])

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-i", "--inFile", dest="inFile", help="input file")
	(options, args) = parser.parse_args()
	
	smallAnalyze(inFile = options.inFile)
