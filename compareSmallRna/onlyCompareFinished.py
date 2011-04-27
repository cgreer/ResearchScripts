#This script will compare the microRNA sequences from intermediate files to Xmer sequences in 
#small Rna libraries.

from bioLibCG import *
from optparse import OptionParser

def compareFinished(baseName = None, merLength = None):
	if not baseName:
		print 'NEED baseName'
		return 1
	if not merLength:
		print 'Assuming merLength is 18'
		merLength = 18
		

	smallPath = '/u/home8/gxxiao/chrisgre/smallRNALibraries'
	merLength = int(merLength)

	print 'Comparing relevent microRNA frames against small library data'
	mirData = {}
	for parseFileName in getDirFiles('./out/%s' % baseName, start = baseName, end = 'collection.intermediate'): #different parseFile than above
		parseFile = open(parseFileName, 'r')

		for line in parseFile:
			hitID = line.strip().split('\t')[0]
			hitFrame = line.strip().split('\t')[1]
			mirData[hitFrame] = hitID
		parseFile.close()

	outFileName = ('./out/%s/' % baseName) + baseName + '.compared'
	outF = open(outFileName, 'w')
	#fast q files have different function for getting sequences
	for smallFile in getDirFiles(smallPath, end = '.fastq'):
		print '  comparing against %s' % smallFile
		for smallSeq in getFastQSeqs(smallFile):
			if len(smallSeq) >= merLength:
				for smallFrame in returnFrames(smallSeq, merLength):
					if smallFrame in mirData:
						outF.write('%s:%s:%s\n' % (smallFile, mirData[smallFrame], smallFrame))

	fileList = getDirFiles(smallPath, end = '.fna')
	fileList.extend(getDirFiles(smallPath, end = '.fa'))
	#fna/fa files have different function for getting sequences
	for smallFile in fileList:
		print '  comparing against %s' % smallFile
		for smallSeq in getFnaSeqs(smallFile):
			if len(smallSeq) >= merLength:
				for smallFrame in returnFrames(smallSeq, merLength):
					if smallFrame in mirData:
						outF.write('%s:%s:%s\n' % (smallFile, mirData[smallFrame], smallFrame))
						
if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-m", "--mlength", dest="merLength", help="length of frame to compare")
	parser.add_option("-b", "--baseName", dest="baseName", help="base name of file")
	(options, args) = parser.parse_args()
	
	compareFinished(baseName = options.baseName, merLength = options.merLength)


