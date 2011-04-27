#This script will compare the microRNA sequences from intermediate files to Xmer sequences in 
#small Rna libraries
#WITH COUNTS FILES.

from bioLibCG import *
from optparse import OptionParser

def compareCounts(baseName = None, merLength = None, smallPath = None):
	if not baseName:
		print 'NEED baseName'
		return 1
	if not merLength:
		print 'Assuming merLength is 18'
		merLength = 18
	if not smallPath:
		print 'No Small Path'
		return 1
		
	merLength = int(merLength)

	print 'Comparing relevent microRNA frames against small library data'
	#create mirData dictionary
	mirData = {}
	for parseFileName in getDirFiles('./out/%s' % baseName, start = baseName, end = 'collection.intermediate'): #different parseFile than above
		parseFile = open(parseFileName, 'r')
		for line in parseFile:
			hitID = line.strip().split('\t')[0]
			hitFrame = line.strip().split('\t')[1]
			mirData[hitFrame] = hitID
		parseFile.close()
	
	#Have to add addition of nonUniqueKmers
	
	#now compare every counts file against mature frames.
	outFileName = ('./out/%s/' % baseName) + baseName + '.compared'
	outF = open(outFileName, 'w')

	for countFileName in recurseDir(smallPath, end = '.counts'):
		print countFileName
		countFile = open(countFileName, 'r')
		for line in countFile:
			seq = line.strip().split('\t')[0]
			count = int(line.strip().split('\t')[1])
			if len(seq) >= merLength:
				for seqFrame in returnFrames(seq, merLength):
					if seqFrame in mirData:
						i = 0
						while i < count:
							outF.write('%s:%s:%s\n' % (countFileName, mirData[seqFrame], seqFrame))
							i = i +1

						
if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-m", "--mlength", dest="merLength", help="length of frame to compare")
	parser.add_option("-b", "--baseName", dest="baseName", help="base name of file")
	(options, args) = parser.parse_args()
	
	compareFinished(baseName = options.baseName, merLength = options.merLength)


