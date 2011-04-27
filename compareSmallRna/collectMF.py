####This script collects the mature microRNA sequence from the folding frames in the mirPipe directory.  It is generally called from
##matchAgainstSmallLibraries and the files are split into packets...
#some frames have the kmers at the end so the sequence grabbed isn't 18 bp long...


from bioLibCG import *
from optparse import OptionParser
import cgConfig

def collectMatureFrames(baseName = None, merLength = None):
	
	#Defaults, Definitions, and Castings
	if not baseName:
		print 'base name of file family needed!'
		return 1
	if not merLength:
		merLength = 18
		
	merLength = int(merLength)


	###############collect IDS that have the same kmer
	print 'Collecting kmer IDs'
	
	conf = cgConfig.cgConfig()
	idFileName = conf.conf['resultsRaw']
	idFile = open(idFileName, 'r')
	interFileName = ('./out/%s/' % baseName) + baseName + '.collection.intermediate'
	interFile = open(interFileName, 'w')
	
	doneList = []

	for line in idFile:
		#For each kmer: grab id and Xmer frames --> output
		kmerID = line.strip().split('\t')[0].split('.')[0]
		if kmerID not in doneList:
			doneList.append(kmerID)
			matureSeq = line.strip().split('\t')[3]
			for frame in returnFrames(matureSeq, merLength):
				interFile.write('%s\t%s\n' % (kmerID, frame))
				
	interFile.close()
	idFile.close()


if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-d", "--inDirectory", dest="inDirectory", help="Directory where all the folding.frames files are")
	parser.add_option("-b", "--baseName", dest="baseName", help="base name of file (LIThuman-s3k8b17) or other")
	parser.add_option("-n", "--packetNumber", dest="packetNumber", help="packet number of file (starts at 1)")
	parser.add_option("-g", "--genome", dest="genome", help="the genome assembly (e.g. hg18 for human)")
	parser.add_option("-m", "--merLength", dest="merLength", help="length of frames to be stripped (defaults to 18 bp)")

	(options, args) = parser.parse_args()
	
	collectMatureFrames(baseName = options.baseName, packetNumber = options.packetNumber, inDirectory = options.inDirectory, genome = options.genome, merLength = options.merLength)

