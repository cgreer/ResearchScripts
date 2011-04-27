####This script collects the mature microRNA sequence from the folding frames in the mirPipe directory.  It is generally called from
##matchAgainstSmallLibraries and the files are split into packets...
#some frames have the kmers at the end so the sequence grabbed isn't 18 bp long...


from bioLibCG import *
from optparse import OptionParser


def collectMatureFrames(baseName = None, packetNumber = None, inDirectory = None, genome = None, merLength = None):
	
	#Defaults, Definitions, and Castings
	if not inDirectory:
		print 'Input directory needed!'
		return 1
	if not baseName:
		print 'base name of file family needed!'
		return 1
	if not packetNumber:
		print 'the packet number is needed!'
		return 1
	if not merLength:
		merLength = 18
	if not genome:
		print 'please provide genome assembly!'
		return 1
		
	packet = packetNumber
	merLength = int(merLength)

	###############collect IDS that have the same kmer
	print 'Collecting relevent IDs'

	idFileName = inDirectory + '/' + baseName + '.ALL.FINAL.mirs.tsv'
	idFile = open(idFileName, 'r')

	finalIDs = []
	kmerIDs = []

	for line in idFile:
		#strip id and add to list if kmer isn't already present...
		fullID = line.strip().split('\t')[0] + '.' + line.strip().split('\t')[4]
		kmerID = fullID.split('.')[0] + '.' + fullID.split('.')[2]
		frameID = fullID.split('.')[1]
		if (kmerID not in kmerIDs) and (int(frameID) > 3): #frames with .3 or lower won't be long enough for 18mer
			finalIDs.append(fullID)
			kmerIDs.append(kmerID)

	print '  Number of IDs: %s' % len(finalIDs)

	idFile.close()


	################Collect mature MiRNA sequences and their frames from finalIDs
	print 'Collecting frames of all relevent mature microRNAs'

	
	interFileName = ('./out/%s/' % baseName) + baseName + '.' + packet + '.collection.intermediate'
	interFile = open(interFileName, 'w')

	count = 0

	fFile = open('%s%s.%s.folding.frames.tsv' % (inDirectory, baseName, packet), 'r')
	print '  Using folding frame file %s' % packet
	lineID = 'none'
	gen = 'none'
	grabStart = 'none'
	mirSeq = 'none'
	kmer = 'none'
		
	for line in fFile:
		if line.startswith('>'): # parse it
			lineID = line.strip().split(' ')[1].split('.')[0] + '.' + line.strip().split(' ')[1].split('.')[1] + '.' + line.strip().split(' ')[1].split('.')[6]
			gen = line.strip().split(' ')[1].split('.')[2]
			grabStart = abs(int(line.strip().split(' ')[1].split('.')[7]) - int(line.strip().split(' ')[1].split('.')[6]))
			myMirChrom = line.strip().split(' ')[1].split('.')[3]
			myMirStrand = line.strip().split(' ')[1].split('.')[4]
			myMirStart = line.strip().split(' ')[1].split('.')[7]
		elif line == '\n':
			continue
		else: #Collect if valid
			if (lineID in finalIDs) and (gen == genome):
				mirSeq = UtoT(line.strip()[grabStart:grabStart + 22])
				mirFrames = returnFrames(mirSeq, 18)
				if mirFrames == 1: #sequence not long enough -> filter out
					print mirSeq, lineID, gen
				else:
					for frame in mirFrames:
						interFile.write('%s.%s.%s.%s\t%s\n' % (lineID, myMirChrom, myMirStrand, myMirStart, frame))
				count = count + 1
			
					
		
	fFile.close()
	interFile.close()
	print count

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

