###This is the MAIN script to compare small RNA libraries to the mature microRNA sequences.  It first calls
##collectMatureFrames and then calls compareFinished to do the actual comparisons.

from bioLibCG import *
import os
from optparse import OptionParser
import subprocess
from onlyCompareFinished import *

def matchAgainstSmallLibraries(baseName = None, inDirectory = None, genome = None, merLength = None):
	
	#Defaults, Definitions, Castings, and Globals
	if not inDirectory:
		print 'Input directory needed!'
		return 1
	if not baseName:
		print 'base name of file family needed!'
		return 1
	if not merLength:
		merLength = 18
	if not genome:
		print 'please provide genome assembly!'
		return 1

	
	smallPath = '/u/home9/gxxiao/chrisgre/smallRNALibraries'
	merLength = int(merLength)
	
	#Get number of packets, submit a job for each packet
	print 'Submitting packets for ID collection'
	foldingFiles = getDirFiles(inDirectory, start = baseName, end = 'folding.frames.tsv')
	numFiles = len(foldingFiles)
	
	parseFile = open('%s.jobsinfo.txt' % baseName, 'w')
	
	#Make directory for files to go into
	filePath = 'out/%s' % baseName
	if baseName in os.listdir('out/'):
		#Delete all contents
		for file in os.listdir(filePath):
			os.remove('%s/%s' % (filePath,file) )
	else:
		os.mkdir(filePath)
	
	#submit jobs
	for i in range(1, numFiles + 1):
		packetNumber = str(i)
		
		if i < 25:
			subprocess.Popen(args=['qsub','-cwd', '-q', 'gxxiao16.q', '-e', 'JobErrors.txt' , '-o', 'JobErrors.txt', '-N', '%s-%s' % (baseName[0:4], i),
						'-V', 'callCollect.sh', inDirectory, baseName, packetNumber, genome, str(merLength)], stdout = parseFile ).wait()
		else:
			subprocess.Popen(args=['qsub','-cwd', '-q', 'gxxiao32.q', '-e', 'JobErrors.txt' , '-o', 'JobErrors.txt', '-N', '%s-%s' % (baseName[0:4], i),
						'-V', 'callCollect.sh', inDirectory, baseName, packetNumber, genome, str(merLength)], stdout = parseFile ).wait()
	
	parseFile.close()
	
	parseFileName = '%s.jobsinfo.txt' % baseName 
	if queryJobsDone(parseFileName, baseName):
		compareFinished(baseName, merLength)
	

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-d", "--inDirectory", dest="inDirectory", help="Directory where all the folding.frames files are")
	parser.add_option("-b", "--baseName", dest="baseName", help="base name of file (LIThuman-s3k8b17) or other")
	parser.add_option("-g", "--genome", dest="genome", help="the genome assembly (e.g. hg18 for human)")
	parser.add_option("-m", "--merLength", dest="merLength", help="length of frames to be stripped (defaults to 18 bp)")

	(options, args) = parser.parse_args()
	
	matchAgainstSmallLibraries(baseName = options.baseName, inDirectory = options.inDirectory, genome = options.genome, merLength = options.merLength)
