#Stage 2
import subprocess
import time, os
from optparse import OptionParser
import cgConfig as c

def stageTWO(packetNumber = None):
	
	#Caste and defaults
	if not packetNumber:
		print 'need number of packet to run'
		return 1
	else:
		packetNumber = int(packetNumber)
		
		
		
	#####################################################
	#Load CONFIGURATION FILE:
	#####################################################

	conf = c.cgConfig()
	'''
	print '\nConfiguration:'
	for entry in conf.conf:
		print '..%s: %s' % (entry, conf.conf.conf[entry])
		'''

	
	####################################################
	#Files and File Naming
	####################################################

	#scratchfile
	scratchFile = open('./scratch.txt', 'w') #This is just the file that the warnings for the perl script are redirected to
	
	#directories
	outDir = conf.conf['outDirectory']

	#Filenaming
	#!!!THIS IS DIFFERENENT THAN MIRPIPEPARA!!! (pipename is the name of packet...)
	pipeName = conf.conf['runName'] + '-' + 's' + conf.conf['frameStep'] + 'k' + conf.conf['kmerLength'] + 'b' + conf.conf['mirBasePairs'] + '.' + str(packetNumber)
	extractOut = outDir + '/' + pipeName + '.folding.frames.tsv'
	splitDirectory = outDir + '/' + pipeName + '/'
	foldOut = outDir + '/' + pipeName + '.folded.frames.txt'

	####################################################
	#Pipeline
	####################################################

	print "\nSTARTING STAGE TWO\n(packet # %s)" % packetNumber

				
	if int(conf.conf['numSplitFiles']) == 1: #if you want to run on one node
		print '-Folding Frames using RNAfold on SINGLE NODE'
		subprocess.Popen(['./RNAfold.sh', extractOut, foldOut]).wait()
	else: #Else parallize
		print '-Splitting folding frames into %s files (%s)' % (conf.conf['numSplitFiles'], packetNumber)
		subprocess.Popen(['python', './splitdb.py',
						'-i', extractOut,
						'-b', pipeName,
						'-n', conf.conf['numSplitFiles'],
						'-d', splitDirectory]).wait()
						
		print '-Submitting %s seperate jobs to cluster (%s)' % (conf.conf['numSplitFiles'], packetNumber)
		subprocess.Popen(['python', './parFold.py',
						'-b', pipeName,
						'-n', conf.conf['numSplitFiles'],
						'-d', splitDirectory]).wait()

		#Get Job ID's:
		parseFile = open('%s.jobsinfo.txt' % pipeName, 'r') #!!! edit this?
		jobIDs = []
		for line in parseFile:
			jobIDs.append(line.split(' ')[2])
			
		#check if right number were submitted
		if len(jobIDs) == int(conf.conf['numSplitFiles']):
			print '....jobs were submitted correctly'

	
	print 'DONE (%s)' % packetNumber
	

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-n", "--number", dest="packetNumber", help="number of the packet to run")
	(options, args) = parser.parse_args()
	
	#splitdb(fIN = 'out/Step3k8.folding.frames.tsv', directory = 'out/test/', numOut = 32, basename = 'test')
	stageTWO(packetNumber = options.packetNumber)

