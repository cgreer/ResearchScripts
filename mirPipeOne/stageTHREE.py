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

	print '\nConfiguration:'
	for entry in conf.conf:
		print '..%s: %s' % (entry, conf.conf[entry])


	####################################################
	#Files and File Naming
	####################################################

	#scratchfile
	scratchFile = open('./scratch.txt', 'w') #This is just the file that the warnings for the perl script are redirected to
	
	#directories
	outDir = conf.conf['outDirectory']

	#Filenaming
	pipeName = conf.conf['runName'] + '-' + 's' + conf.conf['frameStep'] + 'k' + conf.conf['kmerLength'] + 'b' + conf.conf['mirBasePairs'] + '.' + str(packetNumber)
	foldOut = outDir + '/'  + pipeName + '.folded.frames.txt'
	filterOut = outDir + '/' + pipeName + '.filtered.mirs.tsv'
	finalOut = outDir + '/' + pipeName + '.FINAL.mirs.tsv'


	####################################################
	#Pipeline
	####################################################

	print "\nSTARTING STAGE THREE\n(packet # %s)" % packetNumber

				

	print "-Filtering frames and finding prospective mirs"
	subprocess.Popen(['python', './filter.frames.py',
					'-i', foldOut,
					'-g', conf.conf['genomes'],
					'-b', conf.conf['mirBasePairs'],
					'-m', conf.conf['mirLength'],
					'-o', filterOut]).wait()
	'''		
	print "-Running Conservation filter"
	subprocess.Popen(['perl', './get_percent_identity_list_fix.pl',
					'-g', conf.conf['genomes'],
					'-l', filterOut,
					'-o', finalOut], stderr = scratchFile).wait()
	'''
	
	print "DONE"

if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-n", "--number", dest="packetNumber", help="number of the packet to run")
	(options, args) = parser.parse_args()
	
	#splitdb(fIN = 'out/Step3k8.folding.frames.tsv', directory = 'out/test/', numOut = 32, basename = 'test')
	stageTWO(packetNumber = options.packetNumber)
