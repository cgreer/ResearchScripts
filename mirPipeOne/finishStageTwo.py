import cgConfig as c
import subprocess

def finish(packetNumber):
	#init
	conf = c.cgConfig()

	#directories
	outDir = conf.conf['outDirectory']

	#Filenaming
	pipeName = conf.conf['runName'] + '-' + 's' + conf.conf['frameStep'] + 'k' + conf.conf['kmerLength'] + 'b' + conf.conf['mirBasePairs'] + '.' + str(packetNumber)
	extractOut = outDir + '/' + pipeName + '.folding.frames.tsv'
	splitDirectory = outDir + '/' + pipeName + '/'
	foldOut = outDir + '/' + pipeName + '.folded.frames.txt'


	print '-Stitching %s files back into one (%s)' % (conf.conf['numSplitFiles'], packetNumber)
	subprocess.Popen(['python', './stitchdb.py',
					'-b', pipeName,
					'-n', conf.conf['numSplitFiles'],
					'-o', outDir + '/',
					'-i', splitDirectory]).wait()

if __name__ == "__main__":
	import sys
	
	finish(sys.argv[1])
