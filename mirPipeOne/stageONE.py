#mirPipe using parrallel cpu's
import subprocess
import time, os
import cgConfig as c



startTime = time.time()
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
pipeName = conf.conf['runName'] + '-' + 's' + conf.conf['frameStep'] + 'k' + conf.conf['kmerLength'] + 'b' + conf.conf['mirBasePairs']
conservedOut = outDir + '/' + pipeName + '.ALL.conserved.kmers.tsv'
inputTccs = outDir + '/' + conf.conf['kmerDatabase']

####################################################
#Pipeline
####################################################

print "\nSTARTING STAGE ONE"
print 'Run Name: %s' % pipeName

'''
print "-Retrieving conserved kmers"
subprocess.Popen(args=['perl','./get_conserved_kmers_list.pl',
					'-g', conf.conf['genomes'],
					'-k', conf.conf['kmerLength'],
					'-c', inputTccs,
					'-o', conservedOut], stderr = scratchFile).wait()
			'''
				
print "-Splitting conserved kmers into manageable file sizes"
subprocess.Popen(['python', './splitdbconserved.py',
					'-i', conservedOut,
					'-b', pipeName,
					'-l', conf.conf['numKmersPerFile'],
					'-d', outDir + '/']).wait()


print "-Extracting folding frames"

files = os.listdir(outDir + '/')
for i,file in enumerate(files):
	if file.endswith('split.kmers.tsv') & file.startswith(pipeName):
		packet = file.split('.')[1]
		print '  extracting packet # %s' % packet
		subprocess.Popen(['python', './extract.folding.frames.py',
					'-i', outDir + '/' + file,
					'-g', conf.conf['genomes'],
					'-w', conf.conf['windowLength'],
					'-f', conf.conf['frameLength'],
					'-s', conf.conf['frameStep'],
					'-o', outDir + '/' + pipeName + '.' + packet + '.folding.frames.tsv']).wait()

print "DONE"
