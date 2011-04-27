import subprocess
import time, os


file = open('config.conf', 'r')

conf = {}
for line in file:
	if line[0] != '>':
		continue
	key = line.strip().split('::')[0][1:]
	value = line.strip().split('::')[1]

	conf[key] = value

#Show configuration:
print '\nConfiguration:'
for entry in conf.keys():
	print '..%s: %s' % (entry, conf[entry])

####################################################
#Files and File Naming
####################################################

#scratchfile
scratchFile = open('./scratch.txt', 'w') #This is just the file that the warnings for the perl script are redirected to

#Filenaming
pipeName = conf['runName'] + '-' + 's' + conf['frameStep'] + 'k' + conf['kmerLength'] + 'b' + conf['mirBasePairs']
conservedOut = 'out/' + pipeName + '.ALL.conserved.kmers.tsv'


files = os.listdir('out/')
for i,file in enumerate(files):
	if file.endswith('split.kmers.tsv') & file.startswith(pipeName):
		packet = file.split('.')[1]
		print '  extracting packet # %s' % packet
		subprocess.Popen(['python', './extract.folding.frames.py',
					'-i', 'out/' + file,
					'-g', conf['genomes'],
					'-w', conf['windowLength'],
					'-f', conf['frameLength'],
					'-s', conf['frameStep'],
					'-o', 'out/' + pipeName + '.' + packet + '.folding.frames.tsv']).wait()

print "DONE"
