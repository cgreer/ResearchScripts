#submit parallel jobs for split dbs.

from optparse import OptionParser
import subprocess


def submitJobs(directory = None, numFiles = None, basename = None):
	
	if directory is None:
		print 'No directory for job submission!'
		return 1
	if numFiles is None:
		print 'Need number of files for job submission'
		return 1
	if basename is None: basename = 'NONAME'
	
	#caste
	numFiles = int(numFiles)
	
	parseFile = open('%s.jobsinfo.txt' % basename, 'w')
	#submit jobs
	for i in range(1, numFiles + 1):
		inFile = directory + basename + '.' + str(i) + '.tsv'
		outFile = directory + basename + '.' + str(i) + '.out'
		
		
		subprocess.Popen(args=['qsub','-cwd', '-e', 'folderrors.txt', '-o', 'folderrors.txt', '-N', '%s-%s' % (basename, i),
					'-V', 'RNAfold.sh',
					inFile, outFile,], stdout = parseFile ).wait()
		
	
	parseFile.close()


if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-n", "--numfiles", dest="numFiles", help="number of files to split into")
	parser.add_option("-d", "--directory", dest="directory", help="directory split files will go into")
	parser.add_option("-b", "--basename", dest="basename", help="base of file name")
	(options, args) = parser.parse_args()
	
	#submitJobs(directory = 'out/test/', numFiles = 32)
	submitJobs(directory = options.directory, numFiles = options.numFiles, basename = options.basename)
