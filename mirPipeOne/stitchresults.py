from optparse import OptionParser

def stitchdb(inDirectory = None, numFiles = None, basename = None, outDirectory = None):
	
	if inDirectory is None:
		print 'No inDirectory to save split db files!'
		return 1
	if outDirectory is None:
		print 'No out directory to save split db files!'
		return 1		
	if basename is None: basename = 'test'
	if numFiles is None:
		print 'must specify number of files to stitch'
		return 1
	
	#caste
	numFiles = int(numFiles)
	
	#stitch
	fileList = []
	inFile = None
	currFileName = None
	for i in range(1, numFiles + 1):
		#open file, add lines to linelist
		currFileName = inDirectory + basename + '.' + str(i) + '.filtered.mirs.tsv'
		outFile = open(currFileName, 'r')
		fileList.extend(outFile.readlines())
		outFile.close()
		
	#write list to output file in output directory
	currFileName = outDirectory + basename + '.ALL.FINAL.mirs.tsv'
	outFile = open(currFileName, 'w')
	for line in fileList:
		outFile.write(line)
	outFile.close()
	
if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-b", "--basename", dest="basename", help="base of file name")
	parser.add_option("-i", "--inDirectory", dest="inDirectory", help="inDirectory split files will go into")
	parser.add_option("-o", "--output", dest="outDirectory", help="directory stitched files will go")
	parser.add_option("-n", "--numfiles", dest="numFiles", help="number of files to split into")
	(options, args) = parser.parse_args()
	
	#splitdb(fIN = 'out/Step3k8.folding.frames.tsv', inDirectory = 'out/test/', numOut = 32, basename = 'test')
	stitchdb(inDirectory = options.inDirectory, outDirectory = options.outDirectory, basename = options.basename, numFiles = options.numFiles)
