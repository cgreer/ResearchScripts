#Given tsv database, split it into X number of tsv databases and store
'''
Notes:
The big if else statement in the for loop is to determine if the current line is the first line
of an entry.  For example, if in your database there is a description line starting with a '>' and 
a sequence starting with AAAAetc, the description must go with the sequence and be in the same file.
In other words, DON't start writing to a new file until you hit a description line...

'''

import math, os
from optparse import OptionParser

def splitdb(fIN = None, directory = None, numOut = None, filebase = None):
	
	if fIN is None:
		print 'No Input File for splitting db!'
		return 1
	if directory is None:
		print 'No directory to save split db files!'
	if numOut is None:
		print 'Must define how many pieces db will split into'
	if filebase is None: filebase = 'test'
	
	inFile = open(fIN, 'r')
	numOut = int(numOut)
	
	#put lines in list and get amount of lines
	#STOP PUTTING FILE IN LIST
	#fileList = inFile.readlines()
	#numLines = len(fileList)
	
	#Get amount of lines per file
	#TO DO THIS: iterate over file and count
	numLines = 0
	for line in inFile:
		numLines = numLines + 1
	linesPerFile = int(math.ceil(float(numLines)/float(numOut)))
	
	inFile.seek(0)#go back to beginning of file
	
	
	
	#make directory if it doesn't exist or delete all files in it
	if os.path.exists(directory):
		for f in os.listdir(directory):
			os.remove('%s%s' % (directory, f))
	else:
		os.mkdir(directory)	
	
	#split the db into X files.
	currFileName = None
	switchFileName = None
	outFile = None
	for i,line in enumerate(inFile):
		
		#Switch file names if lines need it...
		if line[0] == '>': #desc line, can switch files if necessary
			filenum = int(math.ceil(float(i+1)/float(linesPerFile))) #plus one is to avoid zero
			switchFileName = directory + filebase + '.' + str(filenum) + '.tsv'
			#print switchFileName
			#print filenum
			if switchFileName != currFileName:
				if outFile: outFile.close()
				#print switchFileName
				outFile = open(switchFileName, 'w')
				currFileName = switchFileName
		
		#write to file
		outFile.write(line)
		
	


if __name__ == "__main__":
	import sys
	parser = OptionParser()
	parser.add_option("-i", "--infile", dest="inFile", help="input file")
	parser.add_option("-b", "--basename", dest="basename", help="base of file name")
	parser.add_option("-n", "--numfiles", dest="numOut", help="number of files to split into")
	parser.add_option("-d", "--directory", dest="directory", help="directory split files will go into")
	(options, args) = parser.parse_args()
	
	#splitdb(fIN = 'out/Step3k8.folding.frames.tsv', directory = 'out/test/', numOut = 32, basename = 'test')
	splitdb(fIN = options.inFile, directory = options.directory, numOut = options.numOut, filebase = options.basename)

