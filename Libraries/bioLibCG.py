#bioLibCG
import os, time, sys
import subprocess
import urllib
import compareData as compare
import cgConfig as c

letters = ['A','T','C','G','U','N']
humanChromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
mouseChromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY', 'chrM']
zebrafishChromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chrM']

acceptableChroms = [('chr' + str(x)) for x in xrange(1,99)]
acceptableChroms.extend(['chrM', 'chrX', 'chrY', 'chrZ'])


def convertDcdToTcc(dcdList):
	tccList = []
	for dcd in dcdList:
		chr = dcd.strip().split(':')[0]
		strand = dcd.strip().split(':')[1]
		if strand == '+':#switch to numbers if they aren't already.
			strand = '1'
		elif strand == '-':
			strand = '-1'
		start = dcd.strip().split(':')[2].split('-')[0]
		end = dcd.strip().split(':')[2].split('-')[1]
		
		tccList.append('%s:%s:%s:%s' % (chr, strand, start, end))
		
	return tccList
			
###JUST USE REPLACE FXN: (replace('U','T'))

def returnFrames(sequence, frameLength):
	'''Given a sequence and length of frame, return a list of frames from sequence'''
	
	#check for idiots
	if frameLength > len(sequence):
		print 'the length of the frame is larger than the sequence itself!!!'
		return 1
	
	frames = []
	i = 0
	
	while (i+frameLength) < (len(sequence)+1):
		frames.append(sequence[i:(i + frameLength)])
		i = i + 1
		
	return frames
	
def returnFramesDict(sequence, frameLength):
	'''Given a sequence and length of frame, return a dict with the frame and number of times it appeared'''
	
	#check for idiots
	if frameLength > len(sequence):
		print 'the length of the frame is larger than the sequence itself!!!'
		return 1
	
	frameDict = {}
	i = 0
	while (i+frameLength) < (len(sequence)+1):
		frame = sequence[i:(i + frameLength)]
		if frame in frameDict:
			frameDict[frame] = frameDict[frame] + 1
		else:
			frameDict[frame] = 1
		i = i + 1
		
	return frameDict

def getFastQSeqs(inFile):
	fastFile = open(inFile, 'r')
	
	seqs = []
	
	for line in fastFile:
		chars = list(line.strip())
		passed = True
		for char in chars: #might want to test only 5 or 6 letters...nah
			if char not in letters:
				passed = False
				break
		
		if passed:
			seqs.append(line.strip())
	
	fastFile.close()
	return seqs
			
def getFnaSeqs(inFile):
	fnaFile = open(inFile, 'r')
	
	seqs = []
	holdingSeq = "none"
	
	for line in fnaFile:
		if line.startswith('>') and (holdingSeq != "none"):
			seqs.append(holdingSeq)
			holdingSeq = ""
		else:
			passed = True
			chars = list(line.strip())
			for char in chars:
				if char not in letters:
					passed = False
					break
			
			if passed: holdingSeq = holdingSeq + line.strip()
	
	return seqs
		
			
def getDirFiles(directory, start = None, end = None):
	'''return a list of Absolute file names. NOTE: absolute path should NOT contain following slash "/"'''
	#get rid of last slash if there
	if directory.endswith('/'):
		directory = directory[0:-1]
		
	files = os.listdir(directory)
	finalFiles = []
	if start or end:
		for i,file in enumerate(files):
			check = True
			if start:
				if not file.startswith(start):
					check = False
			if end:
				if not file.endswith(end):
					check = False
			
			if check:
				finalFiles.append('%s/%s' % (directory, file))
	else:
		for file in files:
			finalFiles.append('%s/%s' % (directory, file))
	
	return finalFiles
	
def recurseDir(directory, start = None, within = None, end = None, finalFiles = [], top = True, addDirs = True):
	'''return a list of Absolute file names in all directories under one passed
	 NOTE: absolute path should NOT contain following slash "/"
	 the top option should be passed when using --> it's a bad fix for a memory problem...'''
        if top:
	        finalFiles = []
	
        #get rid of last slash if there
        if directory.endswith('/'):
                directory = directory[0:-1]
		
        files = os.listdir(directory)
        dirs = []
        dirFiles = []
        for d in files:#split them up so you go through directories first --> else duplicates :(
                if os.path.isdir(os.path.join(directory,d)):
                        dirs.append(d)
                else:
                        dirFiles.append(d)
			
        for d in dirs:#if it's a directory, go through that directory
                

                for f in recurseDir(directory + '/' + d, start, end, finalFiles, top = False):
                        if f not in finalFiles:
                                finalFiles.append(f)
                                #print 'd append', f
				
        for f in dirFiles:
                if start or end:
                        check = True
                        if start:
                                if not f.startswith(start):
                                        check = False
                        if within:
                                if not within in f:
                                        check = False

                        if end:
                                if not f.endswith(end):
                                        check = False
                        if check:
                                finalFiles.append('%s/%s' % (directory, f))
                else:
                        finalFiles.append('%s/%s' % (directory, f))
                        #print 'f append', directory, f
			
        return finalFiles


def cgprint(string, fam = 0, showList = None):
	'''level printing solution'''
	if showList is None:
		print string
	else:
		if fam in showList:
			print string

def queryJobsDone(parseFileName = None, baseName = None):
	
	qFileName = '%s.qinfo.txt' % baseName
	
	parseFile = open(parseFileName, 'r')
	jobIDs = []
	for line in parseFile:
		jobIDs.append(line.split(' ')[2])

	#Get qstat info every 3 seconds and check if jobid is in the queue.
	#Also, you should have just learned about python threads...whoops
	while True:
		flag = True 
		time.sleep(3)
		#retrieve qfile info
		qFile = open(qFileName, 'w')
		subprocess.Popen(['qstat', '-u','chrisgre'], stdout = qFile).wait()
		qFile.close()
		
		#check qfile info
		qFile = open(qFileName, 'r')
		for line in qFile:
			for id in jobIDs:
				if id in line.split(' '):
					flag = False
		qFile.close()
		
		#if flag is true, no ids present, delete qfile and break
		if flag:
			os.remove(qFileName)
			break

	return True

def stripTripleColon(coord = None):
	
	if coord is None:
		print 'failed to pass coordinate'
		return False
		
	coordDict = {}
	coordDict['chromosome'] = coord.strip().split(':')[0]
	coordDict['strand'] = coord.strip().split(':')[1]
	coordDict['start'] = coord.strip().split(':')[2]
	coordDict['end'] = coord.strip().split(':')[3]
	
	return coordDict
	
def getURL(url, path, gunzip=True):
	"Download url to path, extract if needed"
	try:
		if os.path.exists(path):
			os.remove(path)
		urllib.urlretrieve(url, path)
	except:
		print("\tERROR downloading %s: %s\n" % (url, sys.exc_info()[1]))
		return False
	else: #note: this else statement executes if there are no exceptions
		os.chmod(path, 0660)
		if gunzip and path.endswith('.gz'):
			subprocess.Popen(['gunzip','-q',path]).wait()

def parseUCSC(filePath):
	'''given a ucsc (knownGene) filetype, return a gene structure.
	-I'm not sure if this will give transcipts or genes (I think multiple transcipts)
	-Maybe I could have option to give back tripleColonCoords!
	!!!FOR NOW THIS JUST RETURNS A LIST...'''
	file = open(filePath, 'r')
	
	tccList = []
	for line in file:
		chrom = line.strip().split('\t')[1]
		strand = line.strip().split('\t')[2]
		if strand == '+':
			numStrand = '1'
		else:
			numStrand = '-1'
		geneStart = line.strip().split('\t')[3]
		geneEnd = line.strip().split('\t')[4]
		
		tccList.append('%s:%s:%s:%s' % (chrom, numStrand, geneStart, geneEnd))
		
	return tccList
		
def parseEnsemble(filePath):
	'''given ensemble filetype, return a gene structure.
	-I'm not sure if this will give transcipts or genes (I think multiple transcipts)
	-Maybe I could have option to give back tripleColonCoords!
	!!!FOR NOW THIS JUST RETURNS A LIST...'''
	file = open(filePath, 'r')
	
	tccList = []
	for line in file:
		chrom = line.strip().split('\t')[2]
		strand = line.strip().split('\t')[3]
		if strand == '+':
			numStrand = '1'
		else:
			numStrand = '-1'
		geneStart = line.strip().split('\t')[4]
		geneEnd = line.strip().split('\t')[5]
		
		tccList.append('%s:%s:%s:%s' % (chrom, numStrand, geneStart, geneEnd))
		
	return tccList
		
def parseRefSeq(filePath):
	'''given refseq filetype, return a gene structure.
	-I'm not sure if this will give transcipts or genes (I think multiple transcipts)
	-Maybe I could have option to give back tripleColonCoords!
	!!!FOR NOW THIS JUST RETURNS A LIST...'''
	file = open(filePath, 'r')
	
	tccList = []
	for line in file:
		chrom = line.strip().split('\t')[2]
		strand = line.strip().split('\t')[3]
		if strand == '+':
			numStrand = '1'
		else:
			numStrand = '-1'
		geneStart = line.strip().split('\t')[4]
		geneEnd = line.strip().split('\t')[5]
		
		tccList.append('%s:%s:%s:%s' % (chrom, numStrand, geneStart, geneEnd))
		
	return tccList

def checkOverlap(tcc):
	'''This is for when using a mask file'''
	split = tcc.split(':')
	chrom = split[0]
	strand = split[1]
	start = int(split[2])
	end = int(split[3])

	if strand == '1':
		check = '1'
	else:
		check = '2'

	checkFile = open('/u/home8/gxxiao/chrisgre/Masks/Masks/%s.mask' % chrom , 'r')
	checkFile.seek(start, 0)
	checkString = checkFile.read(end - start)
	checkFile.close()	
	
	if check in checkString:
		return True
	else:
		return False

def getChromosomeLength(chromosome, assembly):
	'''Count the number of lines in a fasta file to figure out chromosome length'''
	chromFileName = '/home/chrisgre/apps/genomes/' + assembly + '/CHR/' + chromosome + '.fa'
	chromFile = open(chromFileName, 'r')
	lines = chromFile.readlines()
	chromLength = (len(lines)*50) + len(lines[-1]) - 101 #101 = 50 + 50 + 1 is for last line, first line, break character
	chromFile.close()
	return chromLength

def createHitMap(coordList, multi = True):
	hitMap = {}
	for coord in coordList:
        	coordDict = stripTripleColon(coord)

		#This is for a multidimensional hitmap:
		##That means that for each hit coordinate their could be more than one originating sequence
		###format: hitMap[i] = [tcc,tcc,etc.].  This is different than non-multi: hitMap[i] = tcc
		if multi:
	        	for i in range(int(coordDict['start']), int(coordDict['end'])):
        	        	if i in hitMap: #already has list, add coords to list
                	        	hitMap[i].append(coord)
		                else: # doesn't have list, give it one and add current coord
        		                hitMap[i] = [coord]
		else: #a non-multi hitmap will just overwrite the originating sequence
	        	for i in range(int(coordDict['start']), int(coordDict['end'])):
				hitMap[i] = coord
	return hitMap

'''def simpleOverlap(s1, e1, s2, e2):
	if (int(s2) >= int(s1) and (int(s2) <= int(e1))) or ((int(e2) >= int(s1)) and (int(e2) <= int(e1))):
		return True
	else:
		return False


def simpleOverlap(s1, e1, s2, e2):
	s1 = int(s1)
	e1 = int(e1)
	s2 = int(s2)
	e2 = int(e2)
	if (s1 <= e2) and (e1 >= s2):
		return True
	else:
		return False
'''

def simpleOverlap(s1, e1, s2, e2):
	s1 = int(s1)
	e1 = int(e1)
	s2 = int(s2)
	e2 = int(e2)
	if (s1 <= e2) and (e1 >= s2): # it overlaps, now get the amount it overlaps
		#print '@ Got here'
                overlap = e2 - s1 + 1
		o2 = e1 - s2 + 1
		if o2 < overlap:
			overlap = o2
			
		#if inside then just set it to length
		l1 = e1 - s1 + 1
		l2 = e2 - s2 + 1
		
		if overlap > l1:
			overlap = l1
		if overlap > l2:
			overlap = l2
		
	        #print '@', overlap	
		return overlap #if greater than 0 --> True
	else:
		return 0


def tccOverlap(tcc1, tcc2, amount = False):
	oneCoords = stripTripleColon(tcc1)
	twoCoords = stripTripleColon(tcc2)
	if (oneCoords['chromosome'] == twoCoords['chromosome']) and (oneCoords['strand'] == twoCoords['strand']):
		overlap = simpleOverlap(oneCoords['start'], oneCoords['end'], twoCoords['start'], twoCoords['end'])
		if amount:
			return overlap #0 is False and anything else is True
		else:
			if overlap:
				return True
			else:
				return False
	else:
		if amount:
			return 0
		else:
			return False
		
class cgTimer:
	iTime = 0
	lReport = 0
	def start(self):
		self.iTime = time.time()
		self.lReport = time.time()
		
	def report(self):
		'''Returns time that has expired'''
		nTime = time.time() - self.iTime
		return nTime
	
	def split(self):
		'''returns time elapsed since last split called'''
		nReport = time.time() - self.lReport
		self.lReport = time.time()
		return nReport
		
def sortTccList(tccList):
	'''sort tcc list by lowest coordinate'''
	sortDict = {}
	lowCoordList = []
	for tcc in tccList:
		lowCoord = int(tcc.split(':')[2])
		lowCoordList.append(lowCoord)
		sortDict[lowCoord] = tcc
		
		
	lowCoordList.sort()
	
	sortedList = []
	for coord in lowCoordList:
		sortedList.append(sortDict[coord])
		
	return sortedList

def getTccLength(tcc):
	l = int(tcc.strip().split(':')[3]) - int(tcc.strip().split(':')[2])
	return l
	
def getTccListTotalLength(tccList):
	length = 0
	
	for tcc in tccList:
		length = length + getTccLength(tcc)
	
	return length
	
def convertTccListToBed(tccList):
	'''convert tcc format to bed format'''
	bedList = []
	for tcc in tccList:
		coordDict = stripTripleColon(tcc)
		if coordDict['strand'] == '1':
			strand = '+'
		else:
			strand = '-'
			
		bedList.append('%s\t%s\t%s\t.\t.\t%s\n' % (coordDict['chromosome'], coordDict['start'], coordDict['end'], strand))
	
	return bedList

def convertTccListToGff(tccList):
	'''convert tcc format to gff format'''
	gffList = []
	for tcc in tccList:
		coordDict = stripTripleColon(tcc)
		if coordDict['strand'] == '1':
			strand = '+'
		else:
			strand = '-'
			
		gffList.append('%s\t.\t.\t%s\t%s\t.\t%s\t.\t.\n' % (coordDict['chromosome'], coordDict['start'], coordDict['end'], strand))
		
	return gffList

def convertTccFileToGff(tccFileName):
	tccList = compare.tccFileToList(tccFileName, 0)
	gffList = convertTccListToBed(tccList)
	
	gffFileName = tccFileName + '.gff'
	gffFile = open(gffFileName, 'w')
	for line in gffList:
		gffFile.write(line)
	gffFile.close()
	
def convertTccFileToBed(tccFileName):
	tccList = compare.tccFileToList(tccFileName, 0)
	bedList = convertTccListToBed(tccList)
	
	bedFileName = tccFileName + '.bed'
	bedFile = open(bedFileName, 'w')
	for line in bedList:
		bedFile.write(line)
	bedFile.close()

def convertGffFileToTcc(gffFileName):
	gffFile = open(gffFileName, 'r')
	tccList = []
	for line in gffFile:
		if line.startswith('#'): continue
		chrom = line.strip().split('\t')[0]
		strand = line.strip().split('\t')[6]
		start = line.strip().split('\t')[3]
		end = line.strip().split('\t')[4]
		
		if not chrom.startswith('chr'):
			chrom = 'chr' + chrom
		if not chrom in acceptableChroms: continue
		if strand == '+' or strand == '1':
			strand = '1'
		else:
			strand = '-1'
		
		tccList.append('%s:%s:%s:%s\n' % (chrom, strand, start, end))
	gffFile.close()
	
	tccFileName = gffFileName + '.tcc'
	tccFile = open(tccFileName, 'w')
	for line in tccList:
		tccFile.write(line)
		
def convertBedFileToTcc(bedFileName):
	bedFile = open(bedFileName, 'r')
	tccList = []
	for line in bedFile:
		if line.startswith('#'): continue
		chrom = line.strip().split('\t')[0]
		strand = line.strip().split('\t')[5]
		start = line.strip().split('\t')[1]
		end = line.strip().split('\t')[2]
		
		if not chrom.startswith('chr'):
			chrom = 'chr' + chrom
		if not chrom in acceptableChroms: continue
		if strand == '+' or strand == '1':
			strand = 1
		else:
			strand = '-1'
		
		tccList.append('%s:%s:%s:%s\n' % (chrom, strand, start, end))
	bedFile.close()

        tccFileName = bedFileName + '.tcc'
	tccFile = open(tccFileName, 'w')
	for line in tccList:
		tccFile.write(line)

def getBaseFileName(absFileName, naked = False):
	#Given the absolute file name, return the base of the file
	'''Naked returns the base file name without period additions'''
	
	if naked:
		return absFileName.strip().split('/')[-1].split('.')[0]
	else:
		return absFileName.strip().split('/')[-1] 
	
def getMetaFileDict(fileName):
	'''Meta files are tsv files with a filename as the first column and attributes after that'''
	mFile = open(fileName, 'r')
	
	#strip the basenames of the files, including the periods...
	metaDict = {}
	for line in mFile:
		baseName = line.strip().split('\t')[0].split('.')[0]
		data = line.strip().split('\t')[1:]
		metaDict[baseName] = data
	
	return metaDict

def appendToLine(line, data, position, sep = '\t'):
	'''Position is in 0 based format (start counting from 0)'''
	
	if line.endswith('\n'):
		returnFlag = True
	else:
		returnFlag = False
	
	lineContents = line.strip().split(sep)
	numSlots = len(lineContents) -1 #subtract one for zero order
	
	if position <= numSlots: #just reset the slot already made
		lineContents[position] = str(data)
	else:
		for i in range(numSlots +1, position + 1):
			lineContents.append('.')
		lineContents[position] = str(data)
	
	newString = sep.join(lineContents)
	if returnFlag:
		newString += '\n'
	
	return newString

def ss(string, spacer = '\t'):
	if spacer == 1:
		spacer = ':'
	'''this is a shortcut for stripsplit'''
	return string.strip().split(spacer)

def makeTcc(chrom, strand, start, end):
	
	return '%s:%s:%s:%s' % (chrom, strand, start, end)


def tccSplit(tcc, text = False):
	'''given tcc -> tuple of parameters'''
	
	chrom = tcc.strip().split(':')[0]
	strand = tcc.strip().split(':')[1]
	start = int(tcc.strip().split(':')[2])
	end = int(tcc.strip().split(':')[3])
	
	return (chrom, strand, start, end)

def switchStrandFormat(strand):

        if strand == '+':
                strand = '1'
        elif strand == '-':
                strand = '-1'
        elif strand == '1':
                strand = '+'
        elif strand == '-1':
                strand = '-'

        return strand

def convertToAS(tcc):
	chrom, strand, start, end = tccSplit(tcc)
	if strand == '1':
		strand = '-1'
	else:
		strand = '1'
	
	return makeTcc(chrom, strand, start, end)
	
def returnChromLengthDict(assembly):
	mConf = c.getConfig('Main.conf')
	lenFileName = mConf.conf['chromosomeLengths'] + '/' + assembly
	f = open(lenFileName, 'r')
	lenDict = {}
	for line in f:
		lenDict[line.split('\t')[0]] = int(line.split('\t')[1])
	
	return lenDict

def clearDirectory(dir, overwrite = True):
	'''Clears the directory if it exists and makes it if it doesnt
	Directory DOES NOT HAVE TRAILING SLASH'''
	
	if os.path.exists(dir):
		if not overwrite: return
		for f in os.listdir(dir):
			os.remove('%s/%s' % (dir, f))
	else:
		os.mkdir(dir)	

def submitArgs(fxn, args):
        '''up to 9, worst written function in history of computers... but pretty :P '''
		
	largs = len(args)
	if largs == 2:
		return fxn(args[1])
	elif largs == 3:
		return fxn(args[1], args[2]) 
	elif largs == 4:
		return fxn(args[1], args[2], args[3])
	elif largs == 5:
		return fxn(args[1], args[2], args[3], args[4]) 
	elif largs == 6:
		return fxn(args[1], args[2], args[3], args[4], args[5]) 
	elif largs == 7:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6]) 
	elif largs == 8:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7]) 
	elif largs == 9:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8]) 
	elif largs == 10:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9]) 
	elif largs == 11:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10]) 
	elif largs == 12:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11]) 
	elif largs == 13:
		return fxn(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12]) 

def uniqueList(list):
        uniq = []
        for item in list:
                if item not in uniq:
                        uniq.append(item)

        return uniq                        

class dominantSpotter:

        def __init__(self, dominantList):
                '''dominant list has the items in dominant order with most dominant at pos 0'''
                self.item_spot = {}
                for i, item in enumerate(dominantList):
                        self.item_spot[item] = i
                self.spot_item = dict((y,x) for x,y in self.item_spot.items()) 


        def spotItem(self, itemList):
                
                #create numbered list
                numberedList = [self.item_spot[x] for x in itemList]
                numberedList.sort()
                return self.spot_item[numberedList[0]]

class MyError(Exception):
        def __init__(self, value):
                self.value = value
        def __str__(self):
                return repr(self.value)


def recursePaths(dir, start = None, end = None):

        paths = []
        for walk in os.walk(dir):
                if walk[0].endswith(end): paths.append(walk[0])

        return paths                

def getNumFileLines(fn):
	f = open(fn, 'r')
        for i, l in enumerate(f):
                pass
        f.close()

	return i + 1

def prettyTime(s):
        hours = s/3600
        lH = s % 3600
        minutes = lH / 60
        lM = lH % 60
        return '%sh %sm %.1fs' % (hours, minutes, lM)
