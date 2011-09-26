#given a file, skip to correct line until value is in range or whatever...
import os
import cgSort
import bioLibCG as cg
chromToNum = cgSort.returnChromToNumDict()

def getBytesPerLine(fN):
	
	f = open(fN, 'r')
	f.readline() #could be header...
	line = f.readline()
	f.close()
	
	return len(line)

def getFileSize(fN):
	
	return os.path.getsize(fN)

def getNumLines(fN, fSize, lSize):
	'''This is approximate!!! based on average length of lines
	and filesize'''
	
	return (fSize/lSize)
	
def getBytesPerPart(numParts, totalBytes):
	'''given the total number of bytes, give back the number of bytes per part'''
	
	return (totalBytes/numParts)  + 1

def levelsOverEnd():
	return [0]

def mapSequenceCheckFunction(val, line):
	lineSequence = line.strip().split('\t')[4]
	
	if val < lineSequence:
		return -1
	elif val > lineSequence:
		return 1
	elif val == lineSequence:
		return 0

def mapStartRangeCheckFunction(val, line):
	lineStart = int(line.strip().split('\t')[3])
        lineEnd = lineStart + len(line.strip().split('\t')[4])
        chrom, strand, start, end = cg.tccSplit(val)
        start = int(start)
        end = int(end)

        if cg.simpleOverlap(start, end, lineStart, lineEnd):
                return 0
        else:
                return -1
        
def mapStartCheckFunction(val, line):
	lineStart = int(line.strip().split('\t')[3])
        chrom, strand, start, end = cg.tccSplit(val)
        start = int(start)

	if start < lineStart:
		return -1
	elif start > lineStart:
		return 1
	elif start == lineStart:
		return 0
	
def wigCheckFunction(val, line):
	#print line.strip(), val
	chrom = line.split('\t')[0]
	low = int(line.split('\t')[1])
	high = int(line.split('\t')[2])
		
	valChrom, valStrand, valStart, valEnd = cg.tccSplit(val)
	
	#switch the chroms into numbering...
	valChromNum = chromToNum[valChrom]
	chromNum = chromToNum[chrom]
	
	#first check if chrom is correct or higher or lower
	if valChromNum > chromNum:
		return 1
	elif valChromNum < chromNum:
		return -1
	elif valChromNum == chromNum: #if the same --> check 
		if valStart >= high:
			return 1
		elif valStart < low:
			return -1
		elif low <= valStart < high:
			return 0

def primaryIDCheckFunction(val, line):
        
        id = int(line.split('\t')[0])
        
        if val > id:
                return 1
        elif val < id:
                return -1
        elif val == id:
                return 0


def testSkip(low, high, check):
	'''Just trying to figure out how to skip and cover whole range...'''
	
	i = low + (high -low)/2
	
	if check < i:
		testSkip(low, i, check)
	elif check == i:
                pass
		#print 'found it', i
	else:
		testSkip(i + 1, high, check)
	
	
		
class lineIndex:
	'''return a specific line from an file that is in line form
	a good example of line form is a wig file'''
	
	def __init__(self, fN, header = False):
		
		self.lineLength = getBytesPerLine(fN)
		self.currentByte = 0
		self.file = open(fN, 'r')
		self.fSize = getFileSize(fN)
		self.lSize = getBytesPerLine(fN)
		self.numLines = getNumLines(fN, self.fSize, self.lSize)
		self.numParts = 2
		self.firstLineStart = 0
		
		if header:
			lineOne = self.file.readline()
			a = len(lineOne)
			self.firstLineStart = a + 1
			self.file.seek(0)
		
		
	
	def setCheckVal(self, checkVal):
		'''This is the value that you are looking for in each line'''
		self.checkVal = checkVal
	
	def getLineFromByte(self, i):
		'''same as other but doesn't move the file cursor'''
		#print '  byte:', i
		j = 1
		while (i - j) > 0:
			self.file.seek(i - j)
			char = self.file.read(1)
			#print '  ', str(i-j), char
			if char == '\n':
				self.file.seek(i - j + 1)
				return self.file.readline()
			j += 1
		
		#return first line
		self.file.seek(0)
		return self.file.readline()
	

	def extendUp(self, val):
		'''STILL HAS ZERO ERROR!!!'''
		
		while self.currentByte > 0:
			potentialSeekPoint = self.currentByte
			
			#test one level up
			upByte = potentialSeekPoint - 1
			checkLine = self.getLineFromByte(i = upByte)
			
			check = self.checkFunction(val, checkLine)
			if check == 0: #keep going
				self.setByteToLineStart(upByte)
				continue
			else: # set file pointer to potential
				self.currentByte = potentialSeekPoint
				self.file.seek(self.currentByte)
				break
					
	def setByteToLineStart(self, i = None):
		
		if not i:
			i = self.currentByte
		
		j = 1
		while (i - j) > 0:
			self.file.seek(i - j)
			if self.file.read(1) == '\n':				
				self.currentByte = i - j + 1 #2 because when I read the one byte I was one behind + 1 to get past \n
				self.file.seek(self.currentByte)
				return
			j += 1
		
		#update cursor
		self.currentByte = 0 #beginning of file...
		self.file.seek(self.currentByte)
		
		
	def passCheckFunction(self, function):
		'''this is the function that will be used to check each line... it must pass true, false, higher or lower...
		it must use the checkValue and the line and make a comparison...'''
		self.checkFunction = function
	
	def passOverEndFunction(self, function):
		'''If the value you passed is over the length of the index, do this'''
		
		self.overEnd = function
	
	def binarySearch(self, val, low = None, high = None, checkEnd = True, skipEnd = False):
		'''SkipEnd will use the binarySkipEnd function for checking files without the value'''
		
		if not low:
			low = self.firstLineStart
		if not high:
			high = self.fSize
		
		if checkEnd:
			i = self.fSize - 1 #last line
			checkLine = self.getLineFromByte(i)
			check = self.checkFunction(val,checkLine)
			
			#!!!Should incorporate an overEnd Fxn
			if check == 1:
				raise NameError('Value Passed is not in index!\nValue: %s\nLast Line:%s' % (val, checkLine))
		
		if skipEnd:
			return self.binarySkipEnd(low, high, val)
		else:
			return self.binarySkip(low, high, val)
		
		
	def binarySkip(self, low, high, val):
				
                if low == high:
                        raise NameError(' (%s) is not in file, use skipEnd if you want' % val)
		i = low + (high -low)/2
		checkLine = self.getLineFromByte(i)
		check = self.checkFunction(val,checkLine)
		
		if check == -1:
			return self.binarySkip(low, i, val)
		elif check == 1:
			return self.binarySkip(i + 1, high, val)
		else:
			self.setByteToLineStart(i)
			return checkLine
	
	def binarySkipEnd(self, low, high, val):
		'''Slower, has an extra check to see if the same low and
		high are being checked.  Used for line files in which the
		value might not exactly be in the file'''
		
		i = low + (high -low)/2
                #print low, high, i
		checkLine = self.getLineFromByte(i)
		check = self.checkFunction(val,checkLine)
                #print i, low, high
		#check if at end
		if i == low or i == high:
			self.setByteToLineStart(i)
			return checkLine
		
		if check == -1:
			return self.binarySkipEnd(low, i, val)
		elif check == 1:
			return self.binarySkipEnd(i + 1, high, val)
		else:
			self.setByteToLineStart(i)
			return checkLine
		
	def scanLine(self, check):
		#read lines until it matches...(should go up and down? or only forward?)
		#for now only forward...
		
		#set byte to current line...
		self.setByteToLineStart()
		
		for line in self.file:
			if self.checkFunction(check, line) == 0:
				return line
		
		return None
	
	def close(self):
		self.file.close()
		
