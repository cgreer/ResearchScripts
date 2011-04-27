
def getFastQType(fName, quick = False):
	'''The quick option is for determing if it is +33 or +64 quickly'''
	
	fastFile = open(fName, 'r')
	sangerFlag = False #Is there a character below 58?
	notSangerFlag = False #Is there a character above 73?
	solexaFlag = False #is notSanger true and a character below 64
	illuminaOneCharFlag = False #This just checks for 64, 65 characters
	isQual = False

	for line in fastFile:
		line = line.strip()
		
		if sangerFlag or solexaFlag:
			break
		if quick and notSangerFlag:
			break
			
		if isQual:#This line contains quality scores
			isQual = False
			
			l = len(line)
			for i in range(0,l):#accrue tags
				aNum = ord(line[i]) #Get the ascii number
				if aNum < 58:
					sangerFlag = True
					break
				if aNum > 73:
					notSangerFlag = True
				if notSangerFlag and (aNum < 64):
					print line
					solexaFlag = True
					break
				if (aNum == 64) or (aNum == 65):
					illuminaOneCharFlag = True
			
						
		if line.startswith('+'):#The next line contains quality scores
			isQual = True
	fastFile.close()
		
	#Deduce what type of file it was
	
	if quick:#if you just want to know if 33 or 64 phred scores
		if sangerFlag:
			return 'Sa'
		else:
			return 'I2'
			
	#try to determine the exact file type.
	if sangerFlag:
		return 'Sa'
	elif solexaFlag:
		return 'So'
	elif notSangerFlag and illuminaOneCharFlag:
		return 'I1'
	else: #Likely I2, but could technically be any other if chars are in specific range.
		return 'I2'	

def isValidFastQ(fName):
	'''DOESN'T currently work for multiline fastq, but most ones today are single lined.
	mainly checks to see if there is a qual score for every nt.'''
	fastFile = open(fName, 'r')
	
	#!FIX THIS
	while True:#Goto first @ Assumes there is a first @...
		if fastFile.readline().startswith('@'):
			break
			
	i = 0
	length = 0
	qualLength = 0
	lengthFlag = True
	
	for line in fastFile:
		if i % 4 == 0: #Sequence Length
			length = len(line.strip())
		if (i-2) % 4 == 0: #Quality line
			qualLength = len(line.strip())
		if (i-3) % 4 == 0: #Check line
			if not (length == qualLength):
				lengthFlag = False
				break
		
		i = i + 1
	
	#Do last check, loop will end without checking last entry
	if not (length == qualLength):
				lengthFlag = False
	
	#Check if passed test
	if lengthFlag:
		return True
	else:
		return False
				
	
	
		
