import os

def createCountFileFastQ(fastQFileName):
	try:
		QFile = open(fastQFileName, 'r')
	except IOError:
		print 'ERROR: OPENING FILE: Counts file will NOT be created'
		return 0
	
	#check to see if there is already a counts file for this file:
	if os.path.isfile(fastQFileName + '.counts'):
		print '%s already exists!  Not overwriting...' % (fastQFileName + '.counts')
		return 0
	
	
	##check if first line begins with @, i = 0 will be "2nd" line because I just read a line
	for line in QFile:
		if line.startswith('@'):
			break
	
	#Count the amount of times a read appears
	seqDict = {}
	i = 0
	for line in QFile:
		if i%4 == 0: #the sequence is on every fourth line,
			seq = line.strip()
			#print seq
			
			if seq in seqDict:
				#print 'inside'
				seqDict[seq] = seqDict[seq] + 1
			else:
				seqDict[seq] = 1
				
		i = i + 1
	QFile.close()
	
	#Now save it into a count file
	#Now save it into a count file
	try:
		newFile = open(fastQFileName + '.counts', 'w')
	except IOError:
		print 'Could not open file to write'
	else:
		for seq in seqDict:
			newFile.write('%s\t%s\n' % (seq, seqDict[seq]))
		print 'counts file for %s created' % (fastQFileName)
		newFile.close()

def createCountFileFasta(fastName):
	fastFile = open(fastName, 'r')
	
	#check to see if there is already a counts file for this file:
	if os.path.isfile(fastName + '.counts'):
		print '%s already exists!  Not overwriting...' % (fastName + '.counts')
		return 0
	
	
	##check where first line is
	for line in fastFile:
		if line.startswith('>') or line.startswith('@') or line.startswith(';'):
			break
	
	runningSeq = ''
	seqDict = {}
	for line in fastFile:
		if line.startswith('>'):
			#save sequence in dict
			seq = runningSeq
			runningSeq = '' #Reset running seq for next entry
			if seq in seqDict:
				#print 'inside'
				seqDict[seq] = seqDict[seq] + 1
			else:
				seqDict[seq] = 1			
		else:
			#add line to running sequence
			runningSeq = runningSeq + line.strip()
	fastFile.close()
	
	#Now save it into a count file
	try:
		newFile = open(fastName + '.counts', 'w')
	except IOError:
		print 'Could not open file to write'
	else:
		for seq in seqDict:
			newFile.write('%s\t%s\n' % (seq, seqDict[seq]))
		print 'counts file for %s created' % (fastName)
		newFile.close()

def createCountFileSL(slFileName):
	slFile = open(slFileName, 'r')
	
	#check to see if there is already a counts file for this file:
	if os.path.isfile(slFileName + '.counts'):
		print '%s already exists!  Not overwriting...' % (slFileName + '.counts')
		return 0
	
	seqDict = {}
	for line in slFile:
		seq = line.strip()
		#print seq
		
		if seq in seqDict:
			seqDict[seq] = seqDict[seq] + 1
		else:
			seqDict[seq] = 1
				
	slFile.close()
	
	#Now save it into a count file
	try:
		newFile = open(slFileName + '.counts', 'w')
	except IOError:
		print 'Could not open file to write'
	else:
		for seq in seqDict:
			newFile.write('%s\t%s\n' % (seq, seqDict[seq]))
		print 'counts file for %s created' % (slFileName)
		newFile.close()


			
