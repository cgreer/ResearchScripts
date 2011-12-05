import cgLuckyCharmsFlat as cgLuckyCharms
from copy import copy
import bioLibCG
import cgIndex
import cgFile

def lineUpdate(lineData, data, position):
	'''lineList must NOT contain CR.  data must be string'''

	numSlots = len(lineData)

	#put data in right position
	if position < numSlots:
		lineData[position] = data
	else:
		#update lines that don't exist/have no values
		for i in range(numSlots, position):
			lineData.append('.')
		lineData.append(data)
	
	return lineData

def getNumFileLines(fn):
        print 'getting num file lines'
        for i, l in enumerate(open(fn)):
            pass
	return i + 1

def getIDRange(paraInfo, fn):
	'''numJobs and job number are 1 based'''
        print 'getting ID range (in fxn)'
	numJobs = int(paraInfo[1])
	n = int(paraInfo[0])
	numIDs = getNumFileLines(fn) #one id per line
	sectorLength = float(numIDs)/numJobs #plus one will not work on small numbers

	rStart = int((n-1)*sectorLength) + 1
	rEnd = int((n)*sectorLength)
	if n == numJobs: rEnd = numIDs - 1
	if n == 1: rStart = 0

	return [rStart, rEnd]	

class Field:

	def __init__(self, dataType, dataDefault, dataSlot):
		self.dataType = dataType
                self.dataDefault = dataDefault #prevents list aggregation
		self.dataSlot = dataSlot

def shellNexus(NX, dataFileName, select):

    newNX = Nexus(dataFileName, NX._dataClass, NX.hasIDs)
    newNX._attName_id_value = {}
    newNX._selectedAttNames = select
    newNX.loadTranscriptionInfo(select)
    newNX._rangeSpecified = copy(NX._rangeSpecified)
    newNX._selectedIDs = set()
    newNX._packetInfo = copy(NX._packetInfo)
    newNX._splitRunFlag = copy(NX._splitRunFlag)
    newNX.hasIDs = NX.hasIDs

    return newNX

class Nexus:

	def __init__(self, dataFileName, dataClass, ids = True):
		self._dataFileName = dataFileName
                self._dataClass = dataClass
		self._attName_id_value = {}
		self._attName_casteFromFxn = {}
		self._attName_casteToFxn = {}
		self._attName_columnPosition = {}
		self._attName_defaultValue = {}
		self._selectedAttNames = []
		self._rangeSpecified = []
	        self._selectedIDs = set()
                self._packetInfo = None
                self._splitRunFlag = False
                self.hasIDs = ids
        
        def __str__(self):
            newLines = []
            for id in self.ids:
                newLine = [str(id)]
                
                for attName in self._attName_id_value:
                    newLine.append(str(self._attName_id_value[attName][id]))
                newLines.append('\t'.join(newLine))    

            return '\n'.join(newLines)

        def collapseColumnNumbers(self, attNames):
                
                #make column numbers 0...n
                for i, attName in enumerate(attNames):
                    self._attName_columnPosition[attName] = i
    
	def bindAttributes(self, attributeNames):
		
		#bind data to class attribute for easy access
                for attributeName in attributeNames:
                    setattr(self, attributeName, self._attName_id_value[attributeName])
	
	def loadTranscriptionInfo(self, attNames):
		'''loads caste fxns, column positions, default values for each selected attribute'''		

		for attName in attNames:
			dataField = getattr(self._dataClass, attName)
			self._attName_casteFromFxn[attName] = cgLuckyCharms.getCasteFunction(dataField.dataType, True)
			self._attName_casteToFxn[attName] = cgLuckyCharms.getCasteFunction(dataField.dataType, False)
			self._attName_columnPosition[attName] = dataField.dataSlot
			self._attName_defaultValue[attName] = dataField.dataDefault


        def initializeMasterDict(self):
                #initialize master dict
                for attName in self._selectedAttNames:
                        self._attName_id_value[attName] = {}

        def getNumberOfSlots(self):
                try:
                    f = open(self._dataFileName, 'r')
                    numSlots = len(f.readline().split('\t'))
                    f.close()
                    return numSlots
                except IOError:
                    return 0

        def linkIDsToColumn(self):
                self.ids = eval('self.%s' % self._selectedAttNames[0])

	def load(self, attNames, paraInfo = [None, None]):
                '''paraInfo is [runNumber, numberOfRuns]'''
       
                #t = bioLibCG.cgTimer()
                #stage_cumTime = dict( (x, 0.0) for x in (''))
                #t.start()

                if paraInfo == ['splitRun', 'splitRun']:
                        self._splitRunFlag = True
                        paraInfo = [None, None] # now treat paraInfo as if there was nothing...
                
                if paraInfo != [None, None]: 
                        paraInfo[0] = int(paraInfo[0])
                        paraInfo[1] = int(paraInfo[1])
                        self._packetInfo = cgFile.getPacketInfo(self._dataFileName, paraInfo[1])[paraInfo[0] - 1]
                        
		#if running parallel or specific range, mark range info
		self._selectedAttNames = attNames		
		

		#get casting and column info
		self.loadTranscriptionInfo(attNames)

                #init master dictionaries
                self.initializeMasterDict()

		#get number of slots
                numSlots = self.getNumberOfSlots()
		
                #open file and binary skip to correct line if packet            
                dataFile = cgFile.cgFile(self._dataFileName)
                if self._packetInfo:
                        dataFile.seekToLineStart(self._packetInfo[0])

                #print 'before loop', t.split()
                #transcribe values
                currentID = 0
                for line in dataFile.file:

			ls = line.strip().split('\t')

                        #get ID
                        if self.hasIDs:
                            id = int(ls[0]) #id is always first slot
                        else:
                            id = currentID
                            currentID += 1
		
                        #stop if at end of range
                        if self._packetInfo:
                                if id == self._packetInfo[1]:
                                        break

			#transcribe
                        #Note lots of copying is SLOW (10x)
                        #only copy if list?
			for attName in attNames:
                                colPosition = self._attName_columnPosition[attName]
				if colPosition < numSlots:
                                        if ls[colPosition] != '.':
                                                self._attName_id_value[attName][id] = self._attName_casteFromFxn[attName](ls[colPosition])
                                        else:
					        self._attName_id_value[attName][id] = copy(self._attName_defaultValue[attName])
				else:
					self._attName_id_value[attName][id] = copy(self._attName_defaultValue[attName])
                #print 'after loop', t.split()
                dataFile.file.close()
                
		#bind attribute names to dictionaries
                self.bindAttributes(attNames)

                #bind id attribute to first attribute, they all have the same ids...
                self.linkIDsToColumn()
                #print 'finishing stuff', t.split()

        def save(self, outFN = None):
		
		if outFN == None: outFN = self._dataFileName

                if self._packetInfo:
			outFN += '.range.%s.%s' % (self._packetInfo[0], self._packetInfo[1]) 

                
                dataFile = cgFile.cgFile(self._dataFileName)
                if self._packetInfo:
                        dataFile.seekToLineStart(self._packetInfo[0])
                
                #create new file contents
                currentID = 0
		newLines = []
                for line in dataFile.file:
			ls = line.strip().split('\t')
			
                        if self.hasIDs:
                            id = int(ls[0])
                        else:
                            id = currentID
                            currentID += 1
                        
                        if self._packetInfo:
                                if id == self._packetInfo[1]: break

                        #save the rest
                        #TODO: lineUpdate with multiple injections
			for attName in self._selectedAttNames:
				newVal = self._attName_casteToFxn[attName](self._attName_id_value[attName][id])
				ls = lineUpdate(ls, newVal, self._attName_columnPosition[attName])

                        #only one newLine no matter the amount of attributes updated	
			newLines.append('%s\n' % '\t'.join(ls))
                dataFile.file.close()

		#output file
                newLines = ''.join(newLines) #might cause less clogging if there is only one write operation...
		f = open(outFN, 'w')
		f.write(newLines)
		f.close()

		#exit signal for parallel processes
                if self._packetInfo or self._splitRunFlag:
                        f = open(outFN + '.exitSignal', 'w')
                        f.write('DONE')
                        f.close()

def quickTable(*args):
    '''make a class for nexus to use as table template
    on the fly
    args are 4tuples (attName, type, defVal, colPosition)'''
    
    class QuickTable: pass

    for col in args:
        setattr(QuickTable, col[0], Field(col[1], col[2], col[3]))
    
    return QuickTable


