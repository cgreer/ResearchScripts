import cgLuckyCharmsFlat as cgLuckyCharms
from copy import copy
import bioLibCG
import cgIndex

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
	with open(fn) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def getIDRange(paraInfo, fn):
	'''numJobs and job number are 1 based'''

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

class Nexus:

	def __init__(self, dataFileName, dataClass):
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
                self._conditions = {}

	def bindAttribute(self, attributeName):
		
		#bind data to class attribute for easy access
		setattr(self, attributeName, self._attName_id_value[attributeName])
	
	def loadTranscriptionInfo(self, attNames):
		'''loads caste fxns, column positions, default values for each selected attribute'''		

		for attName in attNames:
			dataField = getattr(self._dataClass, attName)
			self._attName_casteFromFxn[attName] = cgLuckyCharms.getCasteFunction(dataField.dataType, True)
			self._attName_casteToFxn[attName] = cgLuckyCharms.getCasteFunction(dataField.dataType, False)
			self._attName_columnPosition[attName] = dataField.dataSlot
			self._attName_defaultValue[attName] = dataField.dataDefault

	def load(self, attNames, paraInfo = [None, None], idRange = [], conditions = {}):
                '''paraInfo is [runNumber, numberOfRuns].  First parallel is checked, and then idRange'''
                self._conditions = conditions

		#if running a parallel job, split it into the right ids...
		if paraInfo != [None, None]:
			idRange = getIDRange(paraInfo, self._dataFileName)
	
		#if running parallel or specific range, mark range info
		self._selectedAttNames = attNames		
		if idRange:
			self._rangeSpecified = [idRange[0], idRange[-1]]
			
		#get casting and column info
		self.loadTranscriptionInfo(attNames)

		#initialize master dict
		for attName in attNames:
			self._attName_id_value[attName] = {}

		#get number of slots
		f = open(self._dataFileName, 'r')
		numSlots = len(f.readline().split('\t'))
		f.close()
		
                loadTime = 0.0
                stripTime = 0.0
                idTime = 0.0
                tranTime = 0.0
                conditionTime = 0.0

                #skip to start of specified range
                if self._rangeSpecified:
                        fIndex = cgIndex.lineIndex(self._dataFileName)
                        fIndex.passCheckFunction(cgIndex.primaryIDCheckFunction)
                        fIndex.binarySearch(self._rangeSpecified[0])
                        f = fIndex.file
                else:
                        f = open(self._dataFileName, 'r')
	        

                #transcribe values
                for line in f:

			ls = line.strip().split('\t')
			id = int(ls[0]) #id is always first slot
		

			#only transcribe selected range!
			if idRange:
				if id > idRange[1]:
		                        break

			#transcribe
			for attName in attNames:
				if self._attName_columnPosition[attName] < numSlots:
                                        if ls[self._attName_columnPosition[attName]] != '.':
                                                self._attName_id_value[attName][id] = self._attName_casteFromFxn[attName](ls[self._attName_columnPosition[attName]])
                                        else:
					        self._attName_id_value[attName][id] = copy(self._attName_defaultValue[attName])
				else:
					self._attName_id_value[attName][id] = copy(self._attName_defaultValue[attName])

                        #do conditions
                        if conditions:
                                for attName in conditions:
                                        if attName == 'ID':
                                                if conditions['ID'](id):
                                                        self._selectedIDs.add(id)
                                                else:            
                                                        for aName in attNames:
                                                                del self._attName_id_value[aName][id]
                                                
                                        else:
                                                if conditions[attName](self._attName_id_value[attName][id]):
                                                        self._selectedIDs.add(id)
                                                else:            
                                                        for aName in attNames:
                                                                del self._attName_id_value[aName][id]
                f.close()

                

		#bind attribute names to dictionaries
		for attName in attNames:
			self.bindAttribute(attName)

                print 'done loading', self._dataFileName


	def save(self, outFN = None):
		
		if outFN == None: outFN = self._dataFileName
		if self._rangeSpecified:
			outFN += '.range.%s.%s' % (self._rangeSpecified[0], self._rangeSpecified[1]) 
                
                
                #skip to start of specified range
                if self._rangeSpecified:
                        fIndex = cgIndex.lineIndex(self._dataFileName)
                        fIndex.passCheckFunction(cgIndex.primaryIDCheckFunction)
                        fIndex.binarySearch(self._rangeSpecified[0])
                        f = fIndex.file
                else:
                        f = open(self._dataFileName, 'r')
		
                #create new file contents
		newLines = []
		for line in f:
			ls = line.strip().split('\t')
			id = int(ls[0])
                        
                        #skip those that weren't selected and in range
                       
                        ''' #I'm not sure I want to only write the ones with conditions...
                        if self._conditions:
                                if id not in self._selectedIDs: continue
			'''

                        #stop checking for ids once out of range
                        if self._rangeSpecified:
				if id > self._rangeSpecified[1]: break
                       
                        #save the rest
			for attName in self._selectedAttNames:
				newVal = self._attName_casteToFxn[attName](self._attName_id_value[attName][id])
				ls = lineUpdate(ls, newVal, self._attName_columnPosition[attName])
			
			newLines.append('%s\n' % '\t'.join(ls))
		f.close()

		#output file
		f = open(outFN, 'w')
		f.writelines(newLines)
		f.close()

		#exit signal for parallel processes
                if self._rangeSpecified:
                        f = open(outFN + '.exitSignal', 'w')
                        f.write('DONE')
                        f.close()

