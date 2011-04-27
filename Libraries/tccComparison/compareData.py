import math
import bioLibCG


def createStep(sortedList, complexity):
	'''Might want to check if the complexity is greater than the list size
	Given the complexity (or number of bins requested), devise the best spacer or "step"'''
	step = int(math.ceil(float((sortedList[-1] + 1 - sortedList[0]))/complexity))
	return step
	
def getStepIndex(givenNumber, minNum, stepSize):
	'''Given a number and index parameters, return the # of the bin that the given
	number belongs to for index retrieval purposes'''
	spacer = math.floor(float((int(givenNumber) - int(minNum)))/stepSize)
	index = int(int(minNum) + spacer*stepSize)
	
	return index
	
def returnTccIndex(tccList, complexity):
	'''Given a list of tcc coordinates and the number of bins, return an index dictionary 
	containing all coordinates placed in their respective bins.  The fxn also returns the 
	initial bin key and the "step," or space between each bin which is used to cycle through 
	the bins for index retrieval'''
	sortedList = [int(x.strip().split(':')[2]) for x in tccList]
	sortedList.extend([int(x.strip().split(':')[3]) for x in tccList])#last coords
	sortedList.sort()
	#print sortedList
	step = createStep(sortedList, complexity)
	#print step

	indexDict = {}
	for coord in tccList:
		#print coord
		numCheckLow = int(coord.strip().split(':')[2])
		numCheckHigh = int(coord.strip().split(':')[3])
		indexCheckLow = getStepIndex(numCheckLow, sortedList[0], step)
		indexCheckHigh = getStepIndex(numCheckHigh, sortedList[0], step)

		indexChecks = range(indexCheckLow, indexCheckHigh + 1, step)

		#print indexChecks
		for indexCheck in indexChecks:
			if indexCheck in indexDict:
				indexDict[indexCheck].append(coord)
			else:
				indexDict[indexCheck] = [coord]

	return (indexDict, sortedList[0], step)
			
def tccFileToList(tccFileName, column):
	'''creates list of tcc coordinates from tab sep file'''
	tccFile = open(tccFileName, 'r')
	
	tccList = [line.strip().split('\t')[column] for line in tccFile]
	
	return tccList
			
def compareTwoTcc(tccListOne, tccListTwo, order = 0, complexity = None):
	'''Checks overlaps between two TCC lists
	Indexed --> Runs quicker.
	Complexity is roughly # of bins --> defaults to length of shortest list
	Returns list with overlapping sequences from BOTH lists --> think shared space of vinn diagram
	Takes longer than returning single list --> see other function
	ORDER DECIDES WHICH LIST IS TO BE INDEX --> NON INDEXED LIST IS RETURNED
	0 = shortest is indexed, 1 = first list is returned, 2 = second list is returned'''
	
	#decide which list is to be indexed
	if order == 0: #shortest list is indexed
		if len(tccListOne) < len(tccListTwo):
			indexList = tccListOne
			otherList = tccListTwo
		else:
			indexList = tccListTwo
			otherList = tccListOne
	elif order == 1:
		indexList = tccListTwo
		otherList = tccListOne
	elif order == 2:
		indexList = tccListOne
		otherList = tccListTwo
		
	#make sure complexity is set, should set max complexity here...
	if complexity == None:
		complexity = len(indexList)
	#make the index
	(tccIndex, indexMin, indexStep) = returnTccIndex(indexList, complexity)
	#print tccIndex, complexity
	
	#compare vs sequences
	reducedList = []
	for coord in otherList:
		#print coord
		numCheckLow = int(coord.strip().split(':')[2])
		numCheckHigh = int(coord.strip().split(':')[3])
		indexCheckLow = getStepIndex(numCheckLow, indexMin, indexStep)
		indexCheckHigh = getStepIndex(numCheckHigh, indexMin, indexStep)
		indexChecks = range(indexCheckLow, indexCheckHigh + 1, indexStep)
		#print '  ',numCheckLow, indexCheckLow
		#print '  ',numCheckHigh, indexCheckHigh
		
		for indexCheck in indexChecks:
			if coord in reducedList: break #Already overlapped from another indexCheck...
			if indexCheck in tccIndex:
				#print 'Got Here'
				for indexCoord in tccIndex[indexCheck]:
					#print 'Checking Overlap', indexCoord, coord
					if bioLibCG.tccOverlap(coord, indexCoord):
						#print 'overlapped',indexCoord,coord
						reducedList.append(coord)
						break #already found this coords overlap --> don't check other indexCoords
	return reducedList

def compareTwoTccALL(tccList1, tccList2, complexity = None):
	overlaps = compareTwoTcc(tccList1, tccList2, 1, complexity) #get overlaps from first data set
	overlaps.extend(compareTwoTcc(tccList1, tccList2, 2, complexity)) #"" second data set
	
	return overlaps
	
		
def makeConnectionsDict(tccList, complexity = None):
	'''A connections dictionary gives each coordinates connections to other coordinates in the list
	tcc : [tcc1, tcc2, etc]'''
	#make sure complexity is set
	if complexity == None:
		complexity = len(tccList)

	#make the index
	(tccIndex, indexMin, indexStep) = returnTccIndex(tccList, complexity)

	#compare vs sequences
	connectionsDict = {} # format is coord: [tcc, tcc, etc].  It means X coord is connected to other coords
	for coord in tccList:
		#Add coord to connections dict --> every coord is a key
		connectionsDict[coord] = []
		
		numCheckLow = int(coord.strip().split(':')[2])
		numCheckHigh = int(coord.strip().split(':')[3])
		indexCheckLow = getStepIndex(numCheckLow, indexMin, indexStep)
		indexCheckHigh = getStepIndex(numCheckHigh, indexMin, indexStep)
		indexChecks = range(indexCheckLow, indexCheckHigh + 1, indexStep)

		for indexCheck in indexChecks:
			if indexCheck in tccIndex:
				for indexCoord in tccIndex[indexCheck]:
					if indexCoord == coord:
						pass
					elif bioLibCG.tccOverlap(coord, indexCoord):
						#print 'overlapped',indexCoord,coord
						if indexCoord not in connectionsDict[coord]:
							connectionsDict[coord].append(indexCoord)
						#break #wouldn't breaking here prevent finding all overlaps?
	return connectionsDict
