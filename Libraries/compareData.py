import math
import bioLibCG
import cgClusters
import copy

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
	
	tccFile.close()
	
	return tccList
			
def compareTwoTcc(tccListOne, tccListTwo, order = 0, complexity = None, amount = False):
	'''Checks overlaps between two TCC lists
	Indexed --> Runs quicker.
	Complexity is roughly # of bins --> defaults to length of shortest list
	Returns list with overlapping sequences from BOTH lists --> think shared space of vinn diagram
	Takes longer than returning single list --> see other function
	ORDER DECIDES WHICH LIST IS TO BE INDEX --> NON INDEXED LIST IS RETURNED
	0 = shortest is indexed, 1 = first list is returned, 2 = second list is returned
	amount refers to if you want to return the amount of overlap for each transcript [coord, amount]'''
	
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
					overlap = bioLibCG.tccOverlap(coord, indexCoord, True)
					if overlap:#if there is any overlap...
						#print 'overlapped',indexCoord,coord
						if amount:
							reducedList.append([coord, overlap])
						else:
							reducedList.append(coord)
						break #already found this coords overlap --> don't check other indexCoords
	return reducedList

def getIndividualOverlaps(tccListOne, tccListTwo, order, complexity = None):
	'''Checks overlaps between two TCC lists
	Indexed --> Runs quicker.
	Complexity is roughly # of bins --> defaults to length of shortest list
	Returns list with overlapping sequences from BOTH lists --> think shared space of vinn diagram
	Takes longer than returning single list --> see other function
	0 = shortest is indexed, 1 = first list is returned, 2 = second list is returned'''
	
	#decide which list is to be indexed
	if order == 1:
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
        individualOverlaps = {} # tcc : [tcc, tcc]
	for coord in otherList:

                individualOverlaps[coord] = []

		#print coord
		numCheckLow = int(coord.strip().split(':')[2])
		numCheckHigh = int(coord.strip().split(':')[3])
		indexCheckLow = getStepIndex(numCheckLow, indexMin, indexStep)
		indexCheckHigh = getStepIndex(numCheckHigh, indexMin, indexStep)
		indexChecks = range(indexCheckLow, indexCheckHigh + 1, indexStep)
		#print '  ',numCheckLow, indexCheckLow
		#print '  ',numCheckHigh, indexCheckHigh
		
		for indexCheck in indexChecks:
			if indexCheck in tccIndex:
				for indexCoord in tccIndex[indexCheck]:
					overlap = bioLibCG.tccOverlap(coord, indexCoord, True)
					if overlap:#if there is any overlap...
                                                if indexCoord not in individualOverlaps[coord]:
                                                        individualOverlaps[coord].append(indexCoord)
        return individualOverlaps                     


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

def unionTwoTccLists(tccListOne, tccListTwo):
	overlapping = compareTwoTccALL(tccListOne, tccListTwo)
	
	print len(tccListOne), len(tccListTwo), len(overlapping)
	#Take out overlapping from both tccLists
	for tcc in overlapping:
		if tcc in tccListOne:
			tccListOne.remove(tcc)
		if tcc in tccListTwo:
			tccListTwo.remove(tcc)
	
	print len(tccListOne), len(tccListTwo), len(overlapping)
	
	clusters = cgClusters.cgClusters()
	clusters.makeClusters(overlapping)
	
	
	newSeqs = []
	for cluster in clusters.getClusters():
		newSeqs.append(cgClusters.returnTccListRange(cluster))
	print len(tccListOne), len(tccListTwo), len(overlapping), len(newSeqs)
	
	newSeqs.extend(tccListOne)
	newSeqs.extend(tccListTwo)
	
	print len(tccListOne), len(tccListTwo), len(overlapping), len(newSeqs)
	
	return newSeqs
		
def checkIfOverlaps(tccList):
	'''Check if there are any sequences in the list that overlap each other.'''
	if not len(tccList) > 0:
		return False
	
	connectionsDict = makeConnectionsDict(tccList)	
	
	oFlag = False
	for tcc in connectionsDict:
		if len(connectionsDict[tcc]) > 1:
			oFlag = True
			break
	
	return oFlag

def collapseOverlaps(tccList):
	'''If the list has overlapping sequences, combine the sequences that overlap.
	This makes a cluster and then returns the start and end point of the cluster'''
		
	clusters = cgClusters.cgClusters()
	clusters.makeClusters(tccList)
	
	newSeqs = []
	for cluster in clusters.getClusters():
		newSeqs.append(cgClusters.returnTccListRange(cluster))
	
	return newSeqs

def filterOutTccs(tccList, tccDirectory, getFilteredOut = False):
	'''Given a directory of known tcc files, filter out sequences from a tcc list that
	overlap with any sequences in the known directory'''
	overlapped = []
	for filename in bioLibCG.getDirFiles(tccDirectory, end = '.tcc'):
		tccListKnown = tccFileToList(filename, 0)
		for tcc in compareTwoTcc(tccListKnown, tccList, 2):
			if tcc not in overlapped:
				overlapped.append(tcc)
	
	if getFilteredOut:
		return overlapped
	else:
		for tcc in overlapped:
			tccList.remove(tcc)
		return tccList

def subtractTwoTcc(keeperTcc, otherTcc):
	'''returns list of tcc coordinates after subtracting the second from the first.'''
	chromK = keeperTcc.split(':')[0]
	strandK = keeperTcc.split(':')[1]
	startK = keeperTcc.split(':')[2]
	endK = keeperTcc.split(':')[3]
	
	chromO = otherTcc.split(':')[0]
	strandO = otherTcc.split(':')[1]
	startO = otherTcc.split(':')[2]
	endO = otherTcc.split(':')[3]
	
	if (chromK == chromO) and (strandK == strandO):
		
		#put them all in list
		linearList = []
		linearList.extend([int(startK), int(endK)])
		if int(startO) not in linearList: linearList.append(int(startO))
		if int(endO) not in linearList: linearList.append(int(endO))
		linearList.sort()
		
		if startO == startK:
			same = True
		else:
			same = False
			
		#initalize if keeping or cutting first.
		if same:
			i = 0
			keep = False
		elif linearList[0] == int(startK):
			i = 0
			keep = True
			#print 'keep'
		else:
			i = 1 #if 0 wasn't the known start, 1 is
			keep = False
			#print 'cut'
		
		#print linearList
		
		returnList = []
		while True:
			if keep: #keep this section of coords
				newStart = linearList[i]
				newEnd = linearList[i + 1]
				if not str(newStart) == startK:# Isn't initial beginning, add one
					newStart = newStart + 1
				if str(newEnd) != endK:# the end isn't the final end, subtract one
					newEnd = newEnd - 1
					
				returnList.append('%s:%s:%s:%s' % (chromK, strandK, newStart, newEnd))
				if str(linearList[i+1]) == endK:
					break
				keep = False
				
			else: #skip this section of coords
				if str(linearList[i+1]) == endK:
					break
				keep = True #flip the keep/cut switch
			
			i = i + 1
			
		return returnList

def recurseSubtract(sList, otherList):
	'''When a tcc is subtracted, it's subtraction must be subtracted 
	with all the other tcc's in otherlist.  The best way to do this is
	to break it up in a recursive function'''
	#print 'here recurse'
	
	totalOverlap = False
	subtractList = []
	
	#print 'Lists (keep, other):', sList, otherList
	
	for tccK in sList: #
		overlap = False
		#print 'known check', tccK
		for tccO in otherList:
			#print ' other check', tccO
			if bioLibCG.tccOverlap(tccK, tccO):
				totalOverlap = True
				overlap = True
				for sTcc in subtractTwoTcc(tccK, tccO):
					#print '  adding', sTcc, 'to list'
					if sTcc not in subtractList: subtractList.append(sTcc)
					#print '  list', subtractList
				break
		if not overlap:
			#print 'tccK did not overlap anything, adding to list'
			subtractList.append(tccK)
			#print subtractList
	
	if totalOverlap:
		rList = recurseSubtract(subtractList, otherList)
		return rList #once all other recursions have ended return the final list
	else:
		return subtractList #no overlaps, return all the way to top...
		
	

def subtractTwoTccLists(tccListKeep, tccListOther):
	#Both Lists should already be COLLAPSED!!!
	
	if checkIfOverlaps(tccListKeep): print 'THE KEEPER LIST HAS OVERLAPS (SUBTRACTION)'
	if checkIfOverlaps(tccListOther): print 'THE OTHER LIST HAS OVERLAPS (SUBTRACTION)'
	
	overlapKeep = compareTwoTcc(tccListKeep, tccListOther, 1)
	overlapOther = compareTwoTcc(tccListKeep, tccListOther, 2)	
	
	subList = recurseSubtract(overlapKeep, overlapOther)
	#print 'subList:',  subList
	#for those that didn't overlap return them
	
	newSeqs = []# make a new list so as to not overwrite the other.
	newSeqs.extend(tccListKeep)
	
	for tcc in overlapKeep:
		newSeqs.remove(tcc)
		
	newSeqs.extend(subList)
	
	return newSeqs
				
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getIndividualOverlaps, sys.argv)

