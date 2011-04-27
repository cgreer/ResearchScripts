#puts final hits into clusters...
##Clusters are based off of overlapping neighbors, if you have an overlapping neighbor than you are part of that cluster.
import bioLibCG as cg
import subprocess
import compareData as compare
import cgConfig


#Start Timer
timer = cg.cgTimer()
timer.start()

#Get list of mature tccs
conf = cgConfig.returnConfDict()
finalMirFileName = '/u/home8/gxxiao/chrisgre/projects/PipeRuns/LanderHuman/out/LanderHuman-s3k8b17.ALL.FINAL.mirs.tsv'
finalMirFileName = conf['resultsRaw']
matureTccs = compare.tccFileToList(finalMirFileName, 1) # list of all mature micro in tcc
print 'List getting', timer.split()


#make connections dict
matureConnections = compare.makeConnectionsDict(matureTccs)
print 'Make connections:', timer.split()

#Now have to define Clusters...
clusters = []
addedList = []

#I don't think python passes by reference? also I think this function is in the middle because it uses a global variable :P
def createClusters(item = None, mode = None):
		
	if item in addedList:
		return 0
	elif mode == "top":
		clusters.append([item])
		addedList.append(item) ##creates new cluster with the item already stored in it
		for connectedItem in matureConnections[item]:
			createClusters(connectedItem, "neighbor")
	elif mode == "neighbor":
		clusters[-1].append(item) #add this item to the last cluster created
		addedList.append(item)
		for connectedItem in matureConnections[item]:
			createClusters(connectedItem, "neighbor")
	
for tcc in matureTccs:
	createClusters(tcc, "top")

print 'Make Clusters', timer.split()


#Sort Clusters.
sortedClusters = []

for cluster in clusters:
	sortedClusters.append(cg.sortTccList(cluster))

print 'Sort Clusters:', timer.split()


#Output sorted cluster file
clusterFileName = conf['sortedClusters']
clusterFile = open(clusterFileName, 'w')
for cluster in sortedClusters:
	for hit in cluster:
		clusterFile.write('%s,' % hit)
	clusterFile.write('\n')
clusterFile.close()

'''
#re-create sortedClusters list:
clusterFileName = 'sortedClusters.data'
clusterFile = open(clusterFileName, 'r')
sortedClusters = []


for line in clusterFile:
	sortedClusters.append([])
	line = line.strip()[0:-1] #take off last comma ;P
	for hit in (line.strip().split(',')):
		sortedClusters[-1].append(hit)
'''


print 'Store intermediate data:', timer.split()


#output hitsAround file
outputFile = open(conf['hitsPerFrame'], 'w')
outputFile.write("FrameCoord\thits\n")

frameLength = 200
frameShift = 1
for cluster in sortedClusters:
	#grab first and last coordinate from cluster, for each cluster deduce how many theoretical microRNAs were in hitScope
	clusterChrom = cluster[0].split(":")[0]
	clusterStrand = cluster[0].split(":")[1]
	firstCoord = int(cluster[0].split(":")[2])
	#print cluster[-1]
	lastCoord = int(cluster[-1].split(":")[3])
	
	
	startCoord = firstCoord
	while startCoord < lastCoord:
		#count how many hits there are in this range
		rangeStart = startCoord - (frameLength/2)
		rangeEnd = startCoord + (frameLength/2)
		rangeTcc = '%s:%s:%s:%s' % (clusterChrom, clusterStrand, rangeStart, rangeEnd)
		overlappedList = compare.compareTwoTcc([rangeTcc], cluster, 2)
		hitCount = len(overlappedList) 
		
		#output 
		outputFile.write('%s\t%s\n' % (rangeTcc, hitCount))
		startCoord = startCoord + frameShift #check overlap with range
outputFile.close()

print 'Output Hits per Frame:', timer.split()
print 'Overall Time:', timer.report()


