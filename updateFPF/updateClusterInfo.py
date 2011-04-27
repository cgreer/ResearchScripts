import bioLibCG as cg
import cgConfig

def updateClusterInfo(cName = None):
	#init
	conf = cgConfig.getConfig(cName)
	
	pFileName = conf.conf['results']
	sortedClustersFileName = conf.conf['sortedClusters']
		
	#add cluster IDs to tccs
	sortedFile = open(sortedClustersFileName, 'r')
	idDict = {} # format:  tcc : clusterID
	i = 0
	for line in sortedFile:
		clusterID = str(i)
		i = i + 1
		
		#put cluster into list
		cluster = line.strip().split(',')
		cluster.remove('') #pesky last comma
		
		for tcc in cluster:
			idDict[tcc] = clusterID
	sortedFile.close()
	
	#get tccs that are overlapping with known sequences
	overDict = {}
	predFile = open(pFileName, 'r')
	for line in predFile:
		if line.strip().split('\t')[6] == '1':
			clusterID = idDict[line.strip().split('\t')[1]]
			overDict[clusterID] = 1
	predFile.close()
		
	#now remake file
	newLines = []
	predFile = open(pFileName, 'r')
	for line in predFile:
		clusterID = idDict[line.strip().split('\t')[1]]
		if clusterID in overDict:
			newLine = line.strip().split('\t')
			newLine[7] = str(clusterID)
			newLine[8] = str(1)
			newLine = '\t'.join(newLine) + '\n'
			newLines.append(newLine)
		else:
			newLine = line.strip().split('\t')
			newLine[7] = str(clusterID)
			newLine[8] = str(0)
			newLine = '\t'.join(newLine) + '\n'
			newLines.append(newLine)
	predFile.close()
	
	#write new File
	newFile = open(pFileName, 'w')
	for line in newLines:
		newFile.write(line)
	
if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		updateClusterInfo(sys.argv[1])
	else:
		updateClusterInfo()

