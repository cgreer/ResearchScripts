'''
This Script Appends cluster values for each prediction produced by mirPipe

mirPipe --> defineClusters --> THIS SCRIPT --> (output) final mirs + cvals

foreach predicted:
	Check which blocks it overlaps with
	Get highest cVal for predicted.
	Append Final mirs file to add new cVals for each prediction
'''
import bioLibCG
import cgConfig

def updateDensity(cName = None):
	#Create hitmap for blocks, cValdict for block
	conf = cgConfig.getConfig(cName)
	blockFileName = conf.conf['hitsPerFrame']# created in defineCluster script folder
	blockFile = open(blockFileName, 'r')
	blocksList = []
	cValBlockDict = {}
	
	for line in blockFile:
		blocksList.append(line.strip().split('\t')[0])
		cValBlockDict[line.strip().split('\t')[0]] = int(line.strip().split('\t')[1])
	blockFile.close()
	blockHitmap = bioLibCG.createHitMap(blocksList)
	
	#Now append the cVal for each predicted line:
	
	predictedFileName = conf.conf['results']
	predictedFile = open(predictedFileName, 'r')
	
	newFileList = []
	counter = 0
	for line in predictedFile:
		counter = counter + 1
		#print counter
		cVal = 0
		#what blocks does this prediction overlap?
		tccPrediction = line.strip().split('\t')[1] #This should be mature?
		coordsPrediction = bioLibCG.stripTripleColon(tccPrediction)
		for i in range(int(coordsPrediction['start']), int(coordsPrediction['end'])):
			if i in blockHitmap:
				for block in blockHitmap[i]:
					if bioLibCG.tccOverlap(tccPrediction, block):
						if cValBlockDict[block] > cVal:
							cVal = cValBlockDict[block]
		newLine = line.strip().split('\t')
		newLine[5] = str(cVal)
		newLine = '\t'.join(newLine) + '\n'
		newFileList.append(newLine)
	predictedFile.close()
	
	
	newFileName = conf.conf['results']
	newFile = open(newFileName, 'w')
	for line in newFileList:
		newFile.write(line)
	
	newFile.close()

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		updateDensity(sys.argv[1])
	else:
		updateDensity()
	



