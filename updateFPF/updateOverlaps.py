import bioLibCG as cg
import cgConfig

def updateOverlaps(cName = None):
	#init
	conf = cgConfig.getConfig(cName)
	pFileName = conf.conf['results']
	overlapsFileName = conf.conf['matureOverlaps']
	
	
		
	#put overlapping sequences in list:
	overlaps = []
	overlapFile = open(overlapsFileName, 'r')
	for tcc in overlapFile:
		overlaps.append(tcc.strip())
	overlapFile.close()
	
	#check each line of pred file for overlap, add 1 for overlap and 0 for non
	predFile = open(pFileName, 'r')
	newFileLines = []
	for line in predFile:
		mTcc = line.strip().split('\t')[1]
		if mTcc in overlaps:
			newLine = line.strip().split('\t')
			newLine[6] = str(1)
			newLine = '\t'.join(newLine) + '\n'
			newFileLines.append(newLine)
		else:
			newLine = line.strip().split('\t')
			newLine[6] = str(0)
			newLine = '\t'.join(newLine) + '\n'
			newFileLines.append(newLine)
	predFile.close()
	
	#write new File
	newFile = open(pFileName, 'w')
	for line in newFileLines:
		newFile.write(line)

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		updateOverlaps(sys.argv[1])
	else:
		updateOverlaps()


