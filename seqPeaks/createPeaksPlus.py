#using the continuos blocks from in the small RNA lib file, identify all the peaks...
import bioLibCG as cg
import cgConfig as c

mConf = c.cgConfig('Main.conf')

fileNames = cg.recurseDir(mConf.conf['wigMouse'],  end = '.wig')

for fN in fileNames:
	file = open(fN, 'r')
	file.readline() #header
	
	#get all points in midpoint form
	pointsDict = {}
	
	for line in file:
		start = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		point = start + (end-start)/2 #midpoint
		pointsDict[point] = int(cg.ss(line)[3].split('.')[0])
	file.close()
	
	
	#determine peaks based off of neighbors of each point
	lowest = pointsDict.keys()
	lowest.sort()
	peaks = []
	span = 2 #must be > 0
	for i in range(span + 1,len(lowest) - span - 1):
			
		val = pointsDict[lowest[i]]
		if val < 5: #minimum
			continue
			
		#grab all values of density around peak
		checkList = []
		j = 1
		while j <= span:
			checkList.extend([pointsDict[lowest[i - j]], pointsDict[lowest[i + j]]])
			j += 1
		
		
		#check if lower than val
		peakFlag = True
		for check in checkList:
			if val > check:
				pass
			else:
				peakFlag = False
			
		#output
		if peakFlag:
			#print checkList, val
			peaks.append('%s\t%s\n' % (lowest[i], val))
		
	
	#output new file with new ending
	outFile = open(fN + '.peaks', 'w')
	outFile.writelines(peaks)
	outFile.close()
