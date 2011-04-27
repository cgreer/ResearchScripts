#get results that are only noncoding

import bioLibCG as cg
import compareData as compare
predName = '/home/chrisgre/projects/NoncodingMouse/results/NCmouse-s3k8b17.bothNCandC.results'
keepList = compare.tccFileToList('keepNoncoding.tcc', 0)
predList = compare.tccFileToList(predName, 1)

keepers = compare.compareTwoTcc(predList, keepList, 1)
print len(keepers)

#now go back through pred file and create a new file with only lines that have noncoding in them

predFile = open(predName, 'r')
outFile = open('NCmouse.noncoding.results', 'w')

predLines = predFile.readlines()
predFile.close()
newLines = {}
for keeper in keepers:
	for line in predLines:
		if keeper in line:
			newLines[line] = 1

for line in newLines:
	outFile.write(line)
