import bioLibCG as cg
import os

direc = '/home/chrisgre/apps/projects/small.rna.libs'
metaFileName = direc + '/' + 'small.meta'

##make data of already made file so there aren't any duplicates:
fileDict = {} #filename...
metaFile = open(metaFileName, 'r')

#add new entries
countFiles = cg.recurseDir(direc, end = '.fastq')

for file in countFiles:
	fileName = file.strip().split('/')[-1]
	if len(fileName.split('.')) > 2: #has to specifically end in fastq...
		continue
	fileName = cg.getBaseFileName(file, naked = True)
	dir = file.strip().split('/')[-2]
		
	org = 'NONE'
	if 'human' in dir:
		org = 'human'
	if 'mouse' in dir:
		org = 'mouse'
	if 'pig' in dir:
		org = 'pig'
	if 'dog' in dir:
		org = 'dog'
	if 'rat' in dir:
		org = 'rat'
	if 'zebrafish' in dir:
		org = 'zebrafish'
	
	#if specific meta for this file exists, use it, else use meta...
	adap = 'NONE'
	
	#Check if it was Trimmed...
	
	if os.path.exists(file + '.meta'):
		adapDir = file + '.meta'
	else:
		adapDir = '/'.join(file.strip().split('/')[:-1]) + '/adapter.meta'
	#print adapDir
	if os.path.isfile(adapDir):
		adapFile = open(adapDir, 'r')
		adap = adapFile.readline().strip()
		adapFile.close()
	
	#tissue
	tissueList = ['cancer', 'embryo', 'tumor', 'THP', 'H1', 'hela', 'Bcell', 'ESC', 'ES', 'fibroblast', 'liver', 'MCF', 'melanoma', 'ovary', 'PBMC', 'imr90', 'heart', 'ES', 'brain', 'immune', 'ovary', 'testes', 'testis', 'uterus', 'oocytes', 'ovarian']
	useTissue = "NONE"
	for tissue in tissueList:
		if tissue in dir: useTissue = tissue
	
	#check if not overwriting...
	if fileName not in fileDict:
		fileDict[fileName] = dir + '\t' + org + '\t' + adap + '\t' + useTissue

metaFile = open(metaFileName, 'w')
for file in fileDict:
	metaFile.write('%s\t%s\n' % (file, fileDict[file]))
metaFile.close()

#check how complete the meta file is:
print 'Completeness of Metafile:'
metaFile = open(metaFileName, 'r')
i = 0
j = 0
for line in metaFile:
	i = i + 1
	if line.strip().split('\t')[2] == 'NONE':
		j = j + 1
		print '  ', line.strip()
metaFile.close()

print 'incomplete (organism): (%s/%s)' % (j,i)

metaFile = open(metaFileName, 'r')
i = 0
j = 0
for line in metaFile:
	i = i + 1
	if line.strip().split('\t')[3] == 'NONE':
		j = j + 1
		print '  ', line.strip()
metaFile.close()

print 'incomplete (adapter): (%s/%s)' % (j,i)


	
