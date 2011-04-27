

print 'Populating CID/Exon distribution data'
f = open('noiseLevelsExons.data', 'r')
exonDists = {} #cid: [exon dist]
header = f.readline()
order = {} # num:CID
for i,CID in enumerate(header.strip().split('\t')):
	order[i] = CID
	exonDists[CID] = []
	
for line in f:
	data = line.strip().split('\t')
	for i, dataPoint in enumerate(data):
		if dataPoint == 'NA' or dataPoint == '':
			continue
		else:
			dataPoint = float(dataPoint)
			CID = order[i]
			exonDists[CID].append(dataPoint)

print exonDists['1298']
