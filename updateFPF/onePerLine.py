file = open('/home/chrisgre/projects/NoncodingMouse/results/NCmouse-s3k8b17.results.sorted', 'r')

CIDlist = []
for line in file:
	CID = line.strip().split('\t')[7]
	if CID not in CIDlist:
		#make colon dash:
		tcc = line.strip().split('\t')[2]
		dblCol = '%s:%s-%s' % (tcc.split(':')[0], tcc.split(':')[2], tcc.split(':')[3])
		print line.strip() + '\n' + dblCol
		CIDlist.append(CID)
