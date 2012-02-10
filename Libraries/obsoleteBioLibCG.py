def convertDcdToTcc(dcdList):
	tccList = []
	for dcd in dcdList:
		chr = dcd.strip().split(':')[0]
		strand = dcd.strip().split(':')[1]
		if strand == '+':#switch to numbers if they aren't already.
			strand = '1'
		elif strand == '-':
			strand = '-1'
		start = dcd.strip().split(':')[2].split('-')[0]
		end = dcd.strip().split(':')[2].split('-')[1]
		
		tccList.append('%s:%s:%s:%s' % (chr, strand, start, end))
		
	return tccList
def stripTripleColon(coord = None):
	
	if coord is None:
		print 'failed to pass coordinate'
		return False
		
	coordDict = {}
	coordDict['chromosome'] = coord.strip().split(':')[0]
	coordDict['strand'] = coord.strip().split(':')[1]
	coordDict['start'] = coord.strip().split(':')[2]
	coordDict['end'] = coord.strip().split(':')[3]
	
	return coordDict
