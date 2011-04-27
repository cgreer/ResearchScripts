import compareData as compare

tccList = compare.tccFileToList('snos.tcc', 0)

collapsed = compare.collapseOverlaps(tccList)

for tcc in collapsed:
	print tcc
	
