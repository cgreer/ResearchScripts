#given a coordinate, return wig value
import bioLibCG as cg
import cgConfig as c

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()
wigDir = mConf.conf['wigMouse']

def getWigValueLINE(chrom, strand, point):
	'''Old ONE -> Use "Byte" one'''
	
	#no strand specification in sequencing
	point = int(point)
	if int(strand) == 1:
		fN = wigDir + '/Merge.mouse.1.wig.%s.wig' % chrom
		fNindex = wigDir + '/Merge.mouse.1.wig.%s.wig.index' % chrom
	else:
		fN = wigDir + '/Merge.mouse.-1.wig.%s.wig' % chrom
		fNindex = wigDir + '/Merge.mouse.-1.wig.%s.wig.index' % chrom
	
	#get line in index file
	iFile = open(fNindex, 'r')
	startLine = 0
	for line in iFile:
		beg = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		
		if beg <= point <= end:
			startLine = int(cg.ss(line)[0])
			break
	iFile.close()
	
	#grab value
	f = open(fN, 'r')
	i = 0
	while i < startLine:
		f.readline() #skip header and lines till indexed line...
		i += 1
	wigValue = 0
	for line in f:
		beg = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		
		if beg <= point < end:
			wigValue += float(cg.ss(line)[3].split('.')[0])
			break
	f.close()

	return wigValue

def getWigValue(chrom, strand, point):
	'''Uses Byte Indexes'''
	
	#no strand specification in sequencing
	point = int(point)
	if int(strand) == 1:
		fN = wigDir + '/Merge.mouse.1.wig.%s.wig' % chrom
		fNindex = wigDir + '/Merge.mouse.1.wig.%s.wig.index' % chrom
	else:
		fN = wigDir + '/Merge.mouse.-1.wig.%s.wig' % chrom
		fNindex = wigDir + '/Merge.mouse.-1.wig.%s.wig.index' % chrom
	
	#get line in index file
	iFile = open(fNindex, 'r')
	startByte = 'None'
	for line in iFile:
		beg = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		
		if beg <= point < end:
			startByte = int(cg.ss(line)[0])
			#print 'INDEX', line.strip()
			break
	iFile.close()
	
	#grab value
	f = open(fN, 'r')
	f.seek(startByte, 0)
	'''
	s = ""
	i = 0
	while i < 20:
		s += f.read(1)
		i += 1
	print s
	'''
	wigValue = 0
	for line in f:
		#print 'Line:', line.strip()
		beg = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		
		if beg <= point < end:
			wigValue += float(cg.ss(line)[3].split('.')[0])
			break
	f.close()

	return wigValue

if __name__ == "__main__":
	import sys
	
	getWigValue(sys.argv[1], sys.argv[2], sys.argv[3])
