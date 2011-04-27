#given a coordinate, return wig value
import bioLibCG as cg
import cgConfig as c

#init
mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()
wigDir = mConf.conf['wigMouse']

def getWigValue(chrom, point, fN):
	
	#no strand specification in sequencing
	point = int(point)
	
	#get line in index file
	
	#grab value
	f = open(fN, 'r')
	f.readline()
	wigValue = 0
	for line in f:
		beg = int(cg.ss(line)[1])
		end = int(cg.ss(line)[2])
		fChrom = cg.ss(line)[0]
		
		if beg <= point < end:
			if chrom == fChrom:
				wigValue += float(cg.ss(line)[3].split('.')[0])
				break
	f.close()

	return wigValue

def getAll(chrom, strand, point):
	fNs = cg.recurseDir(mConf.conf['smallPath'], end = '.wig')
	for file in fNs:
		if 'WIG' in file:
			fNs.remove(file)
		elif file == '/home/chrisgre/smallLibs/WIGS/Merge.mouse.1.wig' or '/home/chrisgre/smallLibs/WIGS/Merge.human.1.wig':
			fNs.remove(file)
	
	for fN in fNs:
		fStrand = cg.ss(fN, '.')[-2]
		if str(fStrand) == str(strand):
			val = getWigValue(chrom, point, fN)
			if val > 0:
				print fN, val
		
		
	

if __name__ == "__main__":
	import sys
	
	getWigValue(sys.argv[1], sys.argv[2], sys.argv[3])
