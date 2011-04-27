#get the hairpin values
#can also reuse elements from this script to get CID elements...
import bioLibCG as cg
from bioLibCG import ss
import cgConfig as c

mConf = c.cgConfig('Main.conf')
conf = c.cgConfig()

def getHairpins(fN = None):
	
	predFile = open(fN, 'r')

	#populate CID:hairpin range
	cHairs = {}
	for line in predFile:
		#get cluster ID
		CID = ss(line)[7]
		hairpin = ss(line)[2]
		
		if CID in cHairs:
			#check if the starts and ends need to be stretched
			hStart = int(ss(cHairs[CID], ':')[2])
			hEnd = int(ss(cHairs[CID], ':')[3])
			
			start = int(ss(hairpin, ':')[2])
			end = int(ss(hairpin, ':')[3])
			
			if start < hStart:
				hStart = start
			if end > hEnd:
				hEnd = end
			
			cHairs[CID] = '%s:%s:%s:%s' % (ss(hairpin, 1)[0], ss(hairpin, 1)[1], hStart, hEnd)
		else:
			cHairs[CID] = hairpin
	predFile.close()
	
	return cHairs
	
if __name__ == "__main__":
	import sys
	
	getHairpins()

