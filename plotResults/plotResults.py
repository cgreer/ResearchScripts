import bioLibCG as cg
import stepVectorScan as svs
import cgPeaks
import cgConfig as c
import getHairpins
import os
import cgPlot

def plotResults(fN, cName = None):
		
	cHairs = getHairpins.getHairpins(fN) #CID: HAIRPIN
	
	directory = cg.getBaseFileName(fN)
	cg.clearDirectory(directory)
	
	#change the directory before plotting
	cwd = os.getcwd()
	os.chdir(directory)
	
	for CID in cHairs:
		print 'plotting:', CID
		cgPlot.plotASProfile(cHairs[CID], cName)
	
	os.chdir(cwd)
		
	
if __name__ == "__main__":
	import sys
	
	#take care of extra arg...
	sys.argv.append(None) #this way if only one the others are None...
	plotResults(sys.argv[1], sys.argv[2])
	
