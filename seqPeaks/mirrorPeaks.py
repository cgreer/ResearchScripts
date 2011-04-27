import cgPeaks
import compareData as compare
import math
import bioLibCG as cg

knowns = compare.tccFileToList('mouseKnownMirs.tcc', 0)

eLevels = []
for known in knowns:
	
	chrom, strand, start, end = cg.tccSplit(known, True) #text...
	if strand == '1':
		strand = '-1'
	else:
		strand = '1'
	oppTcc = cg.makeTcc(chrom, strand, start, end)
	
	knownStretch = cgPeaks.stretch(known)
	knownStretch.createPeaks(1,20)
	kPos = knownStretch.getHighestPeak()
	if kPos: eLevels.append(knownStretch.profile[kPos])
	
	oppStretch = cgPeaks.stretch(oppTcc)
	oppStretch.createPeaks(1,20)
	oPos = oppStretch.getHighestPeak()
	
	if oPos and kPos:
		#determine if they are close enough to be considered mirrored...
		if math.fabs(int(kPos) - int(oPos)) < 12:
			print known, oPos, kPos, oppStretch.profile[oPos], knownStretch.profile[kPos]


print eLevels
