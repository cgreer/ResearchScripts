import cgPeaks
import compareData as compare
import stepVectorScan as svs
import cgPlot as plot
import os
import bioLibCG as cg

def profileTargetsHisto(tccList, cName, name = 'boxplot'):
		
	histDict = {} # {coord: []}
	for tcc in tccList:
		
		chrom, strand, start, end = cg.tccSplit(tcc)
		#Get highest peak
		tccStretch = cgPeaks.stretch(tcc, cName)
		tccStretch.createPeaks(span = 2)
		highestCoord = tccStretch.getHighestPeak()
		if highestCoord == None: continue
		
		#profile around point
		zPoint = cg.makeTcc(chrom, strand, highestCoord, end)
		cProfile = svs.profileAroundPoint(zPoint, 200, cName, ratio = True)
		
		for coord in cProfile:
			try:
				histDict[coord].append(cProfile[coord])
			except: #quicker way to initialize
				histDict[coord] = [cProfile[coord]]
	
	
	plot.boxPlotHisto(histDict, name = name)

def profileTargetsHistoAS(tccList, cName, name = 'boxplot'):
	
	range = 50
	histDict = {} # {coord: []}
	histDictAS = {}
	for tcc in tccList:
		
		chrom, strand, start, end = cg.tccSplit(tcc)
		#Get highest peak (sense)
		tccStretch = cgPeaks.stretch(tcc, cName)
		tccStretch.createPeaks(span = 2)
		highestCoord = tccStretch.getHighestPeak()
		if highestCoord == None: continue
		
		#AS
		tccAS = cg.convertToAS(tcc)
		tccStretch = cgPeaks.stretch(tccAS, cName)
		tccStretch.createPeaks(span = 2)
		highestCoordAS = tccStretch.getHighestPeak()
		if highestCoordAS == None: continue
		
		#profile around point (Sense)
		zPoint = cg.makeTcc(chrom, strand, highestCoord, end)
		cProfile = svs.profileAroundPoint(zPoint, range, cName, ratio = True)
		
		for coord in cProfile:
			try:
				histDict[coord].append(cProfile[coord])
			except: #quicker way to initialize
				histDict[coord] = [cProfile[coord]]
	
		#profile around point (AS)
		zPoint = cg.convertToAS(zPoint)
		cProfile = svs.profileAroundPoint(zPoint, range, cName, ratio = True, ratioCoord = highestCoordAS)
		
		for coord in cProfile:
			try:
				histDictAS[coord].append(cProfile[coord])
			except: #quicker way to initialize
				histDictAS[coord] = [cProfile[coord]]
	
	plot.boxPlotHistoAS(histDict, histDictAS, name = name)

def profileTargets(tccList, cName, dir = None, min = 0):
	if not dir:
		dir = 'out'
		
	cg.clearDirectory(dir) #clear and/or make
		
	for tcc in tccList:
		plot.plotASProfile(tcc, cName, directory = dir, min = min)
		
		
		
		
	
	
