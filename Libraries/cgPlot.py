import bioLibCG as cg
import stepVectorScan as svs
import cgPeaks
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import cgConfig as c
import os
import rpy2.rlike.container as rlc

def plotHistogram(valuesList, directory = None, name = 'histo'):
	
	if not directory:
		fN = name + '.png'
	else:
		fN = directory + '/' + name + '.png'
	
	
	#Convert to R acceptable form.
	values = tuple(valuesList)
	vals = robjects.IntVector(values)
	
	gDevice = importr('grDevices')
	gDevice.png(file=fN, width=1680, height=1050)
	#r.hist(vals, breaks=100, xlab = "X", ylab = "Frequency", main=name)
	r.hist(vals, breaks=100, xlab = "# Targets", ylab = "# MicroRNAs", main=name)
	gDevice.dev_off()
	
def plotProfile(tcc, cName, directory = None, min = 0):
	if not directory:
		fN = tcc + '.png'
	else:
		fN = directory + '/' + tcc + '.png'
		
	tccStretch = cgPeaks.stretch(tcc, cName)
	highest = tccStretch.getHighestLevel()
	if highest < min:
		return 0
		
	sortedX = tccStretch.profile.keys()
	sortedX.sort()
	
	sortedY = []
	for X in sortedX:
		sortedY.append(tccStretch.profile[X])
	
	gDevice = importr('grDevices')
	gDevice.png(file=fN, width=1680, height=1050)
	r.plot(sortedX, sortedY, xlab = "Coordinates", ylab = "Expression Level")
	r.lines(sortedX, sortedY, type = "b")
	gDevice.dev_off()

def plotSmallDeg(tcc, smallCName, degCName, outDir = None, description = "None", nameNum = "0"):
	
        if not outDir:
		fN = nameNum + "." + tcc + '.png'
	else:
		fN = outDir + '/' + nameNum + "." + tcc + '.png'
	
        #Get deg Profile
	tccStretch = cgPeaks.stretch(tcc, degCName)
		
	sortedX = tccStretch.profile.keys()                                                                                                     
	sortedX.sort()
	
	sortedY = []
	for X in sortedX:
		sortedY.append(tccStretch.profile[X])
	
	#Get small
	tccStretchSmall = cgPeaks.stretch(tcc, smallCName)
		
	sortedXAS = tccStretchSmall.profile.keys()
	sortedXAS.sort()
	
	sortedYAS = []
	for X in sortedXAS:
		sortedYAS.append(tccStretchSmall.profile[X])
	
	#Plot them
	gDevice = importr('grDevices')
	gDevice.png(file=fN, width=1680, height=1050)
	r('split.screen(c(2,1))')
	r('screen(1)')
	r.plot(sortedX, sortedY, xlab = "Coordinates", ylab = "Degradome Expression" )
	r.lines(sortedX, sortedY, type = "b")
	r('screen(2)')
	r.plot(sortedXAS, sortedYAS, xlab = description, ylab = "Small Expression")
	r.lines(sortedXAS, sortedYAS, type = "b")
	gDevice.dev_off()

def plotASProfile(tcc, cName, directory = None, min = 0, extra = "0"):
	if not directory:
		fN = extra + '.' + tcc + '.png'
	else:
		fN = directory + '/' + extra + '.' + tcc + '.png'
	
	#Get S Profile
	tccStretch = cgPeaks.stretch(tcc, cName)
	highest = tccStretch.getHighestLevel()
	if highest < min:
		return 0
		
	sortedX = tccStretch.profile.keys()
	sortedX.sort()
	
	sortedY = []
	for X in sortedX:
		sortedY.append(tccStretch.profile[X])
	
	#Get AS Profile
	chr, strand, start, end = tcc.strip().split(':')
	if strand == '1':
		strand = '-1'
	else:
		strand = '1'
	tcc = cg.makeTcc(chr, strand, start, end)
	
	tccStretchAS = cgPeaks.stretch(tcc, cName)
	highest = tccStretchAS.getHighestLevel()
	if highest < min:
		return 0 #AS can have minimum I guess...
		
	sortedXAS = tccStretchAS.profile.keys()
	sortedXAS.sort()
	
	sortedYAS = []
	for X in sortedXAS:
		sortedYAS.append(tccStretchAS.profile[X])
	
	#Plot them
	gDevice = importr('grDevices')
	gDevice.png(file=fN, width=1680, height=1050)
	r('split.screen(c(2,1))')
	r('screen(1)')
	r.plot(sortedX, sortedY, xlab = "Coordinates", ylab = "(Syn) Expression Level" )
	r.lines(sortedX, sortedY, type = "b")
	r('screen(2)')
	r.plot(sortedXAS, sortedYAS, xlab = "Coordinates", ylab = "(Anti) Expression Level")
	r.lines(sortedXAS, sortedYAS, type = "b")
	gDevice.dev_off()

def boxPlotHisto(histoDict, name = 'boxplot'):
	
	gDevice = importr('grDevices')
	gDevice.png(file='%s.png' % name, width=1680, height=1050)
	
	# plotting code here
		
	sortedCoords = histoDict.keys()
	sortedCoords.sort()
        
        #Add NAs to rows that don't have same length as max row
        maxLen = 0
        for key in sortedCoords:
                l = len(histoDict[key])
                if l > maxLen: maxLen = l

        for key in sortedCoords:
                l = len(histoDict[key])
                while l < maxLen:
                        histoDict[key].append('NA')
                        l = len(histoDict[key])

	od = rlc.OrdDict() #an ordered dictionary, keeps track of how you add stuff to it.
	for coord in sortedCoords:
		values = tuple(histoDict[coord])
		coord = str(coord)
		od[coord] = robjects.FloatVector(values)
        
        print od

	#lim = robjects.IntVector((0,1))
	df = robjects.DataFrame(od)
	#r.boxplot(df, ylim=lim, outline = False)
	r.boxplot(df, outline = True)
	
	gDevice.dev_off()

def boxPlotHistoAS(histoDict, histoDictAS, name = 'boxplot'):
	
	gDevice = importr('grDevices')
	gDevice.png(file='%s.png' % name, width=1680, height=1050)
	r('split.screen(c(2,1))')
	
	#sense
	sortedCoords = histoDict.keys()
	sortedCoords.sort()
	
	od = rlc.OrdDict() #an ordered dictionary, keeps track of how you add stuff to it.
	for coord in sortedCoords:
		values = tuple(histoDict[coord])
		coord = str(coord)
		od[coord] = robjects.FloatVector(values)
	
	lim = robjects.IntVector((0,1))
	df = robjects.DataFrame(od)
	r('screen(1)')
	r.boxplot(df, ylim=lim, outline = False)
	
	#AS
	sortedCoords = histoDictAS.keys()
	sortedCoords.sort()
	
	od = rlc.OrdDict() #an ordered dictionary, keeps track of how you add stuff to it.
	for coord in sortedCoords:
		values = tuple(histoDictAS[coord])
		coord = str(coord)
		od[coord] = robjects.FloatVector(values)
	
	lim = robjects.IntVector((0,1))
	df = robjects.DataFrame(od)
	r('screen(2)')
	r.boxplot(df, ylim=lim, outline = False)
	
	gDevice.dev_off()


def plotPie(pieDict, title = "My Super Cool Pie Chart"):
	gDevice = importr('grDevices')
	gDevice.png(file='%s.png' % name, width=1680, height=1050)
	# plotting code here
	
	pieDict = {'cool': 200, 'neato': 150}
	
	df = robjects.DataFrame(pieDict)
	r.pie(df, main = title)
	
	
	gDevice.dev_off()
	
	
if __name__ == "__main__":
	import sys
	
	#take care of extra arg...
	sys.argv.append(None) #this way if only one the others are None...
	plotProfile(sys.argv[1], sys.argv[2])
	
