import bioLibCG as cg
import cgConfig as c
import cgPeaks
import stepVectorScan
import subprocess
import os

def rangePoints(low, high, numPoints):
	'''for ints only'''
	points = []
	
	step = (high - low)/numPoints
	points.append(low)
	i = low
	while i < high:
		i += step
		points.append(i)
	#points.append(high -50)#!!!BECAUSE I NEED TO REMAKE TRACKS FOR MOUSE + HUMAN
	return points
	
	

def makePeakInput(cName, minExpression = 2000):
	
	mConf = c.getConfig('Main.conf')
	conf = c.getConfig(cName)
	
	assembly = conf.conf['assembly']
	
	tccList = []
	
	chromLens = cg.returnChromLengthDict(assembly)
	f = open('peakData.%s' % minExpression, 'w')
	for chrom in chromLens:
		if chrom not in cg.acceptableChroms: continue
		for strand in ['1', '-1']:
			print 'Getting Peaks for ', chrom, strand
			prevI = 0
			endCheck = 0
			for i in rangePoints(1, chromLens[chrom], 1000):
				if i == 1:
					prevI = i
					continue
				
				start = prevI
				end = i
				prevI = i
				
				tcc = cg.makeTcc(chrom, strand, start, end)
				#print 'scanning range', tcc
				peaks = cgPeaks.stretch(tcc, cName)
				peaks.createPeaks(span = 3, minVal = minExpression)
				
				for x in peaks.peaks:
					
					if x < endCheck:
						continue
				
					#scan a 30 bp range around this point and find the best roof...
					pRange = 30
					rTcc = cg.makeTcc(chrom, strand, x, x + 1)
					
	
					#now make profile for roof...
					cProfile = stepVectorScan.profileAroundPoint(rTcc, pRange, cName, ratio = True)
					
					
					
					#now get highest stretch length and the rNext coord.
					minVal = .80
					highest = 0
					stretch = 0
					startCurrent = None
					startFinal = None
					endFinal = None
					for i in range(1 - pRange, pRange):
						if cProfile[i] > minVal:
							stretch += 1
							if startCurrent == None:
								startCurrent = i
						else:
							if stretch > 0:
								if stretch > highest: #stretch ended and was higher than previous
									highest = stretch
									endFinal = i - 1
									startFinal = startCurrent
									startCurrent = None
								else:
									startCurrent = None
							stretch = 0
					
					#get +/- 4 value...
					val = [1.0, 1.0]
					if (startFinal) and (endFinal):
						low = startFinal - 4
						high = endFinal + 4
						if low > (1 - pRange) and high < pRange:
								val[0] = float(cProfile[startFinal - 4])
								val[1] = float(cProfile[endFinal + 4])
						else:
							continue
					else:
						continue
					
					endCheck = x + high
					
					#filter out peaks that look a certain way.
					if 14 < highest < 26: #rooflength
						if val[0] < 0.2 and val[1] < .2: #drop values
							goodTcc = cg.makeTcc(chrom, strand, x + low, x + high)
							#print goodTcc
							f.write('%s\n' % goodTcc)
	f.close()
	
def makePeakInputQ(cName, minExpression = 2000):
	'''Uses shell script and qsub to get peaks quickly'''
	
	mConf = c.getConfig('Main.conf')
	conf = c.getConfig(cName)
	
	assembly = conf.conf['assembly']
	
	tccList = []
	
	chromLens = cg.returnChromLengthDict(assembly)
	
	for chrom in chromLens:
		if chrom not in cg.acceptableChroms: continue
		for strand in ['1','-1']:
			print 'Getting Peaks for ', chrom, strand
			prevI = 0
			for i in rangePoints(1, chromLens[chrom], 30):
				if i == 1:
					prevI = i
					continue
				
				start = prevI
				end = i
				prevI = i
				
				tcc = cg.makeTcc(chrom, strand, start, end)
								
				log = 'logs/o-' + str(start)
				elog = 'logs/e-%s-%s-%s-%s' % (chrom, strand, start, end)
				subprocess.Popen(['qsub', '-V', '-cwd', '-e', elog, '-o', log, '-l', 'mem=3G', '-l', 'rt=3600', 'q.sh', tcc, cName, str(minExpression)]).wait()
				
					
def mergeInputs(cName, eLevel):
	
	conf = c.getConfig(cName)
	assembly = conf.conf['assembly']
	ending = '%s.%s' % (eLevel, assembly)
	
	print 'merging all files with ending', ending
	
	newLines = []
	for fN in cg.recurseDir('out', end = ending):
		print os.getcwd(), fN
		fN = os.getcwd() + '/' + fN
		f = open(fN, 'r')
		newLines.extend(f.readlines())
		f.close()
	
	f = open('peakData.%s.%s' % (eLevel, assembly), 'w')
	f.writelines(newLines)
	f.close()
	
							

if __name__ == "__main__":
	import sys
		
	makePeakInput(sys.argv[1])
