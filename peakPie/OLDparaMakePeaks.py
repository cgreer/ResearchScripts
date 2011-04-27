import bioLibCG as cg
import cgConfig as c
import cgPeaks
import stepVectorScan
import subprocess
	
def parallelMakePeaks(tcc, cName, minExpression):
	conf = c.getConfig(cName)
	f = open('out/peakData.%s.%s.%s' % (tcc, minExpression, conf.conf['assembly']), 'w')
	print 'scanning range', tcc
	chrom, strand, start, end = cg.tccSplit(tcc)
	peaks = cgPeaks.stretch(tcc, cName)
	
	
	
	#print 'getting peaks'
	peaks.createPeaks(span = 1, minVal = int(minExpression))
	
	print 'len peaks', len(peaks.peaks)
	endCheck = 0
	for x in peaks.peaks:
		print x, endCheck
                
                '''
		if x < endCheck:
                        print 'endChecked'
			continue
	        '''

		#scan a 30 bp range around this point and find the best roof...
		pRange = 40
		rTcc = cg.makeTcc(chrom, strand, x, x + 1)
		

		#now make profile for roof...
		cProfile = stepVectorScan.profileAroundPoint(rTcc, pRange, cName, ratio = True)
		
		#now get highest stretch length and the rNext coord.
		minVal = .70
		highest = 0
		stretch = 0
		startCurrent = None
		startFinal = None
		endFinal = None
		for i in range(1 - pRange, pRange):
                        print ' ', x + i, cProfile[i] 
			if cProfile[i] > minVal:
				print '  extending stretch'
                                stretch += 1
				if startCurrent == None:
					startCurrent = i
			else:
				if stretch > 0:
					print 'end of stretch'
                                        if stretch > highest: #stretch ended and was higher than previous
						highest = stretch
						endFinal = i - 1
						startFinal = startCurrent
						startCurrent = None
					else:
						startCurrent = None
				stretch = 0
		
		#get +/- extend value...
		val = [1.0, 1.0]
                extend = 1
		if (startFinal) and (endFinal):
			low = startFinal - extend
			high = endFinal + extend
			if low > (1 - pRange) and high < pRange:
					val[0] = float(cProfile[startFinal - extend])
					val[1] = float(cProfile[endFinal + extend])
			else:
                                print 'out of range'
				continue
		else:
                        print 'no start and end of peak'
			continue
	        print low, high, x, endFinal
		endCheck = x + endFinal
		
                #avg expression around peak check...
                #get total expression before peak
                noiseExpression = 0
                lowRange = range(1 - pRange, low)
                highRange = range(high + 1, pRange) 
                totalLength = len(lowRange) + len(highRange)
                for i in lowRange:
                        noiseExpression += cProfile[i]
                for i in highRange:
                        noiseExpression += cProfile[i]
                avgNoise = noiseExpression/float(totalLength)


		#filter out peaks that look a certain way.
                print highest, val[0], val[1], avgNoise
		if 0 < highest < 5: #rooflength 14/26
			if val[0] < 0.20 and val[1] < .20: #drop values
                                if avgNoise < .3:
                                        goodTcc = cg.makeTcc(chrom, strand, x + low, x + high)
				        print '*KEEPER'
				        f.write('%s\n' % goodTcc)
	

	f.close()
	print 'DONE', tcc

if __name__ == "__main__":
	import sys
	
	parallelMakePeaks(sys.argv[1], sys.argv[2], sys.argv[3])
