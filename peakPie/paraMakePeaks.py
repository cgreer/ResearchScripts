import bioLibCG as cg
import cgConfig as c
import cgPeaks
import stepVectorScan
import subprocess

def roofPeakTest(tcc, pRange, minRoofVal, maxAvgNoise, maxDropVal, extend, minPeakLength, maxPeakLength, cName):  	
	        '''Note: extend does not extend the coordinates into the final tcc, just 
                used for declaring peak'''

                chrom, strand, peakPosition, end = cg.tccSplit(tcc)
		cProfile = stepVectorScan.profileAroundPoint(tcc, pRange, cName, ratio = True)
		
                #extend this peak left and right
                leftRange = range(1-pRange, 0)
                rightRange = range(1, pRange)
                leftRange.reverse() #going from the middle outward

                #left
                startFinal = leftRange[-1]
		for i in leftRange:
			if cProfile[i] > minRoofVal:
				print ' extending stretch'
			else:
                                print ' end of stretch L'
                                startFinal = i + 1
                                break
                #right
                endFinal = rightRange[-1] #this only holds if it extends to the end of the range...
                for i in rightRange:
                        if cProfile[i] > minRoofVal:
                                print ' extending stretch'
                        else:
                                print ' end of stretch R'
                                endFinal = i - 1
                                break

	        peakLength = endFinal - startFinal + 1

	        extend = int(extend)	
		val = [0.0, 0.0]
                low = startFinal - extend
                high = endFinal + extend
                
                if low > (1 - pRange) and high < pRange:
                                val[0] = float(cProfile[startFinal - extend])
                                val[1] = float(cProfile[endFinal + extend])
                                if not (val[0] < maxDropVal and val[1] < maxDropVal):
                                        return False
                                        
                else:
                        print 'out of range'
                        return False
                
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
		if (minPeakLength < peakLength < maxPeakLength) and (avgNoise < maxAvgNoise):
                        goodTcc = cg.makeTcc(chrom, strand, peakPosition + startFinal, peakPosition + endFinal)
                        print '*KEEPER'
                        return goodTcc
                else:
                        return False

def extendPeakTest(tcc, pRange, minVal, maxAvgNoise, minPeakLength, maxPeakLength, cName):  	
	        
                chrom, strand, peakPosition, end = cg.tccSplit(tcc)
		cProfile = stepVectorScan.profileAroundPoint(tcc, pRange, cName, ratio = True)
		
                #extend this peak left and right
                leftRange = range(1-pRange, 0)
                rightRange = range(1, pRange)
                leftRange.reverse() #going from the middle outward

                #left
                startFinal = leftRange[-1]
		for i in leftRange:
			if cProfile[i] > minVal:
				print ' extending stretch'
			else:
                                print ' end of stretch L'
                                startFinal = i + 1
                                break
                #right
                endFinal = rightRange[-1]
                for i in rightRange:
                        if cProfile[i] > minVal:
                                print ' extending stretch'
                        else:
                                print ' end of stretch R'
                                endFinal = i - 1
                                break

	        peakLength = endFinal - startFinal + 1

		
                #avg expression around peak check...
                #get total expression before peak
                low = startFinal
                high = endFinal
                noiseExpression = 0
                lowRange = range(1 - pRange, low)
                highRange = range(high + 1, pRange)
                totalLength = len(lowRange) + len(highRange)
                print totalLength, pRange, low, high, lowRange, highRange
                for i in lowRange:
                        noiseExpression += cProfile[i]
                for i in highRange:
                        noiseExpression += cProfile[i]
                try:
                        avgNoise = noiseExpression/float(totalLength)
                except:
                        return False

		#filter out peaks that look a certain way.
		if (minPeakLength < peakLength < maxPeakLength) and (avgNoise < maxAvgNoise):
                        goodTcc = cg.makeTcc(chrom, strand, peakPosition + startFinal, peakPosition + endFinal)
                        print '*KEEPER'
                        return goodTcc
                else:
                        return False

def parallelMakePeaks(tcc, cName, minExpression):
	conf = c.getConfig(cName)
	f = open('out/peakData.%s.%s.%s' % (tcc, minExpression, conf.conf['assembly']), 'w')
	print 'scanning range', tcc
	chrom, strand, start, end = cg.tccSplit(tcc)
	peaks = cgPeaks.stretch(tcc, cName)
		
	#print 'getting peaks'
	peaks.createPeaks(span = 1, minVal = int(minExpression))
	
	print 'len peaks', len(peaks.peaks)
	for x in peaks.peaks:
		print x
                
                newTcc = cg.makeTcc(chrom, strand, x, x + 1)
                #testedPeak = extendPeakTest(newTcc, 20, .2, .05, 0, 6, cName) 
                testedPeak = roofPeakTest(newTcc, 30, .85, .9, .2, 2, 17, 24, cName)

                if testedPeak:
                        f.write('%s\n' % testedPeak)
	

	f.close()
	print 'DONE', tcc

if __name__ == "__main__":
	import sys
	
	parallelMakePeaks(sys.argv[1], sys.argv[2], sys.argv[3])
