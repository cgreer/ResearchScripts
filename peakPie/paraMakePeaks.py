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
			        pass
                        else:
                                startFinal = i + 1
                                #print 'found left', startFinal, cProfile[i]
                                break
                #right
                endFinal = rightRange[-1] #this only holds if it extends to the end of the range...
                for i in rightRange:
                        if cProfile[i] > minRoofVal:
                                pass
                        else:
                                endFinal = i - 1
                                #print 'found right', endFinal, cProfile[i]
                                break

	        peakLength = endFinal - startFinal + 1

	        extend = int(extend)	
                low = startFinal - extend
                high = endFinal + extend

                print peakLength
                if low > (1 - pRange) and high < pRange:
                        dropPassL = False
                        dropPassR = False

                        #find if any of the values in the extended range are below drop range
                        leftDrop = [float(cProfile[startFinal - x]) for x in range(1, extend + 1)]
                        rightDrop = [float(cProfile[endFinal + x]) for x in range(1, extend + 1)]
                        leftDropPass = [True if x < maxDropVal else False for x in leftDrop]
                        rightDropPass = [True if x < maxDropVal else False for x in rightDrop]
                       
                        if True in leftDropPass:
                                dropPassL = True
                        if True in rightDropPass:
                                dropPassR = True
                        
                        if (not dropPassL) or (not dropPassR):
                                #print 'dropVal Fail', 'dropLeft', 'dropRight'
                                #print startFinal, endFinal, leftDrop, rightDrop, leftDropPass, rightDropPass 
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
		#if (minPeakLength < peakLength < maxPeakLength):
                        goodTcc = cg.makeTcc(chrom, strand, peakPosition + startFinal, peakPosition + endFinal)
                        #print '*KEEPER', goodTcc, peakLength, avgNoise
                        return goodTcc
                else:
                        #print 'bad peak', chrom, strand, peakPosition + startFinal, peakPosition + endFinal
                        #print 'reason', peakLength, avgNoise
                        return False

def extendPeakTest(tcc, pRange, minVal, maxAvgNoise, minPeakLength, maxPeakLength, cName):  	
	        
                chrom, strand, peakPosition, end = cg.tccSplit(tcc)
		cProfile = stepVectorScan.profileAroundPoint(tcc, pRange, cName, ratio = True)
                cProfile2 = stepVectorScan.profileAroundPoint(tcc, pRange, cName, ratio = False) #just for debugging
 
		
                #extend this peak left and right
                leftRange = range(1-pRange, 0)
                rightRange = range(1, pRange)
                leftRange.reverse() #going from the middle outward

                #left
                startFinal = leftRange[-1]
		for i in leftRange:
			if cProfile[i] > minVal:
                                pass
                                #print ' extending stretch'
			else:
                                #print ' end of stretch L'
                                startFinal = i + 1
                                break
                #right
                endFinal = rightRange[-1]
                for i in rightRange:
                        if cProfile[i] > minVal:
                                pass
                                #print ' extending stretch'
                        else:
                                #print ' end of stretch R'
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
                print totalLength, pRange, low, high, lowRange, highRange, cProfile2[0]
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
	chrom, strand, start, end = cg.tccSplit(tcc)
        peaks = cgPeaks.stretch(tcc, cName)
		
	print 'getting peaks'
	peaks.createPeaks(span = 1, minVal = int(minExpression))
        	
	for x in peaks.peaks:

		print x
                
                newTcc = cg.makeTcc(chrom, strand, x, x + 1)
                #testedPeak = extendPeakTest(newTcc, 20, .2, .05, 0, 6, cName) 
                testedPeak = roofPeakTest(newTcc, 30, .85, .9, .2, 6, 17, 24, cName)

                if testedPeak:
                        f.write('%s\n' % testedPeak)
	

	f.close()

if __name__ == "__main__":
	import sys
	
	parallelMakePeaks(sys.argv[1], sys.argv[2], sys.argv[3])
