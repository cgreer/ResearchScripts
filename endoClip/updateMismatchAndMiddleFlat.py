import bioLibCG
import cgNexusFlat
import cgAlignmentFlat
import cgWig

def markCenterExpression(aFN, wigDir, rn = None, tn = None):

        extend = 25
        
        timer = bioLibCG.cgTimer()
        timer.start()

        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['centerExpression', 'tTcc', 'tStart', 'sLength'], [rn, tn])
        
        #load expression of degradome
        wigDict = cgWig.loadWig(wigDir)
        
        for aID in aNX.centerExpression:
                aNX.centerExpression[aID] = [0.0, 0.0, 0.0]      
                chrom, strand, start, end = bioLibCG.tccSplit(aNX.tTcc[aID])
                offset = aNX.tStart[aID]
                sLen = aNX.sLength[aID]

                if strand == '1':
                        start = start - extend + offset
                        end = start + sLen
                else:
                        end = end + extend - offset
                        start = end - sLen

                scanRange = bioLibCG.makeTcc(chrom, strand, start, end)
                stretch = cgWig.getExpressionProfile(scanRange, wigDict)
                expressionSum = sum(stretch.values())
                sortedKeys = stretch.keys()
                sortedKeys.sort()

                if strand == '-1':
                        sortedKeys.reverse()
                

                if expressionSum != 0:

                        sumE = 0.0
                        for key in sortedKeys[8:12]:
                                sumE += stretch[key]
                        aNX.centerExpression[aID][0] = sumE/expressionSum

                        sumE = 0.0
                        for key in sortedKeys[7:13]:
                                sumE += stretch[key]
                        aNX.centerExpression[aID][1] = sumE/expressionSum

                        sumE = 0.0
                        for key in sortedKeys[6:14]:
                                sumE += stretch[key]
                        aNX.centerExpression[aID][2] = sumE/expressionSum
        
        aNX.save()


                        
def markMismatchedPairs(aFN, rn = None, tn = None):
    
        #make mismatchDict
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['mismatchStatus', 'mismatchPositions'], [rn, tn])
        
        lowRange = range(8,12) # remember the small locations are 0-based, so 10 is 9
        midRange = range(7,13)
        highRange = range(6,14)
        for aID in aNX.mismatchStatus:

                aNX.mismatchStatus[aID] = [False, False, False]
                #check mismatches
                for i in lowRange:
                        if i in aNX.mismatchPositions[aID]:
                                aNX.mismatchStatus[aID][0] = True
                                break

        
                for i in midRange:
                        if i in aNX.mismatchPositions[aID]:
                                aNX.mismatchStatus[aID][1] = True
                                break

                for i in highRange:
                        if i in aNX.mismatchPositions[aID]:
                                aNX.mismatchStatus[aID][2] = True
                                break
                
        aNX.save()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
        
