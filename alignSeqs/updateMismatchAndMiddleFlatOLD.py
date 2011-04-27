import bioLibCG
import cgNexusFlat
import cgPeaks
import cgAlignmentFlat

def markCenterExpression(aFN, cName, rn = None, tn = None):
        
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['centerExpression', 'tTcc', 'tStart', 'sLength'], [rn, tn])

        for aID in aNX.centerExpression:
                aNX.centerExpression[aID] = [0.0, 0.0, 0.0]      
                chrom, strand, start, end = bioLibCG.tccSplit(aNX.tTcc[aID])
                offset = aNX.tStart[aID]
                sLen = aNX.sLength[aID]

                if strand == '1':
                        start = start - 19 + offset
                        end = start + sLen
                else:
                        end = end + 19 - offset
                        start = end - sLen

                scanRange = bioLibCG.makeTcc(chrom, strand, start, end)
                
                stretch = cgPeaks.stretch(scanRange, cName)
                expressionSum = stretch.getSumOfLevels()
                sortedKeys = stretch.profile.keys()
                sortedKeys.sort()

                if strand == '-1':
                        sortedKeys.reverse()
                
                if expressionSum != 0:

                        sum = 0.0
                        for key in sortedKeys[8:12]:
                                sum += stretch.profile[key]
                        aNX.centerExpression[aID][0] = sum/expressionSum

                        sum = 0.0
                        for key in sortedKeys[7:13]:
                                sum += stretch.profile[key]
                        aNX.centerExpression[aID][1] = sum/expressionSum

                        sum = 0.0
                        for key in sortedKeys[6:14]:
                                sum += stretch.profile[key]
                        aNX.centerExpression[aID][2] = sum/expressionSum

        aNX.save()


                        
def markMismatchedPairs(aFN, rn = None, tn = None):
    
        #make mismatchDict
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['mismatchStatus', 'mismatchPositions'], [rn, tn])
        
        for aID in aNX.mismatchStatus:

                aNX.mismatchStatus[aID] = [False, False, False]
                lowRange = range(9,13)
                midRange = range(8,14)
                highRange = range(7,15)
                
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
                        if i in aNX.mismatchPositions:
                                aNX.mismatchStatus[aID][2] = True
                                break

        aNX.save()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
        
