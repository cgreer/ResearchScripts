import bioLibCG
import cgPeaks
import cgDB
import cgAlignment

def markCenterExpression(aDir, cName):
        
        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()


        for alignment in id_alignment.values():
                alignment.centerExpression = [0.0, 0.0, 0.0]      
                chrom, strand, start, end = bioLibCG.tccSplit(alignment.tTcc)
                offset = alignment.tStart
                sLen = alignment.sLength

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
                        alignment.centerExpression[0] = sum/expressionSum

                        sum = 0.0
                        for key in sortedKeys[7:13]:
                                sum += stretch.profile[key]
                        alignment.centerExpression[1] = sum/expressionSum

                        sum = 0.0
                        for key in sortedKeys[6:14]:
                                sum += stretch.profile[key]
                        alignment.centerExpression[2] = sum/expressionSum

        aDC.commit(id_alignment)

def markCenterExpressionOLD(smallFN, targetFN, alignmentFN, cName, outFN):
        
        #print 'making target dict'
        #make targetDict
        f = open(targetFN, 'r')
        targetDict = {} # tID: tLoc
        for line in f:
                ls = line.strip().split('\t')
                targetDict[int(ls[0])] = ls[1]
        f.close()

        #print 'making alignment dict'
        #make alignmentDict
        alignDict = {} # sid: {target: offset}
        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                sID = int(ls[0])
                tID = int(ls[1])
                offset = int(ls[4])
                if not sID in alignDict:
                           alignDict[sID] = {}

                alignDict[sID][tID] = offset #assumes one source to target...
        f.close()

        f = open(smallFN, 'r')
        fOut = open(outFN, 'w')

        for line in f:
                ls = line.strip().split('\t')
                sID = int(ls[0])
                sLoc = ls[1]
                sLen = len(sLoc) #This is the sequence for simulated reads... 
                #sLen = bioLibCG.getTccLength(sLoc) #off by one?
                tIDs = ls[4].split(',')

                for tID in tIDs:
                        tID = int(tID)
                        tLoc = targetDict[tID]
                        chrom, strand, start, end = bioLibCG.tccSplit(tLoc)
                        offset = alignDict[sID][tID]

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
                        

                        lowE = 0.0
                        midE = 0.0
                        highE = 0.0
                        

                        if expressionSum != 0:

                                sum = 0.0
                                for key in sortedKeys[8:12]:
                                        sum += stretch.profile[key]
                                lowE = sum/expressionSum

                                sum = 0.0
                                for key in sortedKeys[7:13]:
                                        sum += stretch.profile[key]
                                midE = sum/expressionSum

                                sum = 0.0
                                for key in sortedKeys[6:14]:
                                        sum += stretch.profile[key]
                                highE = sum/expressionSum

                        
                        fOut.write('%s\t%s\t%s\t%s\t%s\n' % (sID, tID, lowE, midE, highE))


                        
def markMismatchedPairs(aDir):
      
        #make mismatchDict
        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()
        
        for alignment in id_alignment.values():

                alignment.mismatchStatus = [False, False, False]
                lowRange = range(9,13)
                midRange = range(8,14)
                highRange = range(7,15)
                
                #check mismatches
                for i in lowRange:
                        if i in alignment.mismatchPositions:
                                alignment.mismatchStatus[0] = True
                                break

        
                for i in midRange:
                        if i in alignment.mismatchPositions:
                                alignment.mismatchStatus[1] = True
                                break

                for i in highRange:
                        if i in alignment.mismatchPositions:
                                alignment.mismatchStatus[2] = True
                                break

        aDC.commit(id_alignment)

def markMismatchedPairsOLD(smallResultsFN, alignmentFN, oFN):
        
        #make mismatchDict
        alignDict = {} # sid: {target: [mismatches]}
        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                sID = int(ls[0])
                tID = int(ls[1])
                if len(ls) > 9:
                        mismatches = [int(x) for x in ls[9].split(',')]
                else:
                        mismatches = []
                if not sID in alignDict:
                        alignDict[sID] = {}

                alignDict[sID][tID] = mismatches
        f.close()



        f = open(smallResultsFN, 'r')
        fOut = open(oFN, 'w')

        for line in f:
                #print line
                ls = line.strip().split('\t')
                sID = int(ls[0])
                tIDs = ls[4].split(',')
                
                #for every target, get middle mismatch position...
                for tID in tIDs:
                        #print sID, tID
                        tID = int(tID)
                        alignment.mismatchPositions = alignDict[sID][tID]
                        lowRange = range(9,13)
                        midRange = range(8,14)
                        highRange = range(7,15)
                        lowCheck = 0
                        midCheck = 0
                        highCheck = 0


                        #check mismatches
                        for i in lowRange:
                                if i in alignment.mismatchPositions:
                                        lowCheck = 1
                                        break

                
                        for i in midRange:
                                if i in alignment.mismatchPositions:
                                        midCheck = 1
                                        break

                        for i in highRange:
                                if i in alignment.mismatchPositions:
                                        highCheck = 1
                                        break

                        fOut.write('%s\t%s\t%s\t%s\t%s\n' % (sID, tID, lowCheck, midCheck, highCheck))
        f.close()


if __name__ == "__main__":
        import sys
        #bioLibCG.submitArgs(markMismatchedPairs, sys.argv)
        bioLibCG.submitArgs(markCenterExpression, sys.argv)
        
