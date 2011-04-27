import bioLibCG
import cgDB
import cgAlignment
import cgOriginRNA
import dumpObj

def filterSmall(smallFN, inTranscript, inMicro, eLevelMin, outFN):
        eLevelMin = int(eLevelMin)
        
        f = open(smallFN, 'r')
        fOut = open(outFN, 'w')

        for line in f:
                ls = line.strip().split('\t')
                eLevel = ls[3]
                if eLevel == '.':
                        eLevel = 0
                else:
                        eLevel = int(ls[3])

                if inMicro:
                        microStat = ls[5]
                        if not microStat == '1':
                                continue

                if inTranscript:
                        tranStat = ls[6]
                        if not tranStat == '1':
                                continue

                if eLevel < eLevelMin:
                        continue

                #passed all filters --> write
                fOut.write(line)

def filterTargets(oRNADir, aDir, inTranscript, misLevel, centerLevel, minCenterLevel):
        if inTranscript == 'True': inTranscript = True
        if inTranscript == 'False': inTranscript = False
        misLevel, centerLevel, minCenterLevel =  int(misLevel), int(centerLevel), float(minCenterLevel)

        oDC = cgDB.dataController(oRNADir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()

        for oRNA in id_oRNA.values():
                oRNA.filteredTargets = []
                for aID in oRNA.targets:
                        alignment = id_alignment[aID]

                        #transcriptOverlap
                        if inTranscript:
                                if not alignment.transcriptOverlap:
                                        #print 'tOverlap Fail', cgAlignment.pretty#print(alignment)
                                        continue

                        #misLevel
                        if alignment.mismatchStatus[misLevel]:
                                #print 'mismatch Fail', cgAlignment.pretty#print(alignment)
                                continue
                           
                        #centerLevel
                        if alignment.centerExpression[centerLevel] < minCenterLevel:
                                #print 'expression Fail', cgAlignment.pretty#print(alignment)
                                continue
                        
                        oRNA.filteredTargets.append(aID)

        oDC.commit(id_oRNA)                        


def filterOutTargets(resultsFN, centerFN, mismatchFN, targetFN, tranCheck, mPick, cPick, minCenterLevel, inputPosition, updatePosition, outFN):
        '''Pick is the range in which you need the values for. 0 is 4bp around, 1 is 6...'''	

        #make mismatch dict
        mmDict = {} # siD: tID : mmVal
        f = open(mismatchFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                sID = int(ls[0])
                tID = int(ls[1])
                mmVal = int(ls[2 + mPick])
                if not sID in mmDict:
                        mmDict[sID] = {}

                mmDict[sID][tID] = mmVal
        f.close()

        #make center dict
        centerDict = {} # sID: tID: centerVal
        f = open(centerFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                sID = int(ls[0])
                tID = int(ls[1])
                centerVal = float(ls[2 + cPick])
                if not sID in centerDict:
                        centerDict[sID] = {}

                centerDict[sID][tID] = centerVal
        f.close()

        #make transcript target dict
        tranVals = {} # tID : tranValue
        f = open(targetFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                tID = int(ls[0])
                tranVal = ls[3]
                tranVals[tID] = tranVal

        #go through and test all the targets.

	f = open(resultsFN, 'r')
	newLines = []
	for line in f:
                sID = int(line.strip().split('\t')[0])
		targets = line.strip().split('\t')[int(inputPosition)]
		targets = targets.strip().split(',')
                
                newTargetList = []
                for tID in targets:
                        ##print sID, tID
                        tID = int(tID)
                        
                        #Check for inside transcript
                        if tranCheck:
                                if tranVals[tID] == '0':
                                        continue

                        #check mismatches
                        mmVal = mmDict[sID][tID]
                        if mmVal == 1:
                                continue

                        #check center Expression
                        centerVal = centerDict[sID][tID]
                        if centerVal < minCenterLevel:
                                continue
	        	
                        newTargetList.append(str(tID))

                if len(newTargetList) < 1: continue 
                newTargets = ','.join(newTargetList)

		#update newLines
	        newLines.append(bioLibCG.appendToLine(line, newTargets, int(updatePosition)))
                
	f.close()
	
	
	#update file
	f = open(outFN, 'w')
	f.writelines(newLines)
	f.close()



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterTargets, sys.argv)
