import bioLibCG
import cgAlignmentFlat
import cgOriginRNAFlat
import dumpObj
import cgNexusFlat


def filterOrigin(oFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'endContigLength', 'totalContigLength', 'sequenceDuplicate', 'passedFilter', 'entropy'], [rn, tn])

        for oID in oNX.entropy:
                
                oNX.passedFilter[oID] = False
                
                if oNX.entropy[oID] < 1.15:
                        continue
                
                if oNX.endContigLength[oID] > 6:
                        continue

                if oNX.totalContigLength[oID] > 6:
                        continue

                if oNX.sequenceDuplicate[oID]:
                        continue

                if not oNX.filteredTargets[oID]:
                        continue

                oNX.passedFilter[oID] = True

        oNX.save()                
                

def filterTargets(oFN, aFN, inTranscript, misLevel, centerLevel, minCenterLevel, rn = None, tn = None):
        if inTranscript == 'True': inTranscript = True
        if inTranscript == 'False': inTranscript = False
        misLevel, centerLevel, minCenterLevel =  int(misLevel), int(centerLevel), float(minCenterLevel)

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['filteredTargets', 'targets'], [rn, tn])
        
        print inTranscript, misLevel, centerLevel, minCenterLevel, oFN, aFN, rn, tn

        #make selection set
        targets = set()
        for oID in oNX.targets:
                for target in oNX.targets[oID]:
                        targets.add(target)

        c = {'ID' : lambda x: x in targets}
        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['transcriptOverlap', 'mismatchStatus', 'centerExpression'], conditions = c)

        for oID in oNX.filteredTargets:
                oNX.filteredTargets[oID] = []
                for aID in oNX.targets[oID]:

                        #transcriptOverlap
                        if inTranscript:
                                if not aNX.transcriptOverlap[aID]:
                                        #print 'tOverlap Fail', cgAlignment.pretty#print(alignment)
                                        continue

                        #misLevel
                        if aNX.mismatchStatus[aID][misLevel]:
                                #print 'mismatch Fail', cgAlignment.pretty#print(alignment)
                                continue
                           
                        #centerLevel
                        if aNX.centerExpression[aID][centerLevel] < minCenterLevel:
                                #print 'expression Fail', cgAlignment.pretty#print(alignment)
                                continue
                        
                        oNX.filteredTargets[oID].append(aID)

        oNX.save()                       

def filterTargetsInPlace(aFN, inTranscript, misLevel, centerLevel, minCenterLevel, rn = None, tn = None):
        if inTranscript == 'True': inTranscript = True
        if inTranscript == 'False': inTranscript = False
        misLevel, centerLevel, minCenterLevel =  int(misLevel), int(centerLevel), float(minCenterLevel)


        aNX = cgNexusFlat.Nexus(aFN, cgAlignmentFlat.cgAlignment)
        aNX.load(['passed', 'transcriptOverlap', 'mismatchStatus', 'centerExpression'], [rn, tn])

        for aID in aNX.passed:

                aNX.passed[aID] = False

                #transcriptOverlap
                if inTranscript:
                        if not aNX.transcriptOverlap[aID]:
                                #print 'tOverlap Fail', cgAlignment.pretty#print(alignment)
                                continue

                #misLevel
                if aNX.mismatchStatus[aID][misLevel]:
                        #print 'mismatch Fail', cgAlignment.pretty#print(alignment)
                        continue
                   
                #centerLevel
                if aNX.centerExpression[aID][centerLevel] < minCenterLevel:
                        #print 'expression Fail', cgAlignment.pretty#print(alignment)
                        continue
                
                aNX.passed[aID] = True

        aNX.save()                       




if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
