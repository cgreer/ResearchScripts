import bioLibCG
import cgOriginRNAFlat
import cgNexusFlat

def updateSeqDuplicateMultiTcc(oFN, rn = None, tn = None):
        
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['sequence', 'sequenceDuplicate', 'tcc', 'tccs'], [rn, tn])
        
        knownSeq_oID = {}

        for oID in oNX.sequence:
                oNX.sequenceDuplicate[oID] = False
                if oNX.sequence[oID] in knownSeq_oID:
                        oNX.sequenceDuplicate[oID] = True
                        otherOID = knownSeq_oID[oNX.sequence[oID]]
                        oNX.tccs[otherOID].append(oNX.tcc[oID])
                else:
                        knownSeq_oID[oNX.sequence[oID]] = oID

        oNX.save()

def updateSharedTargets(oFN, rn = None, tn = None):
        '''Just because there are duplicate sequences does not mean that
        the genomic position of each results is the correct one.  The 
        targets for each genomic position should be the same as the targets
        for each duplicate sequence

        make set of targets for each oID --> set each oid's targets'''

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['sequence', 'filteredTargets'], [rn, tn])

        knownSeq_targets = {}        

        #create oID groups and target sets.
        for oID in oNX.sequence:
                currSeq = oNX.sequence[oID]
                        

                #add targets to set
                for tID in oNX.filteredTargets[oID]:
                        knownSeq_targets.setdefault(currSeq, set()).add(tID)

                        
        for oID in oNX.sequence:

                currSeq = oNX.sequence[oID]

                newTargets = list(knownSeq_targets.get(currSeq, set()))
                oNX.filteredTargets[oID] = newTargets

        oNX.save()                
        


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateSeqDuplicateMultiTcc, sys.argv)
