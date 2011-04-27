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

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateSeqDuplicateMultiTcc, sys.argv)
