import bioLibCG
import cgDB
import cgAlignment

def probeAlignments(aDir):
        
        probePairs = [[6, 35934]]

        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()

        for alignment in id_alignment.values():
                for sID, tID in probePairs:

                        if alignment.sID == sID and alignment.tID == tID:
                                print alignment.id, alignment.sID, alignment.tID, alignment.centerExpression, alignment.mismatchStatus, alignment.numMismatches, alignment.transcriptOverlap

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(probeAlignments, sys.argv)
