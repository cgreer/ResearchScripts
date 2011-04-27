import cgDB
import cgOriginRNA
import bioLibCG

def probeORNA(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():
                if oRNA.passedFilter:
                        cgOriginRNA.prettyPrint(oRNA)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(probeORNA, sys.argv)
