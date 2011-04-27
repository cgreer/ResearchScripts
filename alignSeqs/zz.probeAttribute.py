import bioLibCG
import cgDB
import cgOriginRNA

def probeMicro(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():

                if oRNA.passedFilter:
                        print oRNA.id, oRNA.sequence, oRNA.tcc, oRNA.tccs


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(probeMicro, sys.argv)
