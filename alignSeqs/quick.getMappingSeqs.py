import cgDB
import cgOriginRNA
import bioLibCG

def getMappingSeqs(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():
                print '>%s' % oRNA.id
                print oRNA.sequence


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getMappingSeqs, sys.argv)
