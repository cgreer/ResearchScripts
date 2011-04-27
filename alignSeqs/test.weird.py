import bioLibCG
import cgDB
import cgOriginRNA

def test(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        print id_oRNA[1].targets
        id_oRNA[1].targets.append(13)
        print id_oRNA[1].targets
        print id_oRNA[2].targets

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(test, sys.argv)
