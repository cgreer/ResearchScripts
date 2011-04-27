import bioLibCG
import cgDB
import cgAlignment
import dumpObj

def testLoad(aDir):
        
        aDC = cgDB.dataController(aDir, cgAlignment.cgAlignment)
        id_alignment = aDC.load()

        dumpObj.dumpObj(id_alignment[2])

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(testLoad, sys.argv)
