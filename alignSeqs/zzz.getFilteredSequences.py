import bioLibCG
import cgDB
import cgOriginRNA

def getSeqs(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for id, oRNA in id_oRNA.items():
                
                if oRNA.sequenceDuplicate:
                        continue
                if oRNA.totalContigLength > 6:
                        continue
                if oRNA.endContigLength > 6:
                        continue

                print '%s' % id

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getSeqs, sys.argv)
