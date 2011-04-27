import bioLibCG
import cgDB
import cgOriginRNA

def updateTotalContig(oDir):

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        for oRNA in id_oRNA.values():
                seq = oRNA.sequence

                highestLength = 1
                cLength = 1
                letters = list(seq)
                for i,letter in enumerate(letters):
                        if i == 0: continue

                        if letters[i] == letters[i-1]:
                                cLength += 1
                                if cLength > highestLength:
                                        highestLength = cLength
                        else:
                                cLength = 1

                oRNA.totalContigLength = highestLength

        oDC.commit(id_oRNA)                
        
def updateEndContig(oDir):

        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()
        
        for oRNA in id_oRNA.values():
                seq = oRNA.sequence

                #5'
                cLength5 = 1
                for i,letter in enumerate(seq):
                        if i == 0: continue

                        if seq[i] == seq[i-1]:
                                cLength5 += 1
                        else:
                                break
                #3'
                cLength = 1
                revSeq = [x for x in reversed(seq)]
                for i,letter in enumerate(revSeq):
                        if i == 0: continue

                        if revSeq[i] == revSeq[i-1]:
                                cLength += 1
                        else:
                                break

                highest = cLength5
                if cLength > cLength5:
                        highest = cLength

                oRNA.endContigLength = highest                        
               
        oDC.commit(id_oRNA)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateEndContig, sys.argv)
        bioLibCG.submitArgs(updateTotalContig, sys.argv)
