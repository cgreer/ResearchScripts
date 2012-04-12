import bioLibCG
import cgOriginRNAFlat
import cgNexusFlat


def updateTotalContig(oFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['sequence', 'totalContigLength'], [rn, tn])

        for oID in oNX.sequence:
                seq = oNX.sequence[oID] 

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

                oNX.totalContigLength[oID] = highestLength

        oNX.save()
        
def updateEndContig(oFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['sequence', 'endContigLength'], [rn, tn])
        
        for oID in oNX.sequence:
                seq = oNX.sequence[oID]

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

                oNX.endContigLength[oID] = highest                        
               
        oNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateEndContig, sys.argv)
        bioLibCG.submitArgs(updateTotalContig, sys.argv)
