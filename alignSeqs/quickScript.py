import bioLibCG
import cgDB
import cgOriginRNA

def quickScript(oDir):
        
        oDC = cgDB.dataController(oDir, cgOriginRNA.OriginRNA)
        id_oRNA = oDC.load()

        failedList = []

        for id, oRNA in id_oRNA.items():

                if oRNA.endContigLength > 7 or oRNA.totalContigLength > 7:
                        failedList.append(id)

        
        for id in failedList:
                print id

def removeAlignments(aFN, idFN, outFN):

        failedIDs = []
        f = open(idFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                failedIDs.append(ls[0])

        fOut = open(outFN, 'w')        
        f = open(aFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                if ls[0] in failedIDs:
                        continue

                fOut.write(line)
                




if __name__ == "__main__":
        import sys
        #bioLibCG.submitArgs(quickScript, sys.argv)
        bioLibCG.submitArgs(removeAlignments, sys.argv)
