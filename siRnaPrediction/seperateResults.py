import bioLibCG
import compareData

def reportDifference(oldFN, newFN):

        f = open(oldFN, 'r')
        oldList = [x for x in f]
        oldCoords = [x.strip().split('\t')[0] for x in oldList]
        f.close()

        f = open(newFN, 'r')
        newList = [x for x in f]
        newCoords = [x.strip().split('\t')[0] for x in newList]
        f.close()
        
        bothList = []
        bothList.extend(newList)
        bothList.extend(oldList)

        bothCoords = compareData.compareTwoTcc(oldCoords, newCoords, 1)

        onlyOld = []
        for x in oldCoords:
                if x not in bothCoords:
                        onlyOld.append(x)

        onlyNew = []
        for x in newCoords:
                if x not in bothCoords:
                        onlyNew.append(x)


        print 'Both'
        knownList = []
        for x in bothCoords:
                for y in bothList:
                        if x in y:
                                print y.strip()
                                knownList.append(y)
                                break

        print 'old'
        for y in oldList:
                if y not in knownList:
                        print y.strip()

        print 'new'
        for y in newList:
                if y not in knownList:
                        print y.strip()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(reportDifference, sys.argv)
