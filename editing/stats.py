import bioLibCG

def countUniqueID(fN):

        f = open(fN, 'r')

        countDict = {}

        for line in f:
                id = line.strip().split('\t')[0]
                countDict[id] = line.strip()
        f.close()

        print len(countDict)
        for i in countDict:
                print countDict[i]


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(countUniqueID, sys.argv)

