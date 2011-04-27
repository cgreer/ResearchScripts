import bioLibCG

def getSignificantIDs(rFN, avgNumFN):
        
        id_avgNum = {}
        f = open(avgNumFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                avgNum = float(ls[1])
                id_avgNum[id] = avgNum

        f = open(rFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                numTargets = len(ls[4].split(','))
                try:
                        avgNum = id_avgNum[id]
                except KeyError:
                        avgNum = .01

                if float(numTargets) > avgNum:
                        print line.strip()



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getSignificantIDs, sys.argv)
