import bioLibCG
import cgNexusFlat

def countRepeats(targetList, alignmentFile, outFN):

        targets = set()
        f = open(targetList, 'r')
        for line in f:
                ls = line.strip().split('\t')
                targets.add(ls[0])
        f.close()

        nonCount = 0
        rCount = 0
        f = open(alignmentFile, 'r')
        for line in f:
                ls = line.strip().split('\t')
                aID = ls[0]
                rStatus = ls[18]

                if aID in targets:
                        if rStatus == 'T':
                                rCount += 1
                        else:
                                nonCount +=1
        f.close()

        f = open(outFN, 'w')
        f.write('repeat: %s, nonRepeat: %s' % (rCount, nonCount))
        f.close()

                                        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
