import bioLibCG

def targetCount(fN):
       
       targetDict = {}
       f = open(fN, 'r')
       for line in f:
               ls = line.strip().split('\t')
               targets = ls[4].split(',')
               for target in targets:
                        if target in targetDict:
                                targetDict[target] += 1
                        else:
                                targetDict[target] = 1

       for target in targetDict:
               print target + '\t' + str(targetDict[target])


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(targetCount, sys.argv)
