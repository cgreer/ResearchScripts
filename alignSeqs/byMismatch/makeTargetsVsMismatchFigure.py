import bioLibCG
import cgPlot

def makeGraph(baseName, lowRange, highRange):
        lowRange = int(lowRange)
        highRange = int(highRange)
        
        
        histDict = {}
        for i in range(lowRange, highRange + 1):
        
                f = open(baseName + str(i), 'r')
                for line in f:
                        targets = line.strip().split('\t')[4]
                        numTargets = len(targets.split(','))
                        try:
                                histDict[i].append(numTargets)
                        except KeyError:
                                histDict[i] = [numTargets]

                f.close()
        print histDict 
        
        for m in histDict:
                cgPlot.plotHistogram(histDict[m], name = str(m) + ' mismatches')
               
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeGraph, sys.argv)

