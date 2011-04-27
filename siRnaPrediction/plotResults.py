import cgPlot
import bioLibCG

def plotResults(rFN, smallCName, degCName):

        f = open(rFN, 'r')
        
        i = 1
        for line in f:
               chrom, strand, start, end = bioLibCG.tccSplit(line.strip().split('\t')[0])
               start = start - 30
               end = end + 30

               newTcc = bioLibCG.makeTcc(chrom, strand, start, end)
               cgPlot.plotSmallDeg(newTcc, smallCName, degCName, 'newResults', line.strip(), str(i))
               i += 1

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotResults, sys.argv)
                                                                                
