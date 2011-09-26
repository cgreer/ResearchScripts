import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import cgOriginRNAFlat
import compareData
from cgAlignmentFlat import cgAlignment

def getCoords(fN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(fN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tcc'], [rn, tn])

        for id in oNX.tcc:
                chrom, strand, start, end = bioLibCG.tccSplit(oNX.tcc[id])

                name = 'None'
                score = '0'
                strand = bioLibCG.switchStrandFormat(strand)
                thickStart = start
                thickEnd = end

                pString = [str(x) for x in [chrom, start, end, name, score, strand, thickStart, thickEnd]]
                print '\t'.join(pString)

def getCoordsAlignments(fN, rn = None, tn = None):

        aNX = cgNexusFlat.Nexus(fN, cgAlignment)
        aNX.load(['tTcc'], [rn, tn])

        for id in aNX.tTcc:
                chrom, strand, start, end = bioLibCG.tccSplit(aNX.tTcc[id])

                name = 'None'
                score = '0'
                strand = bioLibCG.switchStrandFormat(strand)
                thickStart = start
                thickEnd = end

                pString = [str(x) for x in [chrom, start, end, name, score, strand, thickStart, thickEnd, id]]
                print '\t'.join(pString)
        
def getTccs(fN):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                chrom, strand, start, end = ls[0], bioLibCG.switchStrandFormat(ls[5]), ls[1], ls[2]

                print bioLibCG.makeTcc(chrom, strand, start, end)

def compareTccs(humanFN, liftCoords, rn = None, tn = None):

        mouseList = []
        f = open(liftCoords, 'r')
        for line in f:
                ls = line.strip().split('\t')
                mouseList.append(ls[0])

        humanList = []
        DC = cgNexusFlat.Nexus(humanFN, cgOriginRNAFlat.OriginRNA)
        DC.load(['tcc'], [rn, tn])
        for id in DC.tcc:
                humanList.append(DC.tcc[id])

        mouseList = list(set(mouseList))
        humanList = list(set(humanList))

        x = compareData.compareTwoTcc(humanList, mouseList)

        for i in x:
                print i

def compareTccs(humanFN, liftCoords, rn = None, tn = None):
        '''compare if there is an overlap between a mouse alignment and any human''' 
        mouseList = []
        f = open(liftCoords, 'r')
        for line in f:
                ls = line.strip().split('\t')
                mouseList.append(ls[0])

        humanList = []
        DC = cgNexusFlat.Nexus(humanFN, cgOriginRNAFlat.OriginRNA)
        DC.load(['tcc'], [rn, tn])
        for id in DC.tcc:
                humanList.append(DC.tcc[id])

        mouseList = list(set(mouseList))
        humanList = list(set(humanList))

        x = compareData.compareTwoTcc(humanList, mouseList)

        for i in x:
                print i


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
