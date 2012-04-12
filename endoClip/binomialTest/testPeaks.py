import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
from scipy.stats import binom

@autocast
def testPeaks(degFN, dForm, allGeneInfo, gForm, switchStrand = False):

    #load/configure gene Info
    gNX = Nexus(allGeneInfo, gForm)
    gNX.load(['geneName', 'numReads', 'numSpots'])
    
    gName_numReads = {}
    gName_numSpots = {}
    while gNX.nextID():
        gName_numReads[gNX.geneName] = gNX.numReads
        gName_numSpots[gNX.geneName] = gNX.numSpots
   
   
    #load degFN info
    dNX = Nexus(degFN, dForm)
    dNX.load(['tcc', 'eLevel', 'geneNames', 'pValBin'])

    while dNX.nextID():
       
        gNames, readsForPeak = dNX.geneNames, dNX.eLevel
        chrom, strand, start, end = bioLibCG.tccSplit(dNX.tcc)
        if switchStrand:
            strand = -int(strand)
      
        pVals = []
        for gName in gNames:
            
            #may have to change gene name cuz of multiple spans
            try:
                totGeneReads = gName_numReads[gName]
                numSpotsForGene = gName_numSpots[gName]
            except KeyError:

                try:
                    gName = gName + '_RE_%s_%s' % (chrom, strand)
                    totGeneReads = gName_numReads[gName]
                    numSpotsForGene = gName_numSpots[gName]
                except KeyError:
                    print "FIX THIS GENE NAME", gName
                    continue

            #add psuedocount
            totGeneReads += 1
            numSpotsForGene += 1 # not sure whether to do this yet...

            #check for hidden intron gene overlap
            try:
                q = 1.0/numSpotsForGene
            except ZeroDivisionError:
                continue #intron gene

            #add p val
            pVals.append(binom.sf(readsForPeak, totGeneReads, q))

        dNX.pValBin = max(pVals) if pVals else -1.0

    dNX.save()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

