import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit
import cgWig

def updateTargetIDs(oFN, oFF, aFN, aFF):

    NX = Nexus(oFN, oFF)
    NX.load(['filteredTargets'])

    aNX = Nexus(aFN, aFF)
    aNX.load(['sID'])


    while aNX.nextID():

        NX.id = aNX.sID
        NX.filteredTargets.append(aNX.id)

    NX.save()

def updateRepeatStatus(fN, fF, wigDir, chrom, strand):

    #load oRNAs
    NX = Nexus(fN, fF)
    NX.load(['repeat', 'tcc'])
    
    #load wig file for chrom, strand
    coord_value = cgWig.loadSingleWig(wigDir, chrom, strand, 'REPEAT')

    while NX.nextID():
        oChrom, oStrand, start, end = bioLibCG.tccSplit(NX.tcc)
        if oChrom != chrom or oStrand != strand: continue

        NX.repeat = False
        for i in range(start, end + 1):
            if i in coord_value:
                NX.repeat = True
                break

    NX.save()

@autocast
def updateContext(fN, fF, wigDir, chrom, strand, switchStrand = False):
        
    NX = Nexus(fN, fF)
    NX.load(['tcc', 'context'])
    
    if switchStrand:
        strand = str(-int(strand))
    else:
        strand = str(strand)
    
    print 'loading wig'
    coord_contexts = cgWig.loadSingleWigContext(wigDir, chrom, strand, 'context') 
    print 'done loading'

    ds = bioLibCG.dominantSpotter(['C_EXON', 'C_3UTR', 'C_5UTR', 'NC_EXON', 'NC_3UTR', 'NC_5UTR', 'C_INTRON', 'NC_INTRON', 'INTER']) 


    while NX.nextID():

        oChrom, oStrand, start, end = bioLibCG.tccSplit(NX.tcc)
        
        #deg wigs is AS to actual clipping site
        if switchStrand:
            oStrand = str(-int(strand))
        else:
            oStrand = str(oStrand)

        if oChrom == chrom and oStrand == strand:

            contexts = coord_contexts.get(start, 'INTER').split(',')
            NX.context = ds.spotItem(contexts)

    
    NX.save()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

