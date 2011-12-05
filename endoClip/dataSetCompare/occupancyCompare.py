import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
import cgDegPeak

def makeHitMapDegPeak(dFN, switchStrand = False):

    NX = cgNexusFlat.Nexus(dFN, cgDegPeak.Peak)
    NX.load(['tcc', 'eLevel'])

    c_s_coord = {}

    for id in NX.ids:

        chrom, strand, start, end = bioLibCG.tccSplit(NX.tcc[id])
        if switchStrand:
            strand = bioLibCG.switchStrand(strand)

        for i in range(start, end + 1):
            c_s_coord.setdefault(chrom, {}).setdefault(strand, set()).add(i)

    return c_s_coord

def straightOccupancy(dFN1, dFN2, limitContext = "INTER"):
    

    print "loading data sets"
    NX = cgNexusFlat.Nexus(dFN1, cgDegPeak.Peak)
    NX.load(['tcc', 'eLevel', 'context'])
    
    twoChrom_strand_coord = makeHitMapDegPeak(dFN2, True)



    #calculate the # of peaks in each size category
    print "calculating max coverage"
    theMax = 100
    minSize_count = dict( (i, 0) for i in range(theMax + 1))

    for id in NX.ids:
        if NX.context[id] == limitContext:
            for i in range(theMax + 1):
                if NX.eLevel[id] >= i:
                    minSize_count[i] += 1
        
    
    # calculate overlap for each size
    minSize_overlap = dict( (i, 0) for i in range(0, theMax + 1))
    print "calculating % coverage"
    for id in NX.ids:
        if NX.context[id] == limitContext:
            chrom, strand, start, end = bioLibCG.tccSplit(NX.tcc[id])
            for i in range(start, end + 1):
                if i in twoChrom_strand_coord[chrom][strand]:
                    for i in range(theMax + 1):
                        if NX.eLevel[id] >= i:
                            minSize_overlap[i] += 1
           
    # print           
    for i in range(theMax + 1):
        try:
            #print '%s: %s / %s || overall coverage: %s' % (i, minSize_overlap[i], minSize_count[i], float(minSize_overlap[i])/minSize_count[i])
            print '%s\t%s' % (i, float(minSize_overlap[i])/minSize_count[i]), minSize_count[i]
        except ZeroDivisionError:
            print '0\t0'


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
