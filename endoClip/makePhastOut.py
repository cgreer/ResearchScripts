import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat

def makePhastOut(oFN, outFN, rn = None, tn = None):

        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['tcc'], [rn, tn])

        outF = open(outFN, 'w')
        
        for oID in oNX.tcc:
                chrom, strand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                
                outString = '\t'.join([str(x) for x in [oID, chrom, start, end]]) + '\n'
                outF.write(outString)

        outF.close()                
        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
