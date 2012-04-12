import bioLibCG
import cgNexusFlat
import cgOriginRNAFlat
import cgWig

def updateRepeatStatus(oFN, wigDir, chrom, strand):

        #load oRNAs
        oNX = cgNexusFlat.Nexus(oFN, cgOriginRNAFlat.OriginRNA)
        oNX.load(['repeatStatus', 'tcc'])
        
        #load wig file for chrom, strand
        coord_value = cgWig.loadSingleWig(wigDir, chrom, strand, 'REPEAT')


        for oID in oNX.repeatStatus:

                oChrom, oStrand, start, end = bioLibCG.tccSplit(oNX.tcc[oID])
                if oChrom != chrom or oStrand != strand:
                        continue

                oNX.repeatStatus[oID] = False
                for i in range(start, end + 1):
                        if i in coord_value:
                                oNX.repeatStatus[oID] = True
                                break


        oNX.save()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
