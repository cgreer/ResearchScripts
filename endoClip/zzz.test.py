import bioLibCG
import cgNexusFlat
import cgWig

def test(tcc, wigDir):

        chrom, strand, start, end = bioLibCG.tccSplit(tcc)
        print chrom, strand
        coord_eLevel = cgWig.loadSingleWig(wigDir, chrom, strand, 'ALL')

        sKeys = sorted(coord_eLevel.keys())

        for i in range(start, end + 1):
                print i, coord_eLevel.get(i, 0)
        


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
