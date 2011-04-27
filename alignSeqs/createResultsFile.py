import bioLibCG

def createResultsFile(peakFN, outFN):

        f = open(peakFN, 'r')
        peaks = [x.strip() for x in f]
        f.close()

        outF = open(outFN, 'w')
        for i, peak in enumerate(peaks):
                outF.write('%s\t%s\n' % (i, peak))

if __name__ == "__main__":
        import sys

        bioLibCG.submitArgs(createResultsFile, sys.argv)

