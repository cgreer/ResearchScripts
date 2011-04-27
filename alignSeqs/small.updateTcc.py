import bioLibCG

def updateTcc(peakFN, aDir):
        
        fOut = open(aDir + '/a.tcc.string.none', 'w')
        f = open(peakFN, 'r')
        for i,line in enumerate(f):
                ls = line.strip().split('\t')
                tcc = ls[0] #tcc may need to be shortened?
                fOut.write('%s\t%s\n' % (i, tcc))


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(updateTcc, sys.argv)
