import bioLibCG

def convertTab(aFN, outFN):
        
        fOut = open(outFN, 'w')
        f = open(aFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                fOut.write('%s\n' % '\t'.join(ls))
        f.close()
        fOut.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(convertTab, sys.argv)
