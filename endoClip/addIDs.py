import bioLibCG

def addIDs(fN, outFN):
        
        fOut = open(outFN, 'w')
        i = 0
        f = open(fN, 'r')
        for line in f:
                fOut.write('%s\t%s\n' % (i, line.strip()))
                i += 1

        fOut.close()


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(addIDs, sys.argv)
