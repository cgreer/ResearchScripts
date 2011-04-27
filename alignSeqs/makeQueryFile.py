import bioLibCG

def makeQuery(fN, outFN):
        '''Given a file with raw sequences per line, convert to id seq'''
        f = open(fN, 'r')
        fOut = open(outFN, 'w')

        i = 0
        for line in f:
                seq = line.strip()
                fOut.write('%s\t%s\n' % (i, seq))
                i += 1

        f.close()
        fOut.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(makeQuery, sys.argv)

