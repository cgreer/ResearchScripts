import bioLibCG

def formatFile(fN):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                newList = []
                newList.extend(ls[1:11])
                newList.extend(ls[12:])
                print '\t'.join(newList)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(formatFile, sys.argv)

