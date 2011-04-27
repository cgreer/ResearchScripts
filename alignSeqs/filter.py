import bioLibCG

def filter(fN, oFN):

        f = open(fN, 'r')
        fOut = open(oFN, 'w')

        for line in f:
                microStatus = int(line.strip().split('\t')[5])

                if microStatus == 1:
                        fOut.write(line)

        f.close()
        fOut.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filter, sys.argv)



