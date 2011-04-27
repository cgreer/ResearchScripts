import bioLibCG

def dupLines(fN):

        dup = {}
        f = open(fN, 'r')
        for line in f:
                dup[line.strip()] = 1


        for k in dup:
                print k

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(dupLines, sys.argv)
                
