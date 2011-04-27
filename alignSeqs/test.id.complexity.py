import bioLibCG

def getComplexity(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                complexity = ls[7]
                print '%s\t%s' % (id, complexity)



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getComplexity, sys.argv)
