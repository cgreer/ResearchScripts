import bioLibCG

def percentage(contextFN, numSites):

        numSites = int(numSites)
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                numType = int(ls[1])
                print line.strip() + '\t' + str(100 * float(numType)/numSites)

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(percentage, sys.argv)
