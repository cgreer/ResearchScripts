import bioLibCG

def idTargets(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                numTargets = len(ls[4].split(','))

                print id + '\t' + str(numTargets)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(idTargets, sys.argv)
