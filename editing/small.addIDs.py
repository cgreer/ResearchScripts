import bioLibCG

def addIDs(fN):
        
        f = open(fN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                
                print '%s\t%s' % (line.strip(), i)
                i += 1

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(addIDs, sys.argv)
