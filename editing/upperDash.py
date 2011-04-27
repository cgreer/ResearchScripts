import bioLibCG

def upperDash(fN ):

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                field = ls[0]
                field = field.upper()
                
                print '%s\t%s' % (field, ls[1])

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(upperDash, sys.argv)

