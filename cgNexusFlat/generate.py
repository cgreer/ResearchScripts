import bioLibCG

def generate(numRows):
        
        for i in xrange(0,int(numRows)):
                print '%s\t%s\t%s' % (i, 'ATATATATATATA', '0')

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
