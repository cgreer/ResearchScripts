import bioLibCG
import cgGenes3
import dumpObj
import time
import sys


def doit(fN):

        for i in xrange(10):
                sys.stdout.write('\r%s' % i)
                time.sleep(.2)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(doit, sys.argv)
                        



