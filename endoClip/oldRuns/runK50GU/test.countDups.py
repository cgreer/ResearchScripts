import bioLibCG
import gZipEntropy as gz

def countDups(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[3]
                numDups = len(ls[2].split(','))
                numT = len(ls[9].split(','))
                rStatus = ls[21]
                gEnt = gz.gZipEntropy(seq)

                print '\t'.join([str(numDups), seq, str(numT), rStatus, gEnt])
        f.close()                


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
