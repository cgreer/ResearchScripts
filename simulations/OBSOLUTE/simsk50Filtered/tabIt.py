import bioLibCG

def tabIt(fN, outFN):

        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.split(' ')
                fOut.write('\t'.join(ls))
        f.close()
        fOut.close()
                
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
