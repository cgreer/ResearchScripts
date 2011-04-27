import bioLibCG

def filterDups(fN, oFN):

        outF = open(oFN, 'w')
        knownSeqs = [] 
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[1]
                if seq not in knownSeqs:
                        outF.write(line)
                        knownSeqs.append(seq)
        
        f.close()
        outF.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterDups, sys.argv)
