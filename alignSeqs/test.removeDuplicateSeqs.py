import bioLibCG

def removeDups(fN):
        
        uniqSeqs = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[1]
                if seq not in uniqSeqs:
                        uniqSeqs.append(seq)
                        print line.strip()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(removeDups, sys.argv)
