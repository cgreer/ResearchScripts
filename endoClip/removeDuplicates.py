import bioLibCG

def removeDuplicates(fN, outFN):

        f = open(fN, 'r')
        
        seqs = {}
        for seq in f:
                seqs[seq] = 1
        f.close()

        outF = open(outFN, 'w')
        
        for seq in seqs:
                 outF.write(seq)
        
        
        f.close()

def countDuplicated(fN):

        uniqueSeqs = set()
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                uniqueSeqs.add(ls[1])

        print 'number of unique', len(uniqueSeqs)                

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])



