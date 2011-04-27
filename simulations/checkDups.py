import bioLibCG

def checkDups(dir):
        
        id_seqs = {}
        for i in bioLibCG.recurseDir(dir, end = '.simSeqs'):
              
                f = open(i, 'r')
                for line in f:
                        ls = line.strip().split('\t')
                        id, seq = ls[0], ls[1]
                        id_seqs.setdefault(id, []).append(seq)

        
        for id, seqs in id_seqs.items():
                if len(seqs) != len(list(set(seqs))):
                        print 'fail'
                        print id, seqs



if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(checkDups, sys.argv)
