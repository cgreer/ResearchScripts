import bioLibCG

def numTargets(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id, targets = ls[0], ls[4].split(',')
                uniq = []
                for target in targets:
                        if target not in uniq:
                                uniq.append(target)

                print id, ','.join(uniq)                               

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(numTargets, sys.argv)
