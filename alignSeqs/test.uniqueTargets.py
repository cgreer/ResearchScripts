import bioLibCG

def uniqueTargets(fN):
        
        uniqueTargets = set()
        numOrigins = 0
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                numOrigins += 1
                targets = ls[4].split(',')
                for target in targets:
                        if target not in uniqueTargets:
                                uniqueTargets.add(target)
        
        numTargets = len(uniqueTargets)
        print 'targets:',  numTargets
        print 'origin RNAs:', str(numOrigins)
        print 'target/oRNA:', str(float(numTargets)/numOrigins)                                

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(uniqueTargets, sys.argv)
