import bioLibCG

def SNR(resultsFN, simAvgFN):
        
        id_avgNumTargets = {}
        f = open(simAvgFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                numTargets = float(ls[1])
                id_avgNumTargets[id] = numTargets

        f = open(resultsFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                targets = bioLibCG.uniqueList(ls[4].split(','))
                numTargets = len(targets)
                try:
                        avgNum = id_avgNumTargets[id]
                except KeyError:
                        avgNum = .001
                SNR = float(numTargets)/avgNum

                print '%s\t%s\t%s\t%s' % (id, SNR, numTargets, avgNum)


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(SNR, sys.argv)
