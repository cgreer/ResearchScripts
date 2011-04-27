import bioLibCG
import matplotlib.pyplot as plt

def numTargets(histID):
        
        histID = int(histID)

        totalTargets = {}
        for i in range(0, 2084):
                totalTargets[i] = []

        for j in range(1, 100):
                tarDict = {}
                for i in range(0, 2084):
                        tarDict[i] = 0
                f = open('simulations/simulation.%s/simulation.%s.final.results' % (j,j), 'r')
                for line in f:
                        ls = line.strip().split('\t')
                        id = int(ls[0])
                        numTargets = len(ls[4].split(','))

                        tarDict[id] = numTargets

                for id in tarDict:
                        totalTargets[id].append(tarDict[id])
                
                f.close()
        
        #histVals = totalTargets[histID]
        #plt.hist(histVals, 30, facecolor='b', alpha = .75)

        #plt.show()
        

        for id in totalTargets:
                                 
                numTargets = sum(totalTargets[id])
                if numTargets > 0:
                        avgTargets = float(numTargets)/99
                        print '%s\t%s' % (id, avgTargets)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(numTargets, sys.argv)
