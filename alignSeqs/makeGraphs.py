import bioLibCG
import matplotlib.pyplot as plt

def plotEntropyTargets(rOneFN, rTwoFN, compFN):
        
        id_comp = {}
        f = open(compFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id_comp[int(ls[0])] = float(ls[1])


        f = open(rOneFN)

        complexities = []
        seqNumTargets = []

        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                complexity = id_comp[id]
                numTargets = len(ls[4].split(','))
                complexities.append(complexity)
                seqNumTargets.append(numTargets)

        plt.plot(complexities, seqNumTargets, 'ro', label = 'HeLa sRNA', color = 'red')

        f.close()
        
        f = open(rTwoFN)

        complexities = []
        seqNumTargets = []

        for line in f:
                ls = line.strip().split('\t')
                id = int(ls[0])
                complexity = id_comp[id]
                numTargets = float(ls[1])
                complexities.append(complexity)
                seqNumTargets.append(numTargets)

        plt.plot(complexities, seqNumTargets, 'bo', label = 'simulated sRNA', color = 'blue')
        plt.legend()
        plt.ylabel('Number of targets per small RNA (filter: O >0-NoMicroTran, T YesTran4Mis > .55 6bp')
        plt.xlabel('Complexity of smallRNA')
        plt.title('Origin RNA Target simulation')
        plt.show()
        f.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotEntropyTargets, sys.argv)

        
