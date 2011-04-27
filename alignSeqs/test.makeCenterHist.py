import bioLibCG
import matplotlib.pyplot as plt

def centerHist(fN, shift = 0):
        shift = int(shift)

        histVals = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                cLevel = float(ls[2 + shift])
                histVals.append(cLevel)
        plt.title('Degradome Expression in Center of Small RNA')
        plt.xlabel('Percentage of Degradome Expression Inside %s bps from middle' % (4 + 2*shift))
        plt.ylabel('Number of small RNA/degradation sequence pairs')
        plt.hist(histVals, 30, facecolor='b', alpha = .75)

        plt.show()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(centerHist, sys.argv)
