import bioLibCG
import cgNexusFlat
import matplotlib.pyplot as plt

def plotData(fN, type, end):
        
        y = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                num = int(ls[1])
                y.append(num)

                 
        plt.title('%s %s' % (type, end))
        plt.plot(y)
        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
