import bioLibCG
import matplotlib.pyplot as plt

def quickHisto(histVals, xlabel, ylabel = 'number of X', title = 'Histogram'):


        plt.hist(histVals, 30, facecolor='b', alpha = .75)
        plt.title(title)
	plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(name, sys.argv)
