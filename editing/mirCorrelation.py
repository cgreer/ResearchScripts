import bioLibCG
import matplotlib.pyplot as plt
import math

def corr(fN):

        before = []
        after = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                before.append(math.log(float(ls[2]),2))
                after.append(math.log(float(ls[1]),2))

        
        
        plt.title('RPKM of mir-324-3p before/after editing')
        plt.xlabel('log2(rpkm)')
        plt.ylabel('Number of Genes')
        plt.hist(before, 10, cumulative = True, histtype = 'step', normed = True, label = 'Unedited')
        plt.hist(after, 10, cumulative = True, histtype = 'step', normed = True, label = 'Edited')
        plt.legend()
        plt.show()
                
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(corr, sys.argv)
                             
