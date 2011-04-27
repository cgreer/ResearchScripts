import bioLibCG
import matplotlib.pyplot as plt

def countEnd(fN):
        
        count3 = {'A': [], 'T': [], 'G': [], 'C': []}
        count5 = {'A': [], 'T': [], 'G': [], 'C': []}
        countT = {'A': [], 'T': [], 'G': [], 'C': []}

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[1]
                
                letter0 = seq[-1]
                count = 0
                for i in reversed(seq):
                        if i == letter0:
                                count += 1
                        else:
                                break
                count3[letter0].append(count)


                letter = seq[0]
                countAnother = 0
                for i in seq:
                        if i == letter:
                                countAnother += 1
                        else:
                                break
                count5[letter].append(countAnother)
                
                if count > countAnother:
                        countT[letter0].append(count)
                else:
                        countT[letter].append(countAnother)

                
                

        plt.hist(countT['C'], 15, facecolor='r', label='C', alpha = 1.00)
        plt.hist(countT['G'], 15, facecolor='y', label='G', alpha = .55)
        plt.hist(countT['T'], 15, facecolor='g', label='T', alpha = .55)
        plt.hist(countT['A'], 15, facecolor='b', label='A', alpha = .55)
        plt.xlabel('Length of Longest Contiguos End Region')
        plt.ylabel('Number of Origin RNAs')
        plt.legend()
        plt.show()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(countEnd, sys.argv)
