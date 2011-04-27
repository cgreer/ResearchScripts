import bioLibCG
import matplotlib.pyplot as plt
import numpy as np

def plotBarSNR():


        N = 10
        menMeans = (61324, 51882, 45661, 40722, 29312, 32851, 27770, 24162, 21578, 15677)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.17       # the width of the bars

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, menMeans, width, color='r')

        womenMeans = (49344, 43652, 40213, 36490, 25724, 6236, 5237, 4639, 4108, 2774)
        womenStd =   (220,209,202,191,152,13,12,11,11,9)
        rects2 = ax.bar(ind+width, womenMeans, width, color='y', yerr=womenStd)

        # add some
        ax.set_ylabel('Total Number of Targets')
        ax.set_xlabel('descriptor (a:b:c)\n a = # bp from center where NO mismatches allowed\n b = # of bp from center where at least c% of degradome expression must be found')
        ax.set_title('Total SNR')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( ('4.6.30', '4.6.50', '4.6.60', '4.6.70', '4.6.90', '6.6.30', '6.6.50', '6.6.60', '6.6.70', '6.6.90'))
        ax.axis([0,10,0,70000])

        ax.legend( (rects1[0], rects2[0]), ('Observed', 'Simulated') )

        def autolabel(rects):
                # attach some text labels
                for rect in rects:
                        height = rect.get_height()
                        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height), ha='center', va='bottom')

        autolabel(rects1)
        autolabel(rects2)

        plt.show()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotBarSNR, sys.argv)
