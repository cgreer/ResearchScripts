import bioLibCG
import cgNexusFlat
import numpy as np
import matplotlib.pyplot as plt

def byPercentageBar(countFN, snrFN, imgName):

        N = 4
        ind = np.arange(N)    # the x locations for the groups
        width = 0.1       # the width of the bars: can also be len(x) sequence

        
        mm_counts = []
        f = open(countFN, 'r')
        for i, line in enumerate(f):
            ls = line.strip().split('\t')
            mm_counts.append([])
            for count in ls[:4]:
                mm_counts[-1].append(int(count))
        f.close()
        
        mm_snrs = []
        f = open(snrFN, 'r')
        for i, line in enumerate(f):
            ls = line.strip().split('\t')
            mm_snrs.append([])
            for snr in ls[:4]:
                mm_snrs[-1].append(float(snr))
        f.close()


        plt.subplot(211)
        
        for i in [0,1,2,3,4,5]:
            nums = mm_counts[i]
            print nums 
            plt.bar(.1 + ind + i*width, nums, width, color='w')
            
            #mismatch text
            for a, xloc in enumerate(ind):
                    plt.text(.14 + xloc + i*width, nums[a] + .2, str(i))

        plt.title('Filtered Results by Parameters')
        plt.xticks(ind+width/2. + .4, ('50', '60', '70', '80') )
        plt.ylabel('Number of oRNA w/ SNR > 2 by Mismatch')

        plt.subplot(212)
        for i in [0,1,2,3,4,5]:
            nums = mm_snrs[i]
            print nums 
            plt.bar(.1 + ind + i*width, nums, width, color='w')
            
            #mismatch text
            for a, xloc in enumerate(ind):
                    plt.text(.14 + xloc + i*width, nums[a] + .2, str(i))

        plt.xticks(ind+width/2. + .4, ('50', '60', '70', '80') )
        plt.xlabel('Percentage of Middle Degradome Expression')
        plt.ylabel('Averatge SNR of all oRNAs w/ SNR > 2 by Mismatch')

        
        plt.savefig(imgName, bbox_inches='tight', pad_inches=1)



if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
