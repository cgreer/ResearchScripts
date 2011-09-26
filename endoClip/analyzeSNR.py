import bioLibCG
import cgNexusFlat
import numpy as np
import matplotlib.pyplot as plt

def byPercentageBar(fN):

        N = 4
        ind = np.arange(N)    # the x locations for the groups
        width = 0.1       # the width of the bars: can also be len(x) sequence

        zero = (6.09,6.25,5.0,4.62)
        one = (3.72,4.13,4.3,4.31)
        two = (3.44,3.75,3.6,3.53)
        three = (2.88898665493,2.88898665493,3.03225806452,3.02148973986)
        four = (2.65130372379,2.63807890223,2.80141028589,2.82434973363)
        five = (2.46168137483,2.38446547884,2.47972406078,2.44315366404)

        a = (10,9,6,5)
        b = (57, 55, 49, 47)
        c = (103, 100,107,110)
        d = (75,90,88,93)
        e = (35,44,45,50)
        f = (21,30,30,34)

        zero = (0,0,0,0)
        one = (0,0,0,0)
        two = (25.0, 21.42, 25.0, 50.0)
        three = (5.9, 6.3, 8.4, 9.1)
        four = (4.1,4.0,4.6,4.8)
        five = (3.3,3.2,3.2,3.3)

        a = (0,0,0,0) #0(50), 0(60), etc.
        b = (0, 0, 0, 0)
        c = (4, 3, 2, 1)
        d = (21,20,17,16)
        e = (52,50,47,48)
        f = (46,44,45,50) #5
        
        plt.subplot(211)
        
        for i, nums in enumerate([a,b,c,d,e,f]):
                print nums
                plt.bar(.1 + ind + i*width, nums, width, color='w')
                
                #mismatch text
                for a, xloc in enumerate(ind):
                        plt.text(.14 + xloc + i*width, nums[a] + .2, str(i))

        plt.title('Filtered Results by Parameters')
        plt.xticks(ind+width/2. + .4, ('50', '60', '70', '80') )
        plt.ylabel('Number of oRNA w/ SNR > 2 by Mismatch')

        plt.subplot(212)
        for i, nums in enumerate([zero, one, two, three, four, five]):
                plt.bar(.1 + ind + i*width, nums, width, color='w')
                
                for a, xloc in enumerate(ind):
                        plt.text(.14 + xloc + i*width, nums[a] + .2, str(i))

        plt.xticks(ind+width/2. + .4, ('50', '60', '70', '80') )
        plt.xlabel('Percentage of Middle Degradome Expression')
        plt.ylabel('Total SNR of all oRNAs w/ SNR > 2 by Mismatch')

        plt.show()                


def allSNRBar(fN):
        
        N = 4
        
        zero = (6.01, 6.19, 6.45, 7.5)
        one = (4.69, 5.02, 4.93, 5.73)
        two = (3.28, 3.66, 3.42, 3.33)
        three = (2.76, 2.68, 2.77, 2.80)
        four = (2.47, 2.53, 2.56, 2.66)
        five = (2.40, 2.45, 2.51, 2.60)

        a = (.65, .6, .6, .6)
        b = (4.0, 3.9, 3.2, 3.3)
        c = (7.2, 6.4, 6.6, 6.5)
        d = (5.5, 6.0, 6.1, 6.5)
        e = (3.65, 3.66, 4.1, 4.2)
        f = (2.5, 2.52, 3.12, 3.2)

        zero = (6.84210526316, 7.05882352941, 6.33802816901, 6.15384615385)
        one = (3.54515050167, 3.78184713376, 3.91349124614, 4.02792696026)
        two = (3.25350749879,3.45692779481,3.39275445658,3.16441183215)
        three = (2.76123703419,2.78079084903,2.7889818003,2.77069438996)
        four = (2.5743009126,2.54292263115,2.5817447292,2.65668224745)
        five = (2.3650892892,2.33266656409,2.38028632549,2.39996630444)
        
        a = (13, 12, 9, 8)
        b = (70, 65, 55, 53)
        c = (124, 121, 121, 127)
        d = (105,113,110,113)
        e = (48, 57, 61, 63)
        f = (34, 44, 45, 47)
        
        ind = np.arange(N)    # the x locations for the groups
        width = 0.1       # the width of the bars: can also be len(x) sequence
 
        ax1 = plt.subplot(611)
        plt.bar(ind, zero, width, color='k')
        plt.bar(ind + width, a, width, color='w')
        plt.ylabel('0')
        plt.title('Total SNR w/ iSNR > 2')
        plt.ylim(0, 8) 
        ax1.set_xticklabels([])
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )

        plt.subplot(612, sharex = ax1)
        plt.bar(ind, one, width, color='k')
        plt.bar(ind + width, b, width, color='w')
        plt.ylabel('1')
        plt.ylim(0, 8) 
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )

        plt.subplot(613)
        plt.bar(ind, two, width, color='k')
        plt.bar(ind + width, c, width, color='w')
        plt.ylabel('2')
        plt.ylim(0, 8) 
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )
        
        plt.subplot(614)
        plt.bar(ind, three, width, color='k')
        plt.bar(ind + width, d, width, color='w')
        plt.ylabel('3')
        plt.ylim(0, 8) 
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )
        
        plt.subplot(615)
        plt.bar(ind, four, width, color='k')
        plt.bar(ind + width, e, width, color='w')
        plt.ylabel('4')
        plt.ylim(0, 8) 
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )
        
        plt.subplot(616)
        plt.bar(ind, five, width, color='k')
        plt.bar(ind + width, f, width, color='w')
        plt.ylabel('5')
        plt.ylim(0, 8) 
        plt.xticks(ind+width/2., ('50', '60', '70', '80') )
        plt.xlabel('Percentage of Degradome Expression Within 6nt')
        
        

        plt.show()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
