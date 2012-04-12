from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as PLT
import numpy as NP
import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

def plotTotalSNR(fN):
    fig = PLT.figure()
    ax1 = fig.add_subplot(111, projection='3d')

    Xs, Ys, dZs = cgDL.listFromColumns(fN, [0,1,2], ['float', 'float', 'float'], naToZero = True)
    #xpos = [1,2,3]
    #ypos = [1,2,3]
    #zpos = [0,0,0]

    Zs = [0] * len(dZs)
    dZs = [1 if x == 0.0 else x for x in dZs]
    dx = dy = [.2] * len(Zs)
    ax1.bar3d(Xs, Ys, Zs, dx, dy, dZs, color='#8E4585', zsort = 'max')
    PLT.show()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])


