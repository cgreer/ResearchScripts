import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def takeRepeatsAway(aFN, dRNAF, newAFN):
    '''Given dRNA file without repeats, filter alignments
    Easier than filtering aligning db'''
    
    rIDs = set()
    f = open(dRNAF, 'r')
    for line in f:
        ls = line.strip().split('\t')
        rIDs.add(ls[0])
    f.close()        


    outF = open(newAFN, 'w')
    f = open(aFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        if ls[2] not in rIDs: continue

        outF.write(line)            
    f.close()
    outF.close()




if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
