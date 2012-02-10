import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast

@autocast
def grepBC(grepList, inFile, column, word = False):

    grepList = cgDL.listFromColumns(grepList, [column], ['int'])
   
    f = open(inFile, 'r')
    for line in f:
        ls = line.strip().split('\t')
        
        for w in grepList:
            if word:
                if w == ls[column]:
                    print line,
            else:
                if w in ls[column]:
                    print line,
    f.close()
        

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

