import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def checkSims(queryFN, simFN):
    
    id_seq = {}
    f = open(queryFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        id_seq[ls[0]] = ls[1]
    f.close()        

    id_seq2 = {}
    f = open(simFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        id_seq2[ls[0]] = ls[1]

    for id, seq in id_seq.items():
        print id 
        print seq 
        print id_seq2[id]

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
