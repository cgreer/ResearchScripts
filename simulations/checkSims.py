import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast
from cgFile import nextFilePacket

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

def checkSimulationSaturation(simFN, outFN):
    '''simFN is an aggregate of ALL simulations files'''

    with open(simFN, 'r') as f:
        id_seqSet = {} 
        while True:
            nFP = nextFilePacket(f, 3)
            if nFP == None: break
            idS, sequence, blank = nFP 

            id = idS.strip()[1:]
            sequence = sequence.strip()
            id_seqSet.setdefault(id, set()).add(sequence)

    with open(outFN, 'w') as f:  
        for id, seqSet in id_seqSet.items():
            f.write('%s\t%s\n' % (id, len(seqSet))) 

        


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
