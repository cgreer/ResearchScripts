import bioLibCG
import cgNexusFlat

def compare(fN, fN2):
        
        set1 = set()
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                set1.add(ls[0])

        
        set2 = set()
        f = open(fN2, 'r')
        for line in f:
                ls = line.strip().split('\t')
                set2.add(ls[0])


        print set1 - set2                

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

        
