import bioLibCG
import cgNexusFlat

def countContext(targetList, alignmentFile, outFN):

        targets = set()
        f = open(targetList, 'r')
        for line in f:
                ls = line.strip().split('\t')
                targets.add(ls[0])
        f.close()

        type_count = {}
        f = open(alignmentFile, 'r')
        for line in f:
                ls = line.strip().split('\t')
                aID = ls[0]
                type = ls[17]
                
                if aID in targets:
                        type_count[type] = type_count.get(type, 0) + 1
        
        f.close()

        f = open(outFN, 'w')
        for type in type_count:
                f.write('%s: %s\n' % (type, type_count[type]))
        f.close()

                                        

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
