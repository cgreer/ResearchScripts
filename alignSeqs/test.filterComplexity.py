import bioLibCG

def filterComp(rFN, compFN):

        id_comp = {}
        f = open(compFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                comp = float(ls[1])
                id_comp[id] = comp

        f = open(rFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                comp = id_comp[id]
                if comp > 1.2:
                        print line.strip()
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterComp, sys.argv)
