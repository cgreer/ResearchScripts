import bioLibCG

def removeDupTargets(fN):
        
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                targets = bioLibCG.uniqueList(ls[4].split(','))
                ls[4] = ','.join(targets)

                newLine = '\t'.join(ls)
                print newLine

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(removeDupTargets, sys.argv)
