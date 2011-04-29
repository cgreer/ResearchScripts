import bioLibCG

def copyIDs(copyFN, outFN):
        
        f = open(copyFN, 'r')
        newLines = []
        for line in f:
                ls = line.strip().split('\t')
                newLines.append(ls[0] + '\n')
        f.close()

        outF = open(outFN, 'w')
        outF.writelines(newLines)
        outF.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
