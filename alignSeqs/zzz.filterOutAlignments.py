import bioLibCG
import cgDB
import cgAlignment

def filterOut(aFN, oListFN, outFN):
        
        goodIDs = []
        f = open(oListFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                goodIDs.append(line.strip())

        fOut = open(outFN, 'w')
        f = open(aFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                oID = ls[1]

                if oID in goodIDs:
                        fOut.write(line)

                


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterOut, sys.argv)
