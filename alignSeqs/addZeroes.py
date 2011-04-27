import bioLibCG

def addFiller(fN, filler, zeroPosition, outFN):
        filler = str(filler)
        f = open(fN, 'r')

        newLines = []
        for line in f:
                newLines.append(bioLibCG.appendToLine(line, filler, int(zeroPosition)))
        f.close()

        fOut = open(outFN, 'w')
        fOut.writelines(newLines)
        fOut.close()

def addIDs(fN, outFN, startFrom = 0):
        
        i = 0
        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                newLS = [str(i)]
                newLS.extend(ls)
                
                fOut.write(' '.join(newLS) + '\n')
                i += 1

        f.close()
        fOut.close()                


if __name__ == '__main__':
        import sys
        #bioLibCG.submitArgs(addFiller, sys.argv)
        bioLibCG.submitArgs(addIDs, sys.argv)
                
