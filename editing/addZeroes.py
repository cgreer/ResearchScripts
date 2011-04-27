import bioLibCG

def addFiller(fN, filler, zeroPosition, outFN):
        filler = str(filler)
        
        idFlag = False
        if filler == 'ID':
                idFlag = True

        f = open(fN, 'r')
        newLines = []

        i = 0
        for line in f:
                if idFlag:
                        filler = str(i)
                newLines.append(bioLibCG.appendToLine(line, filler, int(zeroPosition)))
                i += 1
        f.close()

        fOut = open(outFN, 'w')
        fOut.writelines(newLines)
        fOut.close()

if __name__ == '__main__':
        import sys
        bioLibCG.submitArgs(addFiller, sys.argv)
                
