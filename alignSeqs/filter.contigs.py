import bioLibCG


def getContigLength(seq):
        highestLength = 1
        cLength = 1
        letters = list(seq)
        for i,letter in enumerate(letters):
                if i == 0: continue

                if letters[i] == letters[i-1]:
                        cLength += 1
                        if cLength > highestLength:
                                highestLength = cLength
                else:
                        cLength = 1

        return highestLength

def filterContigs(fN, outFN):

        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                seq = ls[1]

                if getContigLength(seq) < 7:
                        fOut.write(line)
                        

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterContigs, sys.argv)


