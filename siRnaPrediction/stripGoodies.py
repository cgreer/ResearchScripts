import bioLibCG

def stripGoodies(resultsFN, goodList):

        f = open(goodList, 'r')
        goodNums = [int(x.strip()) for x in f]
        f.close()

        f = open(resultsFN, 'r')
        fLines = [x for x in f]
        

        for num in goodNums:
                print '(' + str(num) + ')\t' +  fLines[num - 1].strip()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(stripGoodies, sys.argv)
