import bioLibCG

def blankIDs(fN, outFN, numIDs = None):
        '''Make a file with X number of blank IDs, or as many as there is lines in a file'''
        
        if numIDs:
                numIDs = int(numIDs)
        else:
                numIDs = bioLibCG.getNumFileLines(fN)

       
        newLines = []
        for i in xrange(0, numIDs):
                newLines.append('%s\n' % i)

        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()
        
def copyIDs(fN, outFN):

        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]

                fOut.write(id + '\n')

        f.close()
        fOut.close()

def addIDs(fN, outFN, header = False):
        '''Add IDs to the first coluimn'''
        if header == 'True':
                header = True

        newLines = []
        f = open(fN, 'r')
        if header:
                f.readline()
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                newLines.append('%s\t%s\n' % (i, '\t'.join(ls)))
                i += 1
        f.close()

        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()
        
def redoIDs(fN, outFN):
        '''Re-add IDs to the first column'''
        
        newLines = []
        f = open(fN, 'r')
        i = 0
        for line in f:
                ls = line.strip().split('\t')
                ls[0] = str(i)
                newLines.append('%s\n' % ('\t'.join(ls)))

                i += 1
        f.close()

        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()



if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])        
