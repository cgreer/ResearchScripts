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

def add_blank_column(fN, outFN, position = None):
    '''add a blank column of "." BEFORE zero based position.  if position is None than
    add it to last column'''
    position = int(position) if position != None else None

    with open(fN, 'r') as f:
        with open(outFN, 'w') as outF:
            for line in f:
                ls = line.strip().split('\t')
                if not position:
                    ls.append('.')
                else:
                    ls.insert(position, '.')
                outF.write('\t'.join(ls) + '\n')
                    


def copyIDs(fN, outFN):

        fOut = open(outFN, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]

                fOut.write(id + '\n')

        f.close()
        fOut.close()

def addIDs(fN, outFN, header = False, startID = 0):
        '''Add IDs to the first coluimn'''
        if header == 'True':
                header = True
        else:
                header = False
        startID = int(startID)

        newLines = []
        f = open(fN, 'r')
        if header:
                f.readline()
        i = startID
        for line in f:
                ls = line.strip().split('\t')
                newLines.append('%s\t%s\n' % (i, '\t'.join(ls)))
                i += 1
        f.close()

        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()

def addIDsMem(fN, outFN, header = False, startID = 0):
        '''Add IDs to the first coluimn'''
        if header == 'True':
                header = True
        else:
                header = False
        startID = int(startID)

        fO = open(outFN, 'w')
        f = open(fN, 'r')
        if header:
                f.readline()
        i = startID
        for line in f:
                ls = line.strip().split('\t')
                fO.write('%s\t%s\n' % (i, '\t'.join(ls)))
                i += 1
        f.close()
        fO.close()

def redoIDs(fN, outFN, startID = 0):
        '''Re-add IDs to the first column'''
        startID = int(startID)

        
        newLines = []
        f = open(fN, 'r')
        i = startID
        for line in f:
                ls = line.strip().split('\t')
                ls[0] = str(i)
                newLines.append('%s\n' % ('\t'.join(ls)))

                i += 1
        f.close()

        f = open(outFN, 'w')
        f.writelines(newLines)
        f.close()

def redoIDsMem(fN, outFN, startID = 0):
        '''Re-add IDs to the first column'''
        startID = int(startID)

        fO = open(outFN, 'w')        
        f = open(fN, 'r')
        i = startID
        for line in f:
                ls = line.strip().split('\t')
                ls[0] = str(i)
                fO.write('%s\n' % ('\t'.join(ls)))
                i += 1

        f.close()
        fO.close()


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])        
