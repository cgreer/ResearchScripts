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

        
if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(blankIDs, sys.argv)
