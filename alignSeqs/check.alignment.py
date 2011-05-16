import bioLibCG

def checkAlignment(alignmentFN, dDB, qFN):

        id_dSeq = {}
        f = open(dDB, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id_dSeq[int(ls[0])] = ls[1]
        f.close()

        id_qSeq = {}
        f = open(qFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id_qSeq[int(ls[0])] = ls[1]
        f.close()


        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                sID = int(ls[0])
                tID = int(ls[1])
                sOffset = int(ls[4])
                try:
                        mismatchPositions = [(int(x) + sOffset) for x in ls[9].split(',')]
                except IndexError:
                        mismatchPositions = []
                ss = id_qSeq[sID]
                ts = id_dSeq[tID]

                sSpaces = ''.join([' ' for i in list(range(0,sOffset))])
                tSpaces = [' ' for i in list(ts)]
                for i in mismatchPositions:
                        tSpaces[i] = 'X'
                

                #graphically display
                print ts
                print ''.join(tSpaces)
                print '%s%s' % (sSpaces, ss)
                print 
                
def checkAlignmentSingle(alignmentFN, ss):


        f = open(alignmentFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                sOffset = int(ls[5])
                ts = ls[20]
                try:
                        mismatchPositions = [(int(x) + sOffset) for x in ls[10].split(',')]
                except IndexError:
                        mismatchPositions = []
                

                sSpaces = ''.join([' ' for i in list(range(0,sOffset))])
                tSpaces = [' ' for i in list(ts)]
                for i in mismatchPositions:
                        tSpaces[i] = 'X'
                

                #graphically display
                print ' '.join(ls[:11])
                print ts
                print ''.join(tSpaces)
                print '%s%s' % (sSpaces, ss)
                print 

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(checkAlignmentSingle, sys.argv)
