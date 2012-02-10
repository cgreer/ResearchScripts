import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast

#TAG::reverse complement,revComp,sequence parsing::
@autocast
def revComp(seq, rev = True):

    orig_new = dict( (i,j) for i,j in ( ('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C'), ('N', 'N') ) )
    cSeq = ''.join([orig_new[i] for i in seq])

    if rev:
        return cSeq[::-1]
    else:
        print cSeq
        return cSeq


@autocast
def revText(text):

    print text[::-1]

@autocast
def makeFAPipeline(inFile, outFN, smallToggle = True, rc = False):

    print 'small toggle/rev comp', smallToggle, rc

    fOut = open(outFN, 'w')
    f = open(inFile, 'r')
    for line in f:
        ls = line.strip().split('\t')
        if smallToggle:
            theID, theSeq = ls[0], ls[3]
        else:
            theID, theSeq = ls[0], ls[2]

        #rev comp for degradome revcomp format
        if rc:
            theSeq = revComp(theSeq)

        #print fa file
        fOut.write('>%s\n%s\n\n' % (theID, theSeq)) 
    f.close()
    fOut.close()

@autocast
def makeFA(inFile, small = True, rc = False):

    f = open(inFile, 'r')
    for line in f:
        x, smallID, degID, s1, y, s2 = line.strip().split('\t')
        if rc:
            s1 = revComp(s1, True)

        if small:
            print '>%s\n%s\n' % (smallID, s2) 
        else:
            print '>%s\n%s\n' % (degID, s1) 
    f.close()
    

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

