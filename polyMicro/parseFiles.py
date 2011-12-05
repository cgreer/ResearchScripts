import bioLibCG
#import cgNexusFlat
from cgAutoCast import autocast

def parseMirbaseMature(matureFN, outFN):
    '''puts mir-id, acc, seq'''

    mID, aID, seq = None, None, None 
    fOut = open(outFN, 'w')
    f = open(matureFN, 'r')
    for line in f:
        ls = line.strip().split()
        
        if ls[0].startswith('>'):
            mID = ls[0][1:]
            aID = ls[1]
        else:
            seq = ls[0]

            #write miR 
            fOut.write('%s\t%s\t%s\n' % (mID, aID, seq) )
            
    f.close()
    fOut.close()
    
def getSmallSeqs(smallFN, outFN):
    '''strip sequences from fastq'''

    fOut = open(outFN, 'w')
    f = open(smallFN, 'r')
    for i, line in enumerate(f):

        if (i - 1) % 4 == 0:
            seq = line.strip()
            fOut.write(seq + '\n')

    f.close()
    fOut.close()

def adapterGuess(fqFile, guess):
    '''Displays info if you split by guess adapter'''

    inSize_count = {}
    noAdapter = 0
    f = open(fqFile, 'r')
    for i, line in enumerate(f):

        if (i - 1) % 4 == 0:
            seq = line.strip()
            if guess in seq:
                inSize = len(seq.split(guess)[0])
                inSize_count[inSize] = inSize_count.get(inSize, 0) + 1
            else:
                noAdapter += 1
    f.close()

    print "insize counts"
    for inSize in sorted(inSize_count.keys()):
        print inSize, inSize_count[inSize]
    print "# with no adapter", noAdapter
    print "total # reads", i

@autocast
def removeAdapters(seqFN, adapterSeq, minSize, maxSize, outFN):
    ''' if adapter not in seq --> discard
    if insert size < 18 ---> discard '''

    fOut = open(outFN, 'w')
    f = open(seqFN, 'r')
    for line in f:

        if adapterSeq in line:
            read = line.split(adapterSeq)[0]
            if minSize <= len(read) <= maxSize:
                fOut.write(read + '\n')
        
    f.close()
    fOut.close()

@autocast
def removeAdaptersFQ(seqFN, adapterSeq, minSize, maxSize, outFN):
    ''' if adapter not in seq --> discard
    if insert size < 18 ---> discard '''

    theRead = []
    fOut = open(outFN, 'w')
    f = open(seqFN, 'r')
    for i, line in enumerate(f):

        if i % 4 == 0:
            theRead = [line]

        if (i - 1) % 4 == 0:
            theRead.append(line)

        if (i - 2) % 4 == 0:
            theRead.append(line)
        
        if (i - 3) % 4 == 0:
            if adapterSeq in theRead[1]:
                read = theRead[1].split(adapterSeq)[0]
                readLen = len(read)
                if minSize <= readLen <= maxSize:
                    theRead[1] = read + '\n'
                    theRead.append(line[:readLen] + '\n')
                    fOut.writelines(theRead)
        
    f.close()
    fOut.close()

@autocast
def trimRead(seqFN, trimAmount, outFN):
    ''' this is from the 3' end right now...'''

    theRead = []
    fOut = open(outFN, 'w')
    f = open(seqFN, 'r')
    for i, line in enumerate(f):

        if i % 4 == 0:
            theRead = [line]

        if (i - 1) % 4 == 0:
            theRead.append(line)

        if (i - 2) % 4 == 0:
            theRead.append(line)
        
        if (i - 3) % 4 == 0:
            read = theRead[1].strip()[:-trimAmount]
            readLen = len(read)

            theRead[1] = read + '\n'
            theRead.append(line[:readLen] + '\n')
            fOut.writelines(theRead)
        
    f.close()
    fOut.close()

def convertToU(fN):

    data = open(fN, 'r').read()

    dataU = data.replace('T', 'U')
    fOut = open(fN, 'w')
    fOut.write(dataU)
    fOut.close()

def convertToT(fN):

    data = open(fN, 'r').read()

    dataT = data.replace('U', 'T')
    fOut = open(fN, 'w')
    fOut.write(dataT)
    fOut.close()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
