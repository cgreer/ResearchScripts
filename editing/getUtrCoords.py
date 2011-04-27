import bioLibCG
import cgGenes3

def getCoords(gFN):

        geneSet = cgGenes3.createGeneSetEditing(gFN)

        for gName in geneSet.set:
                gene = geneSet.set[gName]

                #find the longest utr...
                longestUTR = None
                longestT = None
                longest = 0
                for transcript in gene.transcripts:
                        l = 0
                        if transcript.utr3 is None: continue
                        for utrPair in transcript.utr3:
                                l += utrPair[1] - utrPair[0] + 1
                        if l > longest:
                                longestUTR = transcript.utr3
                                longestT = transcript

                #Get the coordinate ends/etc
                starts, ends = [], []
                if longestUTR is None: continue
                for utrPair in longestUTR:
                        starts.append(utrPair[0])
                        ends.append(utrPair[1])

                starts.sort()
                ends.sort()

                startS = ','.join([str(x) for x in starts])
                endS = ','.join([str(x) for x in ends])

                print '%s\t%s\t%s\t%s\t%s\t%s' % (transcript.id, transcript.parent, transcript.chromosome, transcript.strand, startS, endS)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getCoords, sys.argv)
                
