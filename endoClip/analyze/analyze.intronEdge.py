import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

    

@autocast
def plotInfoNTFrameEnrichment(dir, frameWidth, outFN):

    frameNum_seqs = {}
    for fChrom in bioLibCG.humanChromosomes:
        for fStrand in ('1', '-1'):
            print fChrom, fStrand
            fN = '%s/%s.%s.primeSeqs' % (dir, fChrom, fStrand)
            f = open(fN, 'r')
            for line in f:
                ls = line.strip().split('\t')
                seq = ls[0]
                for i, frame in enumerate(bioLibCG.returnFrames(seq, frameWidth)):
                    frameNum_seqs.setdefault(i, []).append(frame)
            f.close()


    let_frameCounts = {}
    for fNum in sorted(frameNum_seqs.keys()):
        seqs = frameNum_seqs[fNum]
        let_count = getSeqEnrichment(seqs)
        for let, count in let_count.items():
            let_frameCounts.setdefault(let, []).append(count)


    fOut = open(outFN, 'w')
    for let, fCounts in let_frameCounts.items():
        fCounts = ','.join([str(x) for x in fCounts])
        fOut.write('%s\t%s\n' % (let, fCounts))
    fOut.close()

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
