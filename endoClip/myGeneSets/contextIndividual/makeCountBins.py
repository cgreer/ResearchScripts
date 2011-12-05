import bioLibCG
import cgNexusFlat

def makeBins(fN, fOut, typeFilter):
        
        fOut = open(fOut, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                type, tcc = ls[1:3]
                chrom, strand, st, en = bioLibCG.tccSplit(tcc)

                #only take seqs that are long enough
                if en - st < 100: continue
                if not typeFilter in type: continue
                
                tccBins = [] #0 is the first nt from the 3' end
                if strand == '1':
                        for i in range(0, 100):
                                s, e = en - i, en -i
                                tccBins.append(bioLibCG.makeTcc(chrom,strand,s,e)) 

                elif strand == '-1':
                        for i in range(0, 100):
                                s, e = st + i, st + i
                                tccBins.append(bioLibCG.makeTcc(chrom,strand,s,e))


                pString = [id] + tccBins
                fOut.write('\t'.join([str(x) for x in pString]) + '\n')

def makeBins5(fN, fOut, typeFilter):
        
        fOut = open(fOut, 'w')
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                type, tcc = ls[1:3]
                chrom, strand, st, en = bioLibCG.tccSplit(tcc)

                #only take seqs that are long enough
                if en - st < 100: continue
                if not typeFilter in type: continue
                
                tccBins = [] #0 is the first nt from the 3' end
                if strand == '1':
                        for i in range(0, 100):
                                s, e = st + i, st + i
                                tccBins.append(bioLibCG.makeTcc(chrom,strand,s,e)) 

                elif strand == '-1':
                        for i in range(0, 100):
                                s, e = en - i, en - i
                                tccBins.append(bioLibCG.makeTcc(chrom,strand,s,e))


                pString = [id] + tccBins
                fOut.write('\t'.join([str(x) for x in pString]) + '\n')


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
