import bioLibCG

def loadSingleWig(wigDir, chrom, strand, prefix):

        coord_value = {}
        fN = wigDir + '/%s.%s.%s.wig' % (prefix, chrom, strand)
        f = open(fN, 'r')
        f.readline()
        for line in f:
                ls = line.strip().split('\t')
                start, end, expr = int(ls[1]) + 1, int(ls[2]), int(ls[3]) #1BASE

                if expr == 0: continue
                
                for i in range(start, end + 1):
                        coord_value[i] = expr
        f.close()

        return coord_value

def loadSingleWigContext(wigDir, chrom, strand, prefix):

        coord_value = {}
        fN = wigDir + '/%s.%s.%s.wig' % (prefix, chrom, strand)
        f = open(fN, 'r')
        f.readline()
        for line in f:
                ls = line.strip().split('\t')
                start, end, expr = int(ls[1]) + 1, int(ls[2]), ls[3] #1BASE

                if expr == 'INTER': continue
                
                for i in range(start, end + 1):
                        coord_value[i] = expr
        f.close()

        return coord_value

def loadSingleWigTranscript(wigDir, chrom, strand, prefix):

        coord_value = {}
        fN = wigDir + '/%s.%s.%s.wig' % (prefix, chrom, strand)
        f = open(fN, 'r')
        f.readline()
        for line in f:
                ls = line.strip().split('\t')
                start, end, expr = int(ls[1]) + 1, int(ls[2]), ls[3] #1BASE

                if expr == 'None': continue
                
                for i in range(start, end + 1):
                        coord_value[i] = expr
        f.close()

        return coord_value

def loadWigDict(wigDir):
        '''Wig files in a directory must be a certain file format: NAME.chr.strand.wig'''

        chr_strand_coord_expr = {}
        
        for fN in bioLibCG.recurseDir(wigDir, end = '.wig'):
                chrom, strand = fN.split('/')[-1].split('.')[1], fN.split('/')[-1].split('.')[2]

                chr_strand_coord_expr.setdefault(chrom, {})[strand] = {}
                f = open(fN, 'r')
                f.readline() #header
                for line in f:
                        ls = line.strip().split('\t')
                        start, end, expr = int(ls[1]) + 1, int(ls[2]), int(ls[3]) #1BASE

                        if expr == 0: continue
                        
                        for i in range(start, end + 1):
                                chr_strand_coord_expr[chrom][strand][i] = expr
                f.close()

        return chr_strand_coord_expr


def getExpressionProfile(tcc, wigDict):
        '''assume 1 based'''

        chrom, strand, start, end = bioLibCG.tccSplit(tcc)
        coord_value = {}

        for i in range(start, end + 1):
                
                if chrom not in wigDict:
                        #print 'WARNING: %s not in wig files' % chrom
                        coord_value[i] = 0 
                        continue

                if strand not in wigDict[chrom]:
                        #print 'WARNING: %s, %s not in wig files' % (chrom, strand)
                        coord_value[i] = 0
                        continue

                coord_value[i] = wigDict[chrom][strand].get(i, 0)

        return coord_value                

                                

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
