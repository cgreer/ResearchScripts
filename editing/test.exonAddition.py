import bioLibCG

def testExon(fN):
        
        f = open(fN, 'r')

        known = []
        for line in f:
                ls = line.strip().split('\t')
                tss, tes, css, ces = int(ls[3]), int(ls[4]), int(ls[5]), int(ls[6])
                tType = ls[13]
                if '_noncoding' in tType:
                        stat5 = True
                        if css == tss or css == tes:
                                stat5 = False

                        stat3 = True
                        if ces == tss or ces == tes:
                                stat3 = False
                                
                        if stat3 and stat5:           
                                if tType not in known:
                                        print tType
                                        known.append(tType)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(testExon, sys.argv)

