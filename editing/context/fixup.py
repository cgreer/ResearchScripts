import bioLibCG

def fixup(eFN, tableFN):

        coord_gName = {}
        coord_eRatio = {}
        f = open(eFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                coord = '%s:%s' % (ls[0], ls[1])
                gName = ls[3]
                eRatio = ls[6]
                coord_gName[coord] = gName
                coord_eRatio[coord] = eRatio

        f = open(tableFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                coord = 'chr%s:%s' % (ls[0], ls[1])
                ls.append(coord_eRatio[coord])
                if ls[3] == 'NONE':
                        ls[3] = coord_gName[coord]

                print '\t'.join(ls)                        
                                
        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(fixup, sys.argv)
