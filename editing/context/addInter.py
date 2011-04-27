import bioLibCG

def addInter(tableFN, contextFN, eFN, passedInfo, outFN):

        coord_eID = {}
        eID_Info = {}
        f = open(eFN, 'r')
        f.readline()
        for line in f:
                ls = line.strip().split('\t')
                eID = int(ls[13])
                chrom = ls[0]
                coord = ls[1]
                newcoord = '%s:%s' % (chrom, coord)
                nedited = ls[4]
                ntotal = ls[5]
                coord_eID[newcoord] = eID
                eID_Info[eID] = [newcoord, nedited, ntotal]

        
        eid_gname = {}
        f = open(contextFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                eid = int(ls[0])
                gname = ls[1]
                eid_gname[eid] = gname

        known = {}
        f = open(tableFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                coord = '%s:%s' % (ls[0], ls[1])
                known[coord] = line.strip()

        outF = open(outFN, 'w')
        for coord in coord_eID:
                if coord not in known:
                        eID = coord_eID[coord]
                        outF.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (eID_Info[eID][0].split(':')[0], eID_Info[eID][0].split(':')[1], passedInfo, eid_gname[eID], 'INTER', eID_Info[eID][1], eID_Info[eID][2]))

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(addInter, sys.argv)
