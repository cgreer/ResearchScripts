import bioLibCG
import cgNexusFlat

def mergeAnnotations(fN1, fN2, outFN):
        
        print 'getting set info'
        chrom_strand_startSet = {}
        chrom_strand_endSet = {}
        
        f = open(fN1, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                chrom = ls[1]
                strand = ls[2]
                
                starts = ls[8]
                if starts.endswith(','):
                        starts = starts[:-1]

                ends = ls[8]
                if ends.endswith(','):
                        ends = ends[:-1]
                
                starts = [int(x) for x in starts.split(',')]
                ends = [int(x) for x in ends.split(',')]

                
                #add to sets
                for start in starts:
                        chrom_strand_startSet.setdefault(chrom, {}).setdefault(strand, set()).add(start)

                for end in ends:                        
                        chrom_strand_endSet.setdefault(chrom, {}).setdefault(strand, set()).add(end)
        
        f.close()

        newLines = []

        print 'checking new transcripts'
        counter = 0
        f = open(fN2, 'r')
        for line in f:
                ls = line.strip().split('\t')
                
                chrom = ls[1]
                strand = ls[2]
                if chrom not in chrom_strand_startSet: continue
                if chrom not in chrom_strand_endSet: continue

                starts = ls[8]
                if starts.endswith(','):
                        starts = starts[:-1]

                ends = ls[8]
                if ends.endswith(','):
                        ends = ends[:-1]
                
                starts = [int(x) for x in starts.split(',')]
                ends = [int(x) for x in ends.split(',')]

                tranPass = True
                for start in starts:
                        if start in chrom_strand_startSet[chrom][strand]:
                                tranPass = False
                                break

                for end in ends:
                        if end in chrom_strand_endSet[chrom][strand]:
                                tranPass = False
                                break

                if tranPass:
                        newLines.append(line)
                else:
                        counter += 1
        f.close()

        print 'num already present / num adding'
        print counter, len(newLines)
        

        #copy file one, write newLines
        fOut = open(outFN, 'w')
        
        f = open(fN1, 'r')
        for line in f:
                fOut.write(line)
        f.close()

        
        fOut.writelines(newLines)
        fOut.close()

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
