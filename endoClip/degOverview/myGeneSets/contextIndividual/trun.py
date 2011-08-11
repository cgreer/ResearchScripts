import bioLibCG
import cgNexusFlat

def truncate(fN):
       
       tCount = 0
       shortCount = 0
       fOut = open(fN + '.trun', 'w')
       f = open(fN, 'r')
       for line in f:
               
               ls = line.strip().split('\t')
               type, tcc = ls[1], ls[2]

               tCount += 1
               c, s, st, en = bioLibCG.tccSplit(tcc)

               cLen = en - st

               if cLen < 50:
                        shortCount += 1
                        continue

               if s == '1':
                       en = en - 50
               elif s == '-1':
                       st = st + 50
               else:
                       print 'error'
                       return 1

               ls[2] = bioLibCG.makeTcc(c, s, st, en)
               line = '\t'.join(str(x) for x in [ls[0], ls[1], ls[2]]) + '\n'
               fOut.write(line)

       print shortCount, tCount

def truncate5(fN):
       
       tCount = 0
       shortCount = 0
       fOut = open(fN + '.trun5', 'w')
       f = open(fN, 'r')
       for line in f:
               
               ls = line.strip().split('\t')
               type, tcc = ls[1], ls[2]

               tCount += 1
               c, s, st, en = bioLibCG.tccSplit(tcc)

               cLen = en - st

               if cLen < 50:
                        shortCount += 1
                        continue

               if s == '1':
                       st = st + 50
               elif s == '-1':
                       en = en - 50
               else:
                       print 'error'
                       return 1

               ls[2] = bioLibCG.makeTcc(c, s, st, en)
               line = '\t'.join(str(x) for x in [ls[0], ls[1], ls[2]]) + '\n'
               fOut.write(line)

       print shortCount, tCount

if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

