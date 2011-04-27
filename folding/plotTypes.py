import bioLibCG
import annotateFold

def plotTypes(fN):


        allTypes = {}
        #fold is on every third line
        f = open(fN, 'r')
        for i,line in enumerate(f):
                if (i - 2) % 3 == 0:
                        try:
                                annoD = annotateFold.annotateFold(line.strip().split(' ')[0])
                                types = annoD[200]
                                #print types
                                if len(types) > 2:
                                        print types
                                        print line.strip().split(' ')[0]

                                if 'stemLoop' in types:
                                        if 'stem' in types:
                                                allTypes['stem'] = allTypes.get('stem', 0) + 1
                                                #print '..added', 'stem'
                                        elif 'loop' in types:
                                                #print '..added', 'loop'
                                                allTypes['loop'] = allTypes.get('loop', 0) + 1
                                        elif 'bulge' in types:
                                                #print '..added', 'bulgeS'
                                                allTypes['bulgeS'] = allTypes.get('bulgeS', 0) + 1
                                else:
                                        if 'ds' in types:
                                                #print '..added', 'ds'
                                                allTypes['ds'] = allTypes.get('ds', 0) + 1
                                        elif 'ss' in types:
                                                #print '..added', 'ss'
                                                allTypes['ss'] = allTypes.get('ss', 0) + 1
                                        elif 'bulge' in types:
                                                #print '..added', 'bulge'
                                                allTypes['bulge'] = allTypes.get('bulge', 0) + 1

                        except:
                                pass
                                #print 'failed'

        for a in allTypes:
                print a, allTypes[a]

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(plotTypes, sys.argv)


                
                
