import bioLibCG
import cgNexusFlat
from cgAutoCast import autocast

def getCountsForID(fastFN, mappedFN, outFN, maxCount = 51):
        '''chose fasta filename because it will be same
        format...only gawk will be changed'''

        #grab all ids we are tracking...
        allIDs = set()
        f = open(fastFN, 'r')
        for line in f:
                ls = line.strip().split(' ')
                if line.startswith('>'):
                        allIDs.add(ls[1])
        f.close()                        
                
        
        #gather counts that were there
        id_count = {}
        f = open(mappedFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                id_count[id] = id_count.get(id, 0) + 1
        f.close()                

        
        #update "missing" ids to maximum count
        for id in allIDs:
                if id not in id_count:
                        id_count[id] = int(maxCount)
               
        #output id/count               
        f = open(outFN, 'w')
        for id, count in id_count.iteritems():
                f.write('%s\t%s\n' % (id, count))


if __name__ == "__main__":
        import sys
        if sys.argv[1] == "help":
                bioLibCG.gd(sys.argv[0])
        else:
                bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
